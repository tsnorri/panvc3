/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>	// std::make_heap etc.
#include <chrono>
#include <cstdlib>		// std::exit
#include <filesystem>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/sam.hh>
#include <map>
#include <mutex>
#include <panvc3/alignment_input.hh>
#include <panvc3/cigar.hh>
#include <panvc3/compressed_fasta_reader.hh>
#include <panvc3/sam_tag.hh>
#include <panvc3/utility.hh>
#include <range/v3/view/enumerate.hpp>
#include <regex>
#include <syncstream>
#include <thread>
#include <vector>
#include "cmdline.h"

namespace chrono	= std::chrono;
namespace fs		= std::filesystem;
namespace lb		= libbio;
namespace sam		= libbio::sam;
namespace rsv		= ranges::views;


namespace {
	
	// For use with single thread.
	struct sequence_buffer
	{
		typedef std::vector <char>	buffer_type;
		typedef std::uint64_t		sequence_index_type;
	
		buffer_type			sequence;
		sequence_index_type	index{};
		
		sequence_buffer(sequence_index_type const index_):
			index(index_)
		{
		}
	};
	
	
	class sequence_buffer_store
	{
	public:
		typedef sequence_buffer::sequence_index_type	sequence_index_type;
		typedef sequence_buffer::buffer_type			buffer_type;
		constexpr static inline size_t const MAX_SIZE{4};
	
	private:
		struct sequence_buffer_cmp
		{
			template <typename t_iterator>
			bool operator()(t_iterator &&lhs, t_iterator &&rhs) const { return lhs->second.index > rhs->second.index; } // Minimum first.
		};
		
		typedef std::map <std::size_t, sequence_buffer>	sequence_map;
	
	private:
		sequence_map							m_sequences;
		std::vector <sequence_map::iterator>	m_sequences_by_use;
		
	public:
		buffer_type &get_buffer(std::size_t const idx, sequence_index_type const seq_idx);
	};
	
	
	auto sequence_buffer_store::get_buffer(std::size_t const idx, sequence_index_type const seq_idx) -> buffer_type &
	{
		{
			auto const it(m_sequences.find(idx));
			if (m_sequences.end() != it)
			{
				it->second.index = seq_idx;
				std::make_heap(m_sequences_by_use.begin(), m_sequences_by_use.end(), sequence_buffer_cmp{});
				return it->second.sequence;
			}
		}
		
		if (m_sequences.size() < MAX_SIZE)
		{
			auto const pp(m_sequences.emplace(idx, seq_idx));
			libbio_assert(pp.second);
			m_sequences_by_use.emplace_back(pp.first);
			std::push_heap(m_sequences_by_use.begin(), m_sequences_by_use.end(), sequence_buffer_cmp{});
			return pp.first->second.sequence;
		}
		
		// Re-use a buffer.
		auto it(m_sequences_by_use.front());
		auto node(m_sequences.extract(it));
		node.key() = idx;
		node.mapped().sequence.clear();
		node.mapped().index = seq_idx;
		auto const res(m_sequences.insert(std::move(node)));
		libbio_assert(res.inserted);
		m_sequences_by_use.front() = res.position;
		std::push_heap(m_sequences_by_use.begin(), m_sequences_by_use.end(), sequence_buffer_cmp{});
		return res.position->second.sequence;
	}
	
	
	void append_program_info(sam::header &output_header, int const argc, char const * const * const argv)
	{
		panvc3::append_sam_program_info(
			"panvc3.rewrite-cigar.",
			"PanVC 3 rewrite_cigar",
			argc,
			argv,
			CMDLINE_PARSER_VERSION,
			output_header.programs
		);
	}


	inline void fill_ref_n_positions(
		std::vector <char> const &ref,
		std::vector <std::uint8_t> &ref_n_positions,
		std::size_t const ref_base_pos,
		std::size_t const ref_pos,
		std::size_t const count
	)
	{
		// Reserve space.
		ref_n_positions.resize((ref_pos + count + 7) / 8, 0);

		// Mark the N character positions.
		for (std::size_t i(0); i < count; ++i)
		{
			auto const cc(ref[ref_base_pos + ref_pos + i]);
			if ('N' == cc)
			{
				auto const bidx((ref_pos + i) / 8);
				auto const mask(std::uint8_t(1) << ((ref_pos + i) % 8));
				ref_n_positions[bidx] |= mask;
			}
		}
	}
	
	
	void rewrite_cigar_alignment_match(
		std::vector <sam::cigar_run> &cigar_seq, // non-const for swapping.
		panvc3::cigar_buffer_libbio &cigar_output
	)
	{
		using sam::operator""_cigar_operation;
		
		cigar_output.clear();
		
		for (auto const &cigar_run : cigar_seq)
		{
			auto const count(cigar_run.count());
			auto const op(cigar_run.operation());
			auto const op_(sam::to_char(op));
			
			switch (op_)
			{
				case 'H':	// Hard clipping, consumes nothing.
				case 'P':	// Padding (silent deletion from padded reference), consumes nothing.
				case 'I':	// Insertion, consumes query.
				case 'S':	// Soft clipping, consumes query.
				case 'D':	// Deletion, consumes reference.
				case 'N':	// Skipped region, consumes reference. (In SAMv1, this is only relevant in mRNA-to-genome alignments.)
				case 'M':	// Match or mismatch, consumes both.
					cigar_output.push_back(op, count);
					break;
				
				case '=':	// Match, consumes both.
				case 'X':	// Mismatch, consumes both.
					cigar_output.push_back('M'_cigar_operation, count);
					break;
				
				default:
					libbio_fail("Unexpected CIGAR operation “", op_, "”");
					break;
			}
		}
		
		cigar_output.finish();
		cigar_output.swap_buffer(cigar_seq);
	}
	
	
	template <typename t_query_seq>
	void rewrite_cigar_sequence_match(
		std::vector <char> const &ref,
		std::size_t const ref_base_pos,
		t_query_seq const &query_seq,
		std::vector <sam::cigar_run> const &cigar_seq,
		std::vector <sam::cigar_run> &cigar_output,
		std::vector <std::uint8_t> &ref_n_positions_output
	)
	{
		using sam::operator""_cigar_operation;

		cigar_output.clear();
		ref_n_positions_output.clear();
		
		std::size_t query_pos{};
		std::size_t ref_pos{};
		for (auto const &cigar_run : cigar_seq)
		{
			auto const count(cigar_run.count());
			auto const op(cigar_run.operation());
			auto const op_(sam::to_char(op));
			
			switch (op_)
			{
				case 'H':	// Hard clipping, consumes nothing.
				case 'P':	// Padding (silent deletion from padded reference), consumes nothing.
					cigar_output.push_back(cigar_run);
					break;
					
				case 'I':	// Insertion, consumes query.
				case 'S':	// Soft clipping, consumes query.
					cigar_output.push_back(cigar_run);
					query_pos += count;
					break;
				
				case '=':	// Match, consumes both.
				case 'X':	// Mismatch, consumes both.
					cigar_output.push_back(cigar_run);
					fill_ref_n_positions(ref, ref_n_positions_output, ref_base_pos, ref_pos, count);
					query_pos += count;
					ref_pos += count;
					break;
				
				case 'D':	// Deletion, consumes reference.
				case 'N':	// Skipped region, consumes reference. (In SAMv1, this is only relevant in mRNA-to-genome alignments.)
					cigar_output.push_back(cigar_run);
					fill_ref_n_positions(ref, ref_n_positions_output, ref_base_pos, ref_pos, count);
					ref_pos += count;
					break;
				
				case 'M':	// Match or mismatch, consumes both.
				{
					libbio_assert_lte(query_pos + count, query_seq.size());
					fill_ref_n_positions(ref, ref_n_positions_output, ref_base_pos, ref_pos, count);

					std::size_t prev_count(1);
					auto prev_op(query_seq[query_pos] == ref[ref_base_pos + ref_pos] ? '='_cigar_operation : 'X'_cigar_operation);
					for (std::size_t i(1); i < count; ++i)
					{
						auto const curr_op(query_seq[query_pos + i] == ref[ref_base_pos + ref_pos + i] ? '='_cigar_operation : 'X'_cigar_operation);
						if (curr_op == prev_op)
							++prev_count;
						else
						{
							auto &new_item(cigar_output.emplace_back(prev_op, prev_count));
							prev_count = 1;
							prev_op = curr_op;
						}
					}
					
					auto &new_item(cigar_output.emplace_back(prev_op, prev_count));
					ref_pos += count;
					query_pos += count;
					break;
				}
				
				default:
					libbio_fail("Unexpected CIGAR operation “", op_, "”");
					break;
			}
		}
	}
	
	
	std::jthread start_status_output_thread(
		panvc3::timer &status_output_timer,
		std::uint64_t &rec_idx,
		std::uint16_t const status_output_interval
	)
	{
		typedef chrono::steady_clock	clock_type;
		
		if (!status_output_interval)
			return {};
		
		return std::jthread{
			[
				&status_output_timer,
				&rec_idx,
				status_output_interval
			](){
				auto const start_time(clock_type::now());
				while (status_output_timer.wait_for(chrono::minutes(status_output_interval)))
				{
					auto const pp(clock_type::now());
					auto const running_time{pp - start_time};
					
					// Inaccurate b.c. rec_idx is not atomic.
					// FIXME: come up with a better way to report status.
					std::osyncstream cerr(std::cerr);
					lb::log_time(cerr) << "Time spent processing: ";
					panvc3::log_duration(cerr, running_time);
					cerr << "; processed " << (rec_idx - 1) << " records";
					
					if (rec_idx)
					{
						double usecs_per_record(chrono::duration_cast <chrono::microseconds>(running_time).count());
						usecs_per_record /= rec_idx;
						cerr << " (in " << usecs_per_record << " µs / record)";
					}
					
					cerr << ".\n" << std::flush;
				}
			}
		};
	}
	
	
	void process_alignments_alignment_match(
		panvc3::alignment_input &aln_input,
		std::ostream &os,
		std::uint16_t const status_output_interval
	)
	{
		lb::log_time(std::cerr) << "Processing the alignments…\n";
		std::uint64_t rec_idx{};
		panvc3::timer status_output_timer;
		auto status_output_thread(start_status_output_thread(status_output_timer, rec_idx, status_output_interval));
		panvc3::cigar_buffer_libbio cigar_buffer;
		
		aln_input.read_records(
			[
				&aln_input,
				&os,
				&rec_idx,
				&cigar_buffer
			](sam::record &aln_rec){
				++rec_idx;
				rewrite_cigar_alignment_match(aln_rec.cigar, cigar_buffer);
				sam::output_record(os, aln_input.header, aln_rec);
			}
		);
		
		status_output_timer.stop();
		if (status_output_thread.joinable())
			status_output_thread.join();
	}
	
	
	void process_alignments_sequence_match(
		panvc3::alignment_input &aln_input,
		std::ostream &os,
		char const *ref_path,
		sam::tag_type const ref_n_positions_tag,
		std::uint16_t const status_output_interval
	)
	{
		typedef std::vector <std::uint8_t> n_positions_buffer_type;
		
		// Load the reference sequences.
		lb::log_time(std::cerr) << "Loading the reference sequences…\n";
		panvc3::compressed_fasta_reader fasta_reader;
		fasta_reader.open_path(ref_path);
		
		sequence_buffer_store reference_buffer_store;
		auto const &ref_entries(aln_input.header.reference_sequences);
		
		lb::log_time(std::cerr) << "Processing the alignments…\n";
		std::uint64_t rec_idx{};
		panvc3::timer status_output_timer;
		auto status_output_thread(start_status_output_thread(status_output_timer, rec_idx, status_output_interval));
		std::vector <sam::cigar_run> cigar_buffer;
		n_positions_buffer_type n_positions_buffer;
		
		aln_input.read_records(
			[
				&aln_input,
				&os,
				&fasta_reader,
				&reference_buffer_store,
				&ref_entries,
				&rec_idx,
				&cigar_buffer,
				&n_positions_buffer,
				ref_n_positions_tag
			](sam::record &aln_rec){
				++rec_idx;
				
				auto const pos(aln_rec.pos);
				if (sam::INVALID_POSITION == pos)
				{
					sam::output_record(os, aln_input.header, aln_rec);
					return;
				}
				
				auto const ref_id(aln_rec.rname_id);
				if (sam::INVALID_REFERENCE_ID == ref_id)
				{
					sam::output_record(os, aln_input.header, aln_rec);
					return;
				}
				
				// Load the sequence if needed.
				auto &buffer(reference_buffer_store.get_buffer(ref_id, rec_idx));
				if (buffer.empty())
				{
					auto const &ref_entry(ref_entries[ref_id]);
					auto const &ref_name(ref_entry.name);
					lb::log_time(std::osyncstream(std::cerr)) << "(Re-)loading reference sequence ‘’" << ref_name << "’…\n" << std::flush;
					if (!fasta_reader.read_sequence(ref_name, buffer))
					{
						std::osyncstream(std::cerr) << "ERROR: Unable to load sequence ‘" << ref_name << "’ from the input FASTA.\n" << std::flush;
						std::abort();
					}
					lb::log_time(std::osyncstream(std::cerr)) << "Loading complete.\n" << std::flush;
				}
				
				rewrite_cigar_sequence_match(buffer, pos, aln_rec.seq, aln_rec.cigar, cigar_buffer, n_positions_buffer);
				// SAM reference does not state that empty arrays are not allowed (I think) but Samtools reports a parse error anyway.
				if (!n_positions_buffer.empty())
					aln_rec.optional_fields.obtain <n_positions_buffer_type>(ref_n_positions_tag) = n_positions_buffer;
			
				{
					using std::swap;
					swap(aln_rec.cigar, cigar_buffer);
				}
			
				sam::output_record(os, aln_input.header, aln_rec);
			}
		);
		
		status_output_timer.stop();
		if (status_output_thread.joinable())
			status_output_thread.join();
	}
	
	
	void process(gengetopt_args_info const &args_info, int const argc, char * const * const argv)
	{
		// Open the SAM input.
		auto aln_input(panvc3::alignment_input::open_path_or_stdin(args_info.alignments_arg));
		
		// Tags.
		std::regex const tag_regex{"^[XYZ][A-Za-z0-9]$"};
		auto const make_sam_tag([&tag_regex](char const *tag) -> sam::tag_type {
			if (!tag)
				return 0;

			if (!std::regex_match(tag, tag_regex))
			{
				std::cerr << "ERROR: The given tag '" << tag << "' does not match the expected format.\n";
				std::exit(EXIT_FAILURE);
			}

			std::array <char, 2> buffer{tag[0], tag[1]};
			return sam::to_tag(buffer);
		});
		auto const ref_n_positions_tag(make_sam_tag(args_info.ref_n_positions_tag_arg));
		
		// Status output interval
		if (args_info.status_output_interval_arg < 0)
		{
			std::cerr << "ERROR: Status output interval must be non-negative.\n";
			std::exit(EXIT_FAILURE);
		}
		
		auto aln_output_fh{[&] -> lb::file_handle {
			if (args_info.output_path_arg)
				return lb::file_handle(lb::open_file_for_writing(args_info.output_path_arg, lb::writing_open_mode::CREATE));
			else
				return lb::file_handle(STDOUT_FILENO, false);
		}()};
		
		lb::file_ostream os;
		lb::open_stream_with_file_handle(os, aln_output_fh);
		
		// Prepare the output header.
		auto &header(aln_input.header);
		auto const &input_ref_ids(header.reference_sequences);
		append_program_info(header, argc, argv); // Adding our program info does not affect parsing.
		
		if (args_info.output_sequence_match_ops_given)
			process_alignments_sequence_match(aln_input, os, args_info.reference_arg, ref_n_positions_tag, args_info.status_output_interval_arg);
		else if (args_info.output_alignment_match_ops_given)
			process_alignments_alignment_match(aln_input, os, args_info.status_output_interval_arg);
		else
			libbio_fail("Unexpected mode");
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	//lb::setup_allocated_memory_logging();

	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.

	if (args_info.print_invocation_given)
	{
		std::cerr << "Invocation:";
		for (int i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << '\n';
	}

	if (args_info.print_pid_given)
		std::cerr << "PID: " << ::getpid() << '\n';
	
	process(args_info, argc, argv);
	
	return 0;
}
