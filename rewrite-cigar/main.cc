/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>	// std::make_heap etc.
#include <chrono>
#include <cstdlib>		// std::exit
#include <filesystem>
#include <iostream>
#include <libbio/assert.hh>
#include <map>
#include <mutex>
#include <panvc3/compressed_fasta_reader.hh>
#include <panvc3/utility.hh>
#include <range/v3/view/enumerate.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <syncstream>
#include <thread>
#include <vector>
#include "cmdline.h"

namespace chrono	= std::chrono;
namespace fs		= std::filesystem;
namespace lb		= libbio;
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
		node.mapped().index = seq_idx;
		auto const res(m_sequences.insert(std::move(node)));
		libbio_assert(res.inserted);
		m_sequences_by_use.front() = res.position;
		std::push_heap(m_sequences_by_use.begin(), m_sequences_by_use.end(), sequence_buffer_cmp{});
		return res.position->second.sequence;
	}
	
	
	template <typename t_header>
	void append_program_info(t_header const &input_header, t_header &output_header, int const argc, char const * const * const argv)
	{
		output_header.program_infos = input_header.program_infos;
		panvc3::append_sam_program_info(
			"panvc3.rewrite-cigar.",
			"PanVC 3 rewrite_cigar",
			argc,
			argv,
			CMDLINE_PARSER_VERSION,
			output_header.program_infos
		);
	}
	
	
	template <typename t_query_seq>
	void rewrite_cigar(
		std::vector <char> const &ref,
		std::size_t ref_pos,
		t_query_seq const &query_seq,
		std::vector <seqan3::cigar> const &cigar_seq,
		std::vector <seqan3::cigar> &cigar_output
	)
	{
		using seqan3::get;
		using seqan3::operator""_cigar_operation;
		
		std::size_t query_pos{};
		for (auto const &item : cigar_seq)
		{
			auto const count(get <0>(item));
			auto const op(get <1>(item));
			auto const op_(op.to_char());
			
			switch (op_)
			{
				case 'H':	// Hard clipping, consumes nothing.
				case 'P':	// Padding (silent deletion from padded reference), consumes nothing.
					cigar_output.push_back(item);
					break;
					
				case 'I':	// Insertion, consumes query.
				case 'S':	// Soft clipping, consumes query.
					cigar_output.push_back(item);
					query_pos += count;
					break;
				
				case '=':	// Match, consumes both.
				case 'X':	// Mismatch, consumes both.
					cigar_output.push_back(item);
					query_pos += count;
					ref_pos += count;
					break;
				
				case 'D':	// Deletion, consumes reference.
				case 'N':	// Skipped region, consumes reference. (In SAMv1, this is only relevant in mRNA-to-genome alignments.)
					cigar_output.push_back(item);
					ref_pos += count;
					break;
				
				case 'M':	// Match or mismatch, consumes both.
				{
					libbio_assert_lte(query_pos + count, query_seq.size());
					std::size_t prev_count(1);
					auto prev_op(query_seq[query_pos].to_char() == ref[ref_pos] ? '='_cigar_operation : 'X'_cigar_operation);
					for (std::size_t i(1); i < count; ++i)
					{
						auto const curr_op(query_seq[query_pos + i].to_char() == ref[ref_pos + i] ? '='_cigar_operation : 'X'_cigar_operation);
						if (curr_op == prev_op)
							++prev_count;
						else
						{
							auto &new_item(cigar_output.emplace_back());
							new_item = prev_count; // Assigns to one element.
							new_item = prev_op;
							
							prev_count = 1;
							prev_op = curr_op;
						}
					}
					
					auto &new_item(cigar_output.emplace_back());
					new_item = prev_count;
					new_item = prev_op;
				}
				
				default:
					libbio_fail("Unexpected CIGAR operation “", op_, "”");
					break;
			}
		}
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	void process_alignments(t_aln_input &&aln_input, t_aln_output &&aln_output, char const *ref_path, std::uint16_t const status_output_interval)
	{
		typedef chrono::steady_clock	clock_type;
		
		// Load the reference sequences.
		lb::log_time(std::cerr) << "Loading the reference sequences…\n";
		panvc3::compressed_fasta_reader fasta_reader;
		fasta_reader.open_path(ref_path);
		
		sequence_buffer_store reference_buffer_store;
		auto const &ref_names(aln_input.header().ref_ids());
		
		lb::log_time(std::cerr) << "Processing the alignments…\n";
		std::uint64_t rec_idx{};
		auto const start_time(clock_type::now());
		panvc3::timer status_output_timer;
		auto status_output_thread{[&]() -> std::jthread {
			if (!status_output_interval)
				return {};
			
			return std::jthread{[&](){
				while (status_output_timer.wait_for(chrono::minutes(status_output_interval)))
				{
					auto const pp(clock_type::now());
					auto const running_time{pp - start_time};

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
			}};
		}()};
		
		std::vector <seqan3::cigar> cigar_buffer;
		for (auto &aln_rec : aln_input)
		{
			++rec_idx;
			aln_rec.header_ptr() = &aln_output.header();
			
			auto const pos_(aln_rec.reference_position());
			if (!pos_.has_value())
			{
				aln_output.push_back(aln_rec);
				continue;
			}
			auto const pos(*pos_);
			
			auto const &ref_id_(aln_rec.reference_id());
			if (!ref_id_.has_value())
			{
				aln_output.push_back(aln_rec);
				continue;
			}
			auto const ref_id(*ref_id_);
			
			// Load the sequence if needed.
			auto &buffer(reference_buffer_store.get_buffer(ref_id, rec_idx));
			if (buffer.empty())
			{
				auto const &ref_name(ref_names[ref_id]);
				lb::log_time(std::osyncstream(std::cerr)) << "(Re-)loading reference sequence '" << ref_name << "'…\n" << std::flush;
				if (!fasta_reader.read_sequence(ref_name, buffer))
				{
					std::osyncstream(std::cerr) << "ERROR: Unable to load sequence ‘" << ref_name << "’ from the input FASTA.\n" << std::flush;
					std::abort();
				}
				lb::log_time(std::osyncstream(std::cerr)) << "Loading complete.\n" << std::flush;
			}
			
			rewrite_cigar(buffer, pos, aln_rec.sequence(), aln_rec.cigar_sequence(), cigar_buffer);
			
			{
				using std::swap;
				swap(aln_rec.cigar_sequence(), cigar_buffer);
			}
			
			aln_output.push_back(aln_rec);
		}
		
		status_output_timer.stop();
		if (status_output_thread.joinable())
			status_output_thread.join();
	}
	
	
	void process(gengetopt_args_info const &args_info, int const argc, char * const * const argv)
	{
		// Sanity check.
		if (args_info.alignments_arg && args_info.bam_input_flag)
		{
			std::cerr << "ERROR: --bam-input has no effect when reading input from a file.\n";
			std::exit(EXIT_FAILURE);
		}
		
		// Open the SAM input.
		typedef seqan3::sam_file_input <> input_type;
		auto aln_input{[&](){
			auto const make_input_type([&]<typename t_fmt>(t_fmt &&fmt){
				if (args_info.alignments_arg)
				{
					fs::path const path(args_info.alignments_arg);
					return input_type(path);
				}
				else
				{
					return input_type(std::cin, std::forward <t_fmt>(fmt));
				}
			});

			if (args_info.bam_input_flag)
				return make_input_type(seqan3::format_bam{});
			else
				return make_input_type(seqan3::format_sam{});
		}()};
		auto const &input_ref_ids(aln_input.header().ref_ids()); // ref_ids() not const.
		
		typedef std::remove_cvref_t <decltype(input_ref_ids)>			ref_ids_type;
		typedef seqan3::sam_file_output <
			typename input_type::selected_field_ids,
			seqan3::type_list <seqan3::format_sam, seqan3::format_bam>,
			ref_ids_type
		>																output_type;
		
		// Status output interval
		if (args_info.status_output_interval_arg < 0)
		{
			std::cerr << "ERROR: Status output interval must be non-negative.\n";
			std::exit(EXIT_FAILURE);
		}
		
		// Open the alignment output file.
		auto aln_output{[&](){
			// Make sure that aln_output has some header information.
			auto const make_output_type([&]<typename t_fmt>(t_fmt &&fmt){
				ref_ids_type empty_ref_ids;
				if (args_info.output_path_arg)
				{
					return output_type(
						fs::path{args_info.output_path_arg},
						std::move(empty_ref_ids),
						rsv::empty <std::size_t>()	// Reference lengths; the constructor expects a forward range.
					);
				}
				else
				{
					return output_type(
						std::cout,
						std::move(empty_ref_ids),
						rsv::empty <std::size_t>(),	// Reference lengths; the constructor expects a forward range.
						std::forward <t_fmt>(fmt)
					);
				}
			});

			if (args_info.output_bam_flag)
				return make_output_type(seqan3::format_bam{});
			else
				return make_output_type(seqan3::format_sam{});
		}()};
		
		// Prepare the output header.
		aln_output.header().ref_ids() = input_ref_ids;
		aln_output.header().ref_id_info = aln_input.header().ref_id_info;
		aln_output.header().ref_dict = aln_input.header().ref_dict;
		append_program_info(aln_input.header(), aln_output.header(), argc, argv);
		
		process_alignments(aln_input, aln_output, args_info.reference_arg, args_info.status_output_interval_arg);
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
		std::cerr << "PID: " << getpid() << '\n';
	
	process(args_info, argc, argv);
	
	return 0;
}
