/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <iostream>
#include <libbio/dispatch/dispatch_caller.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/file_handling.hh>
#include <panvc3/alignment_projector.hh>
#include <panvc3/spsc_queue.hh>
#include <panvc3/utility.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/iterator/insert_iterators.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/take_exactly.hpp>
#include <range/v3/view/transform.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include "cmdline.h"

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace rsv	= ranges::views;


namespace {
	
	constexpr inline std::size_t	QUEUE_SIZE{512};
	constexpr inline std::size_t	CHUNK_SIZE{32};
	
	
	template <typename t_string>
	panvc3::msa_index::chr_entry_vector::const_iterator
	find_chr_entry_(panvc3::msa_index::chr_entry_vector const &entries, t_string const &chr_id)
	{
		panvc3::msa_index::chr_entry_cmp chr_cmp;
		auto const it(std::lower_bound(entries.begin(), entries.end(), chr_id, chr_cmp));
		if (entries.end() == it || it->chr_id != chr_id)
		{
			std::cerr << "ERROR: Did not find an entry for chromosome ID “" << chr_id << "” in the MSA index." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		return it;
	}
	
	
	template <typename t_string>
	panvc3::msa_index::chr_entry const &
	find_chr_entry(panvc3::msa_index::chr_entry_vector const &entries, t_string const &chr_id)
	{
		return *find_chr_entry_(entries, chr_id);
	}
	
	
	template <typename t_string>
	panvc3::msa_index::sequence_entry_vector::const_iterator
	find_sequence_entry_(panvc3::msa_index::sequence_entry_vector const &entries, t_string const &seq_id)
	{
		panvc3::msa_index::sequence_entry_cmp seq_cmp;
		auto const it(std::lower_bound(entries.begin(), entries.end(), seq_id, seq_cmp));
		if (entries.end() == it || it->seq_id != seq_id)
		{
			std::cerr << "ERROR: Did not find an entry for sequence ID “" << seq_id << "” in the MSA index." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		return it;
	}
	
	
	template <typename t_string>
	panvc3::msa_index::sequence_entry const &
	find_sequence_entry(panvc3::msa_index::sequence_entry_vector const &entries, t_string const &seq_id)
	{
		return *find_sequence_entry_(entries, seq_id);
	}
	
	
	struct alignment_statistics
	{
		std::size_t flags_not_matched{};
		std::size_t ref_id_missing{};
		std::size_t matched_reads{};
	};
	
	
	template <typename t_input_processor>
	class project_task
	{
		friend t_input_processor;
		
	public:
		typedef typename t_input_processor::input_record_type	record_type;
		typedef std::array <record_type, CHUNK_SIZE>			array_type;
		
	protected:
		t_input_processor	*m_input_processor{};
		array_type			m_records;
		std::size_t			m_valid_records{};
		
	public:
		bool is_full() const { return CHUNK_SIZE == m_valid_records; }
		bool empty() const { return 0 == m_valid_records; }
		record_type &next_record() { return m_records[m_valid_records++]; }
		auto alignment_records() { return m_records | rsv::take_exactly(m_valid_records); }
		auto alignment_records() const { return m_records | rsv::take_exactly(m_valid_records); }
		
		void process();
		void output();
		void reset() { m_valid_records = 0; }
	};
	
	
	template <typename t_aln_input, typename t_aln_output>
	class input_processor
	{
		static_assert(panvc3::is_power_of_2(QUEUE_SIZE));
		
	public:
		typedef t_aln_input											input_type;
		typedef t_aln_output										output_type;
		typedef std::remove_cvref_t <
			decltype(
				std::declval <t_aln_output>().header().ref_ids()
			)
		>															output_reference_ids_type;
		typedef typename input_type::record_type					input_record_type;
		typedef project_task <input_processor>						project_task_type;
		typedef std::vector <char>									sequence_vector;
		typedef panvc3::spsc_queue <project_task_type, QUEUE_SIZE>	queue_type;
		typedef lb::dispatch_ptr <dispatch_queue_t>					dispatch_queue_ptr;
		typedef lb::dispatch_ptr <dispatch_semaphore_t>				semaphore_ptr;
		
	protected:
		panvc3::msa_index			m_msa_index;
		output_reference_ids_type	m_output_reference_ids;	// Has to be before m_aln_output since the latter can have a reference to it.
		input_type					m_aln_input;
		output_type					m_aln_output;
		sequence_vector				m_ref_sequence;
		
		dispatch_queue_ptr			m_output_dispatch_queue;
		
		queue_type					m_task_queue{};
		
		alignment_statistics		m_statistics;
		std::string					m_ref_id;
		std::string					m_msa_ref_id;
		std::string					m_output_seq_id;
		std::string					m_ref_id_separator;
		std::int32_t				m_gap_opening_cost{};
		std::int32_t				m_gap_extension_cost{};
		bool						m_should_consider_primary_alignments_only{};
		bool						m_should_use_read_base_qualities{};
		
	public:
		template <
			typename t_ref_id,
			typename t_msa_ref_id,
			typename t_output_seq_id,
			typename t_ref_id_separator
		>
		input_processor(
			panvc3::msa_index			&&msa_index,
			t_aln_input					&&aln_input,
			t_aln_output				&&aln_output,
			sequence_vector				&&ref_sequence,
			output_reference_ids_type	&&output_reference_ids,
			t_ref_id 					&&ref_id,
			t_msa_ref_id				&&msa_ref_id,
			t_output_seq_id				&&output_seq_id,
			t_ref_id_separator			&&ref_id_separator,
			std::int32_t				gap_opening_cost,
			std::int32_t				gap_extension_cost,
			bool						should_consider_primary_alignments_only,
			bool						should_use_read_base_qualities
		):
			m_msa_index(std::move(msa_index)),
			m_output_reference_ids(std::move(output_reference_ids)),
			m_aln_input(std::move(aln_input)),
			m_aln_output(std::move(aln_output)),
			m_ref_sequence(std::move(ref_sequence)),
			m_output_dispatch_queue(dispatch_queue_create("fi.iki.tsnorri.panvc3.project-alignments.output-queue", DISPATCH_QUEUE_SERIAL)),
			m_ref_id(std::forward <t_ref_id>(ref_id)),
			m_msa_ref_id(std::forward <t_msa_ref_id>(msa_ref_id)),
			m_output_seq_id(std::forward <t_output_seq_id>(output_seq_id)),
			m_ref_id_separator(std::forward <t_ref_id_separator>(ref_id_separator)),
			m_gap_opening_cost(gap_opening_cost),
			m_gap_extension_cost(gap_extension_cost),
			m_should_consider_primary_alignments_only(should_consider_primary_alignments_only),
			m_should_use_read_base_qualities(should_use_read_base_qualities)
		{
			for (auto &task : m_task_queue.values())
				task.m_input_processor = this;
		}
		
		void process_input();
		void output_records(project_task_type &task);
		void finish();
		
		panvc3::msa_index &msa_index() { return m_msa_index; }
		panvc3::msa_index const &msa_index() const { return m_msa_index; }
		input_type &alignment_input() { return m_aln_input; }
		input_type const &alignment_input() const { return m_aln_input; }
		dispatch_queue_t output_dispatch_queue() { return *m_output_dispatch_queue; }
		std::string const &reference_id_separator() const { return m_ref_id_separator; }
		std::string const &output_sequence_id() const { return m_output_seq_id; }
		std::int32_t gap_opening_cost() const { return m_gap_opening_cost; }
		std::int32_t gap_extension_cost() const { return m_gap_extension_cost; }
		sequence_vector const &reference_sequence() const { return m_ref_sequence; }
		bool should_use_read_base_qualities() const { return m_should_use_read_base_qualities; }
	};
	
	
	template <typename t_aln_input, typename t_aln_output>
	void input_processor <t_aln_input, t_aln_output>::process_input()
	{
		static_assert(0 < QUEUE_SIZE);
		
		lb::dispatch_ptr <dispatch_group_t> dispatch_group(dispatch_group_create());
		auto parallel_dispatch_queue(dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0));
		
		auto task_idx(m_task_queue.pop_index()); // Reserve one task.
		for (auto &&[rec_idx, aln_rec] : rsv::enumerate(m_aln_input))
		{
			if (0 == (1 + rec_idx) % 10'000'000)
				lb::log_time(std::cerr) << "Processed " << (1 + rec_idx) << " alignments…\n";
			
			auto const flags(aln_rec.flag());
			if (lb::to_underlying(flags & (
				seqan3::sam_flag::unmapped					|
				seqan3::sam_flag::failed_filter				|
				seqan3::sam_flag::duplicate					|
				seqan3::sam_flag::supplementary_alignment
			))) // Ignore unmapped, filtered, duplicate and supplementary.
			{
				++m_statistics.flags_not_matched;
				continue;
			}
			
			// Ignore secondary if requested.
			if (m_should_consider_primary_alignments_only && lb::to_underlying(flags & seqan3::sam_flag::secondary_alignment))
			{
				++m_statistics.flags_not_matched;
				continue;
			}
			
			auto const &ref_id(aln_rec.reference_id());
			if (!ref_id.has_value())
			{
				++m_statistics.ref_id_missing;
				continue;
			}
			
			auto const rec_ref_pos(aln_rec.reference_position());
			if (!rec_ref_pos.has_value())
			{
				++m_statistics.flags_not_matched;
				continue;
			}
			
			if (*rec_ref_pos < 0)
			{
				++m_statistics.flags_not_matched;
				continue;
			}
			
			++m_statistics.matched_reads;
			
			// Check if records can be added to the current task.
			if (m_task_queue[task_idx].is_full())
			{
				// Process the current task.
				auto &current_task(m_task_queue[task_idx]);
				lb::dispatch(current_task).template group_async <&project_task_type::process>(*dispatch_group, parallel_dispatch_queue);
				
				// Get an empty task.
				task_idx = m_task_queue.pop_index();
			}
			
			// Now there is guaranteed to be space in the current task.
			{
				auto &current_task(m_task_queue[task_idx]);
				
				using std::swap;
				swap(aln_rec, current_task.next_record());
			}
		}
		
		// Finish the last task if needed.
		{
			auto &last_task(m_task_queue[task_idx]);
			if (!last_task.empty())
				lb::dispatch(last_task).template group_async <&project_task_type::process>(*dispatch_group, parallel_dispatch_queue);
		}
		
		// When all the work in the group has been completed,
		// the record output blocks have already been inserted to the serial queue.
		lb::dispatch(*this).template group_notify <&input_processor::finish>(*dispatch_group, *m_output_dispatch_queue);
	}
	
	
	template <typename t_input_processor>
	void project_task <t_input_processor>::process()
	{
		typedef typename t_input_processor::input_type			input_type;
		typedef typename input_type::traits_type				input_traits_type;
		typedef typename input_traits_type::sequence_alphabet	sequence_alphabet;
		
		auto const &msa_index(m_input_processor->msa_index());
		auto const &ref_ids(m_input_processor->alignment_input().header().ref_ids()); // ref_ids() not const.
		auto const &ref_id_separator(m_input_processor->reference_id_separator());
		auto const &ref_seq(m_input_processor->reference_sequence());
		auto const gap_opening_cost(m_input_processor->gap_opening_cost());
		auto const gap_extension_cost(m_input_processor->gap_extension_cost());
		auto const should_use_read_base_qualities(m_input_processor->should_use_read_base_qualities());
		panvc3::alignment_projector alignment_projector;
		
		// Process the records.
		// Try to be efficient by caching the previous pointer.
		typedef typename input_type::ref_id_type		ref_id_type;
		ref_id_type prev_ref_id{}; // std::optional.
		panvc3::msa_index::sequence_entry_vector::const_iterator src_seq_entry_it{};
		panvc3::msa_index::sequence_entry_vector::const_iterator dst_seq_entry_it{};
		for (auto &aln_rec : alignment_records())
		{
			auto const &ref_id_(aln_rec.reference_id());
			libbio_assert(ref_id_.has_value());
			
			// Check if we need to find the entries for this sequence.
			if (prev_ref_id != ref_id_)
			{
				prev_ref_id = ref_id_;
				std::string_view const &ref_id(ref_ids[*ref_id_]);
				
				auto const pos(ref_id.find(ref_id_separator));
				if (pos == std::string_view::npos)
				{
					std::cerr << "ERROR: Unable to find the separator “" << ref_id_separator << "” in the RNAME “" << ref_id << "”." << std::endl;
					std::exit(EXIT_FAILURE);
				}
				
				auto const chr_id(ref_id.substr(0, pos));
				auto const src_seq_id(ref_id.substr(1 + pos));
				
				auto const &chr_entry(find_chr_entry(msa_index.chr_entries, chr_id));
				src_seq_entry_it = find_sequence_entry_(chr_entry.sequence_entries, src_seq_id);
				dst_seq_entry_it = find_sequence_entry_(chr_entry.sequence_entries, m_input_processor->output_sequence_id());
			}
			
			auto const &src_seq_entry(*src_seq_entry_it);
			auto const &dst_seq_entry(*dst_seq_entry_it);
			
			// Rewrite the CIGAR.
			auto const src_pos(*aln_rec.reference_position());
			auto const &query_seq(aln_rec.sequence());
			auto const &cigar_seq(aln_rec.cigar_sequence());
			
			auto const dst_pos(alignment_projector.project_alignment(
				src_pos,
				src_seq_entry,
				dst_seq_entry,
				ref_seq,
				query_seq,
				cigar_seq,
				aln_rec.base_qualities(),
				gap_opening_cost,
				gap_extension_cost
			));
		}
		
		// Continue in the output queue.
		lb::dispatch(*this).template async <&project_task::output>(m_input_processor->output_dispatch_queue());
	}
	
	
	template <typename t_input_processor>
	void project_task <t_input_processor>::output()
	{
		// Now we are in the correct thread and also able to pass parameters to functions
		// without calling malloc, since we do not need to call via libdispatch.
		m_input_processor->output_records(*this);
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	void input_processor <t_aln_input, t_aln_output>::output_records(project_task_type &task)
	{
		// Not thread-safe; needs to be executed in a serial queue.
		for (auto &aln_rec : task.alignment_records())
			m_aln_output.push_back(aln_rec); // Needs non-const aln_rec. (Not sure why.)
		
		// Clean up.
		task.reset();
		m_task_queue.push(task);
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	void input_processor <t_aln_input, t_aln_output>::finish()
	{
		std::cout << std::flush;
		lb::log_time(std::cerr) << "Done." << std::endl;
		// FIXME: output the statistics.
		std::exit(EXIT_SUCCESS);
	}
	
	
	auto open_alignment_input_file(fs::path const &path)
	{
		return seqan3::sam_file_input(path);
	}
	
	
	template <typename t_format>
	auto open_alignment_input_stream(std::istream &stream, t_format &&format)
	{
		return seqan3::sam_file_input(stream, std::forward <t_format>(format));
	}	


	template <typename t_aln_input>
	void process_(
		t_aln_input &&aln_input,
		gengetopt_args_info const &args_info
	)
	{
		libbio_assert(args_info.ref_id_separator_arg);
		
		auto &aln_input_header(aln_input.header()); // ref_ids() not const.
		auto const &input_ref_ids(aln_input_header.ref_ids());
		
		// Type aliases for input and output.
		typedef std::remove_cvref_t <decltype(input_ref_ids)>			ref_ids_type;
		typedef std::remove_cvref_t <t_aln_input>						input_type;
		typedef seqan3::sam_file_output <
			typename input_type::selected_field_ids,
			seqan3::type_list <seqan3::format_sam, seqan3::format_bam>,
			ref_ids_type
		>																output_type;
		
		// Set up the reference name.
		ref_ids_type output_ref_ids(input_ref_ids.size());
		auto const *output_seq_id(args_info.output_seq_id_arg ?: args_info.reference_id_arg);
		for (auto &ref_id : output_ref_ids)
			ref_id = output_seq_id;
		
		// Open the alignment output file.
		output_type aln_output(
			std::cout,
			output_ref_ids,
			aln_input_header.ref_id_info | rsv::transform([](auto const &tup){ return std::get <0>(tup); }), // Reference lengths.
			seqan3::format_bam{}
		);
		
		// Load the MSA index.
		panvc3::msa_index msa_index;
		lb::log_time(std::cerr) << "Loading the MSA index…\n";
		{
			lb::file_istream stream;
			lb::open_file_for_reading(args_info.msa_index_arg, stream);
			cereal::PortableBinaryInputArchive archive(stream);
			archive(msa_index);
		}
		
		// Load the reference sequence.
		lb::log_time(std::cerr) << "Loading the reference sequence…\n";
		std::vector <char> ref_sequence;
		lb::read_single_fasta_sequence(args_info.reference_arg, ref_sequence, args_info.reference_id_arg);
		
		// Process the input.
		static input_processor <input_type, output_type> processor(
			std::move(msa_index),
			std::move(aln_input),
			std::move(aln_output),
			std::move(ref_sequence),
			std::move(output_ref_ids),
			args_info.reference_id_arg,
			args_info.reference_msa_id_arg ?: args_info.reference_id_arg,
			output_seq_id,
			args_info.ref_id_separator_arg,
			args_info.gap_opening_cost_arg,
			args_info.gap_extension_cost_arg,
			args_info.primary_only_flag,
			args_info.use_read_base_qualities_flag
		);
		
		processor.process_input();
	}
	
	
	void process(gengetopt_args_info const &args_info)
	{
		// Open the SAM input.
		if (args_info.alignments_arg)
		{
			fs::path const alignments_path(args_info.alignments_arg);
			auto aln_input(open_alignment_input_file(alignments_path));
			process_(aln_input, args_info);
		}
		else
		{
			auto aln_input(open_alignment_input_stream(std::cin, seqan3::format_sam{}));
			process_(aln_input, args_info);
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	process(args_info);
	dispatch_main();
	
	// Not reached.
	return EXIT_SUCCESS;
}
