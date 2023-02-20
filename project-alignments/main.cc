/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <iostream>
#include <libbio/algorithm/sorted_set_union.hh>
#include <libbio/dispatch/dispatch_caller.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <panvc3/alignment_projector.hh>
#include <panvc3/spsc_queue.hh>
#include <panvc3/utility.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/iterator/insert_iterators.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/take_exactly.hpp>
#include <range/v3/view/transform.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include "cmdline.h"

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace ios	= boost::iostreams;
namespace rsv	= ranges::views;


namespace {
	
	constexpr inline std::size_t	QUEUE_SIZE{16};
	constexpr inline std::size_t	CHUNK_SIZE{4};

	constexpr static auto const preserved_sam_tags{[]() constexpr {
		using seqan3::operator""_tag;

		// Consider saving UQ.
		std::array retval{
			"AM"_tag,	// The smallest template-independent mapping quality in the template
			"AS"_tag,	// Alignment score generated by the aligner
			"BC"_tag,	// Barcode sequence identifying the sample
			"BZ"_tag,	// Phred quality of the unique molecular barcode bases in the OX tag
			"CB"_tag,	// Cell identifier
			"CO"_tag,	// Free-text comments
			"CR"_tag,	// Cellular barcode sequence bases (uncorrected)
			"CS"_tag,	// Color read sequence
			"CT"_tag,	// Complete read annotation tag, used for consensus annotation dummy features
			"CY"_tag,	// Phred quality of the cellular barcode sequence in the CR tag
			"E2"_tag,	// The 2nd most likely base calls
			"FZ"_tag,	// Flow signal intensities
			"LB"_tag,	// Library
			"MI"_tag,	// Molecular identifier; a string that uniquely identifies the molecule from which the record was derived
			"ML"_tag,	// Base modification problems
			"MM"_tag,	// Base modifications / methylation
			"OQ"_tag,	// Original base quality
			"OX"_tag,	// Original unique molecular barcode bases
			"PG"_tag,	// Program
			"PU"_tag,	// Platform unit
			"QT"_tag,	// Phred quality of the sample barcode sequence in the BC tag
			"QX"_tag,	// Quality score of the unique molecular identifier in the RX tag
			"RX"_tag,	// Sequence bases of the (possibly corrected) unique molecular identifier
			"TS"_tag	// Transcript strand
		};

		// Make sure that the values are sorted.
		std::sort(retval.begin(), retval.end());
		
		return retval;
	}()};


	// Convert a SeqAn 3 SAM tag to a std::array <char, N> where 2 ≤ N.
	// I don't think SeqAn 3 has this utility function.
	// Compare to operator""_tag() in <seqan3/io/sam_file/sam_tag_dictionary.hpp>.
	template <std::size_t t_size>
	constexpr void from_tag(std::uint16_t const val, std::array <char, t_size> &buffer)
	requires (2 <= t_size)
	{
		char const char0(val / 256); // Narrowed automatically when () (not {}) are used.
		char const char1(val % 256);
		std::get <0>(buffer) = char0;
		std::get <1>(buffer) = char1;
	}
	
	
	// Convert a std::array <char, 2> to a SeqAn 3 SAM tag.
	constexpr std::uint16_t to_tag(std::array <char, 2> const &buffer)
	{
		// The tag needs to match /|A-Za-z][A-Za-z0-9]/ (SAMv1, Section 1.5 The alignment section: optional fields),
		// so the values will not be negative.
		libbio_always_assert_lte(std::get <0>(buffer), 127);
		libbio_always_assert_lte(std::get <1>(buffer), 127);
		std::uint16_t retval(std::get <0>(buffer)); // Narrows when () (not {}) are used.
		retval *= 256;
		retval += std::get <1>(buffer);
		return retval;
	}
	
	
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
	
	
	struct realigned_range
	{
		std::string		qname{};
		panvc3::range	range{};

		realigned_range(panvc3::range const &range_):
			range(range_)
		{
		}

		realigned_range(std::string const &qname_, panvc3::range const &range_):
			qname(qname_),
			range(range_)
		{
		}
		
		bool operator<(realigned_range const &other) const { return range.to_tuple() < other.range.to_tuple(); }
		bool operator==(realigned_range const &other) const { return range.to_tuple() == other.range.to_tuple(); }
	};

	typedef std::vector <realigned_range>	realigned_range_vector;
	
	// Precondition: the given range is sorted by range.location.
	realigned_range_vector::const_iterator find_first_overlapping(
		realigned_range_vector::const_iterator it,
		realigned_range_vector::const_iterator const end
	)
	{
		if (it == end)
			return end;
		
		auto it_(it);
		++it;
		for (; it != end; ++it)
		{
			auto const range_end(it_->range.location + it_->range.length);
			if (it->range.location < range_end)
				return it_;
			else if (range_end < it->range.location + it->range.length)
				it_ = it;
		}
		
		return end;
	}
	
	
	enum class project_task_status
	{
		inactive,
		processing,
		finishing
	};
		

	template <typename t_input_processor>
	class project_task
	{
		friend t_input_processor;
		
	public:
		typedef typename t_input_processor::input_record_type	record_type;
		typedef std::array <record_type, CHUNK_SIZE>			record_array;
		typedef typename t_input_processor::tag_count_map		tag_count_map;

	protected:
		t_input_processor					*m_input_processor{};
		record_array						m_records;
		panvc3::alignment_projector			m_alignment_projector;
		tag_count_map						m_removed_tag_counts;
		realigned_range_vector				m_realigned_ranges;
		std::size_t							m_valid_records{};
		std::size_t							m_task_id{};
		std::atomic <project_task_status>	m_status{}; // FIXME: make this conditional.
		bool								m_should_store_realigned_range_qnames{};
		
	public:
		bool is_full() const { return CHUNK_SIZE == m_valid_records; }
		bool empty() const { return 0 == m_valid_records; }
		record_type &next_record() { libbio_assert(project_task_status::inactive == m_status.load()); libbio_assert_lt(m_valid_records, m_records.size()); return m_records[m_valid_records++]; }
		auto alignment_records() { auto const st(m_status.load()); libbio_assert(project_task_status::processing == st || project_task_status::finishing == st); return m_records | rsv::take_exactly(m_valid_records); }
		auto alignment_records() const { auto const st(m_status.load()); libbio_assert(project_task_status::processing == st || project_task_status::finishing == st); return m_records | rsv::take_exactly(m_valid_records); }
		
		void process();
		void output();
		void reset() { m_valid_records = 0; m_removed_tag_counts.clear(); m_realigned_ranges.clear(); }
	};
	
	
	// Declare a virtual destructor in order to make assigning input_processor to a std::unique_ptr easier.
	struct input_processor_base
	{
		virtual ~input_processor_base() {}
		virtual void process_input() = 0;
	};
	
	
	template <typename t_aln_input, typename t_aln_output>
	class input_processor final : public input_processor_base
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
		typedef std::map <std::uint16_t, std::size_t>				tag_count_map;
		typedef panvc3::spsc_queue <project_task_type, QUEUE_SIZE>	queue_type;
		typedef lb::dispatch_ptr <dispatch_queue_t>					dispatch_queue_ptr;
		typedef lb::dispatch_ptr <dispatch_semaphore_t>				semaphore_ptr;
		typedef	void												(*exit_callback_type)(void);
		
	protected:
		panvc3::msa_index			m_msa_index;
		input_type					m_aln_input;
		output_type					m_aln_output;
		lb::file_handle				m_realn_range_handle; // Needs to be before m_realn_range_output due to deallocation order.
		lb::file_ostream			m_realn_range_output;
		sequence_vector				m_ref_sequence;
		
		dispatch_queue_ptr			m_output_dispatch_queue;
		exit_callback_type			m_exit_cb{};
		
		queue_type					m_task_queue{};
		
		alignment_statistics		m_statistics;
		tag_count_map				m_removed_tag_counts;
		realigned_range_vector		m_realigned_ranges;
		realigned_range_vector		m_realigned_range_buffer;
		std::string					m_ref_id;
		std::string					m_msa_ref_id;
		std::string					m_output_ref_id;
		std::string					m_ref_id_separator;
		std::int32_t				m_gap_opening_cost{};
		std::int32_t				m_gap_extension_cost{};
		std::uint16_t				m_realigned_ranges_tag{};
		bool						m_should_consider_primary_alignments_only{};
		bool						m_should_use_read_base_qualities{};
		bool						m_should_keep_duplicate_realigned_ranges{};
		bool						m_should_output_debugging_information{};
		
	public:
		template <
			typename t_ref_id,
			typename t_msa_ref_id,
			typename t_output_ref_id,
			typename t_ref_id_separator
		>
		input_processor(
			panvc3::msa_index			&&msa_index,
			t_aln_input					&&aln_input,
			t_aln_output				&&aln_output,
			sequence_vector				&&ref_sequence,
			lb::file_handle				&&realn_range_handle,
			t_ref_id 					&&ref_id,
			t_msa_ref_id				&&msa_ref_id,
			t_output_ref_id				&&output_ref_id,
			t_ref_id_separator			&&ref_id_separator,
			std::array <char, 2> const	&realigned_ranges_tag,
			std::int32_t				gap_opening_cost,
			std::int32_t				gap_extension_cost,
			bool						should_consider_primary_alignments_only,
			bool						should_use_read_base_qualities,
			bool						should_keep_duplicate_realigned_ranges,
			bool						should_output_debugging_information,
			exit_callback_type			exit_cb
		):
			m_msa_index(std::move(msa_index)),
			m_aln_input(std::move(aln_input)),
			m_aln_output(std::move(aln_output)),
			m_realn_range_handle(std::move(realn_range_handle)),
			m_realn_range_output(m_realn_range_handle.get(), ios::never_close_handle),
			m_ref_sequence(std::move(ref_sequence)),
			m_output_dispatch_queue(dispatch_queue_create("fi.iki.tsnorri.panvc3.project-alignments.output-queue", DISPATCH_QUEUE_SERIAL)),
			m_exit_cb(exit_cb),
			m_ref_id(std::forward <t_ref_id>(ref_id)),
			m_msa_ref_id(std::forward <t_msa_ref_id>(msa_ref_id)),
			m_output_ref_id(std::forward <t_output_ref_id>(output_ref_id)),
			m_ref_id_separator(std::forward <t_ref_id_separator>(ref_id_separator)),
			m_gap_opening_cost(gap_opening_cost),
			m_gap_extension_cost(gap_extension_cost),
			m_realigned_ranges_tag(to_tag(realigned_ranges_tag)),
			m_should_consider_primary_alignments_only(should_consider_primary_alignments_only),
			m_should_use_read_base_qualities(should_use_read_base_qualities),
			m_should_keep_duplicate_realigned_ranges(should_keep_duplicate_realigned_ranges),
			m_should_output_debugging_information(should_output_debugging_information)
		{
			for (auto &task : m_task_queue.values())
			{
				task.m_input_processor = this;
				task.m_should_store_realigned_range_qnames = m_should_output_debugging_information;
			}
		}
		
		void process_input() override;
		void output_records(project_task_type &task);
		void output_realigned_ranges(realigned_range_vector const &ranges, std::size_t const task_id = 0);
		void finish();
		
		panvc3::msa_index &msa_index() { return m_msa_index; }
		panvc3::msa_index const &msa_index() const { return m_msa_index; }
		input_type &alignment_input() { return m_aln_input; }
		input_type const &alignment_input() const { return m_aln_input; }
		output_type &alignment_output() { return m_aln_output; }
		output_type const &alignment_output() const { return m_aln_output; }
		dispatch_queue_t output_dispatch_queue() { return *m_output_dispatch_queue; }
		std::string const &reference_id_separator() const { return m_ref_id_separator; }
		std::string const &output_reference_id() const { return m_output_ref_id; }
		std::int32_t gap_opening_cost() const { return m_gap_opening_cost; }
		std::int32_t gap_extension_cost() const { return m_gap_extension_cost; }
		std::uint16_t realigned_ranges_tag() const { return m_realigned_ranges_tag; }
		sequence_vector const &reference_sequence() const { return m_ref_sequence; }
		bool should_use_read_base_qualities() const { return m_should_use_read_base_qualities; }
		bool should_keep_duplicate_realigned_ranges() const { return m_should_keep_duplicate_realigned_ranges; }
	};


	template <typename t_aln_input, typename t_aln_output>
	void input_processor <t_aln_input, t_aln_output>::process_input()
	{
		if (m_realn_range_output)
		{
			if (m_should_output_debugging_information)
			{
				if (m_should_keep_duplicate_realigned_ranges)
					m_realn_range_output << "Location\tLength\tTask\tQNAME\n";
				else
					m_realn_range_output << "Location\tLength\tQNAME\n";
			}
			else
			{
				m_realn_range_output << "Location\tLength\n";
			}
		}
		
		static_assert(0 < QUEUE_SIZE);
		
		lb::dispatch_ptr <dispatch_group_t> dispatch_group(dispatch_group_create());
		auto parallel_dispatch_queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0));
		
		auto task_idx(m_task_queue.pop_index()); // Reserve one task.
		std::size_t task_id{};
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
				++task_id;
				auto &current_task(m_task_queue[task_idx]);
				current_task.m_task_id = task_id;
				current_task.m_status = project_task_status::processing;
				lb::dispatch(current_task).template group_async <&project_task_type::process>(*dispatch_group, parallel_dispatch_queue);
				
				// Get an empty task.
				task_idx = m_task_queue.pop_index();
			}
			
			// Now there is guaranteed to be space in the current task.
			{
				auto &current_task(m_task_queue[task_idx]);
				libbio_assert(project_task_status::inactive == current_task.m_status.load());
				
				using std::swap;
				swap(aln_rec, current_task.next_record());

				aln_rec.clear();
			}
		}
		
		// Finish the last task if needed.
		{
			auto &last_task(m_task_queue[task_idx]);
			if (!last_task.empty())
			{
				++task_id;
				last_task.m_task_id = task_id;
				lb::dispatch(last_task).template group_async <&project_task_type::process>(*dispatch_group, parallel_dispatch_queue);
			}
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
		
		auto &output_header(m_input_processor->alignment_output().header());
		auto const &msa_index(m_input_processor->msa_index());
		auto const &ref_ids(m_input_processor->alignment_input().header().ref_ids()); // ref_ids() not const.
		auto const &ref_id_separator(m_input_processor->reference_id_separator());
		auto const &ref_seq(m_input_processor->reference_sequence());
		auto const gap_opening_cost(m_input_processor->gap_opening_cost());
		auto const gap_extension_cost(m_input_processor->gap_extension_cost());
		auto const should_use_read_base_qualities(m_input_processor->should_use_read_base_qualities());
		auto const should_keep_duplicate_realigned_ranges(m_input_processor->should_keep_duplicate_realigned_ranges());
		auto const realn_ranges_tag(m_input_processor->realigned_ranges_tag());
		
		// Process the records.
		// Try to be efficient by caching the previous pointer.
		typedef typename input_type::ref_id_type		ref_id_type;
		ref_id_type prev_ref_id{}; // std::optional.
		panvc3::msa_index::sequence_entry_vector::const_iterator src_seq_entry_it{};
		panvc3::msa_index::sequence_entry_vector::const_iterator dst_seq_entry_it{};
		for (auto &aln_rec : alignment_records())
		{
			libbio_assert(project_task_status::processing == m_status.load());

			auto const &ref_id_(aln_rec.reference_id());
			libbio_assert(ref_id_.has_value());
			
			// Check if we need to find the entries for this sequence.
			if (prev_ref_id != ref_id_)
			{
				prev_ref_id = ref_id_;
				std::string_view const ref_id(ref_ids[*ref_id_]);
				
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
				dst_seq_entry_it = find_sequence_entry_(chr_entry.sequence_entries, m_input_processor->output_reference_id());
			}
			
			auto const &src_seq_entry(*src_seq_entry_it);
			auto const &dst_seq_entry(*dst_seq_entry_it);
			
			// Rewrite the CIGAR and the position.
			auto const src_pos(*aln_rec.reference_position());
			auto const &query_seq(aln_rec.sequence());
			auto const &cigar_seq(aln_rec.cigar_sequence());
			
			m_alignment_projector.reset();
			auto const dst_pos(m_alignment_projector.project_alignment(
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

			// Copy the realigned ranges.
			auto const &realn_ranges(m_alignment_projector.realigned_reference_ranges());
			auto const realn_range_count(realn_ranges.size());
			if (realn_range_count)
			{
				m_realigned_ranges.reserve(m_realigned_ranges.size() + realn_range_count);

				if (m_should_store_realigned_range_qnames)
				{
					for (auto const &range : realn_ranges)
						m_realigned_ranges.emplace_back(aln_rec.id(), range);
				}
				else
				{
					for (auto const &range : realn_ranges)
						m_realigned_ranges.emplace_back(range);
				}
				
				if (!should_keep_duplicate_realigned_ranges)
				{
					// Sort by the range and remove duplicates.
					std::sort(m_realigned_ranges.begin(), m_realigned_ranges.end());
					m_realigned_ranges.erase(std::unique(m_realigned_ranges.begin(), m_realigned_ranges.end()), m_realigned_ranges.end());
				}
			}

			// Update the tags.
			auto &tags(aln_rec.tags());
			
			// Store the original alignment.
			// SeqAn 3 does not yet have a type definition for the OA tag which replaces OC,
			// so we use the latter instead.
			{
				using seqan3::get;
				using seqan3::operator""_tag;

				auto &oc_tag(tags.template get <"OC"_tag>());
				oc_tag.clear();
				for (auto const cc : cigar_seq)
				{
					oc_tag.push_back(get <0>(cc));
					oc_tag.push_back(get <1>(cc).to_char());
				}
			}
			
			// Remove the non-preserved tags.
			{
				auto it(tags.begin());
				auto const end(tags.end());
				
				while (it != end)
				{
					// Check if the current tag should be preserved.
					auto const tag(it->first);
					if (std::binary_search(preserved_sam_tags.begin(), preserved_sam_tags.end(), tag))
					{
						++it;
						continue;
					}
					
					// Remove and increment the count.
					auto const it_(it);
					++it; // Not invalidated when std::map::erase() is called.
					tags.erase(it_);
					++m_removed_tag_counts[tag];
				}
			}
			
			// Store the re-aligned ranges in query co-ordinates.
			if (realn_range_count)
			{
				std::vector <std::uint32_t> output_ranges(2 * realn_range_count, 0);
				auto const &realn_query_ranges(m_alignment_projector.realigned_query_ranges());
				for (auto const &[idx, range] : rsv::enumerate(realn_query_ranges))
				{
					auto const idx_(2 * idx);
					output_ranges[idx_] = range.location;
					output_ranges[idx_ + 1] = range.location + range.length;
				}
				
				tags[realn_ranges_tag] = std::move(output_ranges);
			}
			
			// Finally (esp. after setting OA/OC) update the CIGAR, the reference position, and the header pointer.
			aln_rec.reference_position() = dst_pos;
			aln_rec.cigar_sequence() = m_alignment_projector.alignment();
			aln_rec.header_ptr() = &output_header;
		}
		
		// Continue in the output queue.
		libbio_assert(project_task_status::processing == m_status.load());
		m_status = project_task_status::finishing;
		lb::dispatch(*this).template async <&project_task::output>(m_input_processor->output_dispatch_queue());
	}
	
	
	template <typename t_input_processor>
	void project_task <t_input_processor>::output()
	{
		// Now we are in the correct thread and also able to pass parameters to functions
		// without calling malloc, since we do not need to call via libdispatch.
		libbio_assert(project_task_status::finishing == m_status.load());
		m_input_processor->output_records(*this);
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	void input_processor <t_aln_input, t_aln_output>::output_records(project_task_type &task)
	{
		// Not thread-safe; needs to be executed in a serial queue.
		
		libbio_assert(project_task_status::finishing == task.m_status.load());
		
		for (auto &aln_rec : task.alignment_records())
		{
			libbio_assert_eq(aln_rec.sequence().size(), aln_rec.base_qualities().size());
			m_aln_output.push_back(aln_rec); // Needs non-const aln_rec. (Not sure why.)
		}
		
		// Update the removed tag counts.
		// Could be done in O(m + n) time with a specialised merge instead of O(m log n) that we currently have.
		// (I don't think this is significant in any way.)
		for (auto const &kv : task.m_removed_tag_counts)
			m_removed_tag_counts[kv.first] += kv.second;
		
		// Handle the realigned ranges.
		auto const &task_realigned_ranges(task.m_realigned_ranges);
		if (m_realn_range_output)
		{
			if (m_should_keep_duplicate_realigned_ranges)
				output_realigned_ranges(task_realigned_ranges, task.m_task_id);
			else
			{
				// Merge in O(n + m) time since both vectors are already sorted.
				m_realigned_range_buffer.clear();
				m_realigned_range_buffer.reserve(m_realigned_ranges.size() + task_realigned_ranges.size());
				lb::sorted_set_union(m_realigned_ranges, task_realigned_ranges, std::back_inserter(m_realigned_range_buffer));
				
				using std::swap;
				swap(m_realigned_ranges, m_realigned_range_buffer);
			}
		}
		
		// Clean up.
		task.reset();
		task.m_status = project_task_status::inactive;
		m_task_queue.push(task);
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	void input_processor <t_aln_input, t_aln_output>::output_realigned_ranges(realigned_range_vector const &ranges, std::size_t const task_id)
	{
		if (m_should_output_debugging_information)
		{
			if (m_should_keep_duplicate_realigned_ranges)
			{
				for (auto const &rr : ranges)
					m_realn_range_output << rr.range.location << '\t' << rr.range.length << '\t' << task_id << '\t' << rr.qname << '\n';
			}
			else
			{
				for (auto const &rr : ranges)
					m_realn_range_output << rr.range.location << '\t' << rr.range.length << '\t' << rr.qname << '\n';
			}
		}
		else
		{
			for (auto const &rr : ranges)
				m_realn_range_output << rr.range.location << '\t' << rr.range.length << '\n';
		}
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	void input_processor <t_aln_input, t_aln_output>::finish()
	{
		m_aln_output.get_stream() << std::flush;
		std::cout << std::flush;
		
		// Output the sorted realigned ranges if needed.
		if (m_realn_range_output)
		{
			if (!m_should_keep_duplicate_realigned_ranges)
				output_realigned_ranges(m_realigned_ranges);
			
			m_realn_range_output << std::flush;
		}
		
		lb::log_time(std::cerr) << "Done." << std::endl;
		
		// Output the statistics.
		std::cerr << "Matched reads:     " << m_statistics.matched_reads << '\n';
		std::cerr << "Ref. ID missing:   " << m_statistics.ref_id_missing << '\n';
		std::cerr << "Flags not matched: " << m_statistics.flags_not_matched << '\n';
		
		if (m_realn_range_output && !m_should_keep_duplicate_realigned_ranges)
			std::cerr << "Re-aligned ranges: " << m_realigned_ranges.size() << '\n';
		
		if (m_removed_tag_counts.empty())
			std::cerr << "No tags removed.\n";
		else
		{
			std::cerr << "Removed tags:\n";
			std::array <char, 3> buffer{'\0', '\0', '\0'};
			for (auto const &kv : m_removed_tag_counts)
			{
				from_tag(kv.first, buffer);
				std::cerr << '\t' << buffer.data() << ": " << kv.second << '\n';
			}
		}
		
		// Clean up.
		(*m_exit_cb)();
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


	std::unique_ptr <input_processor_base> s_input_processor;


	void do_exit(void *)
	{
		s_input_processor.reset();
		std::exit(EXIT_SUCCESS);
	}


	void do_exit_async()
	{
		dispatch_async_f(dispatch_get_main_queue(), nullptr, &do_exit);
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
		
		// Re-aligned range tag.
		auto const realn_ranges_tag{[&args_info]() -> std::array <char, 2> {
			auto const *tag{args_info.realigned_ranges_tag_arg};
			if (!tag)
			{
				std::cerr << "ERROR: Re-aligned ranges tag not given.\n";
				std::exit(EXIT_FAILURE);
			}
			
			std::regex const tag_regex{"^[XYZ][A-Za-z0-9]$"};
			if (!std::regex_match(tag, tag_regex))
			{
				std::cerr << "ERROR: The given tag for re-aligned ranges does not match the expected format.\n";
				std::exit(EXIT_FAILURE);
			}
			
			return {tag[0], tag[1]};
		}()};
		
		// Find the reference ID.
		auto const input_ref_seq_idx{[&args_info, &aln_input_header]() -> std::size_t {
			auto const &ref_dict(aln_input_header.ref_dict);
			auto const *ref_id(args_info.reference_msa_id_arg ?: args_info.reference_id_arg);
			
			// At least now ref_dict’s equality comparison operator and hasher are not transparent.
			std::span const key(ref_id, std::strlen(ref_id));
			auto const it(ref_dict.find(key));
			
			if (ref_dict.cend() == it)
			{
				std::cerr << "ERROR: Reference with ID ‘" << ref_id << "’ not found in input alignment header.\n";
				std::exit(EXIT_FAILURE);
			}
			
			auto const retval(it->second);
			libbio_always_assert_lte(0, retval);
			return retval;
		}()};
		auto const expected_ref_length([&]() -> std::size_t {
			auto const retval(std::get <0>(aln_input_header.ref_id_info[input_ref_seq_idx]));
			libbio_always_assert_lte(0, retval);
			return retval;
		}());
		
		// Open the alignment output file.
		auto aln_output{[&](){
			// Make sure that aln_output has some header information.
			auto const make_output_type([&]<typename t_fmt>(t_fmt &&fmt){
				ref_ids_type empty_ref_ids;
				return output_type(
					std::cout,
					std::move(empty_ref_ids),
					rsv::empty <std::size_t>(), // Reference lengths; the constructor expects a forward range.
					std::forward <t_fmt>(fmt)
				);
			});
			
			if (args_info.output_bam_flag)
				return make_output_type(seqan3::format_bam{});
			else
				return make_output_type(seqan3::format_sam{});
		}()};
		
		// Copy the relevant headers to the output.
		auto const *output_ref_id(args_info.output_ref_id_arg ?: args_info.reference_id_arg);
		{
			auto &aln_output_header(aln_output.header());
			aln_output_header.comments = aln_input_header.comments;
			aln_output_header.read_groups = aln_input_header.read_groups;
			aln_output_header.ref_id_info = aln_input_header.ref_id_info;
			aln_output_header.ref_dict = aln_input_header.ref_dict;
			
			aln_output_header.ref_ids() = aln_input_header.ref_ids(); // Copy.
			aln_output_header.ref_ids()[input_ref_seq_idx] = output_ref_id;
		}
		
		// Open the realigned range output file if needed.
		lb::file_handle realigned_range_handle;
		if (args_info.output_realigned_ranges_arg)
			realigned_range_handle = lb::file_handle(lb::open_file_for_writing(args_info.output_realigned_ranges_arg, lb::writing_open_mode::CREATE));
		
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
		
		if (expected_ref_length != ref_sequence.size())
		{
			std::cerr << "ERROR: Reference sequence length does not match that in the input alignment header (" << ref_sequence.size() << " vs. " << expected_ref_length << ").\n";
			std::exit(EXIT_FAILURE);
		}
		
		// Process the input.
		s_input_processor = std::make_unique <input_processor <input_type, output_type>>(
			std::move(msa_index),
			std::move(aln_input),
			std::move(aln_output),
			std::move(ref_sequence),
			std::move(realigned_range_handle),
			args_info.reference_id_arg,
			args_info.reference_msa_id_arg ?: args_info.reference_id_arg,
			output_ref_id,
			args_info.ref_id_separator_arg,
			realn_ranges_tag,
			args_info.gap_opening_cost_arg,
			args_info.gap_extension_cost_arg,
			args_info.primary_only_flag,
			args_info.use_read_base_qualities_flag,
			args_info.keep_duplicate_ranges_flag,
			args_info.debugging_output_flag,
			&do_exit_async
		);
		
		lb::log_time(std::cerr) << "Processing the alignments…\n";
		s_input_processor->process_input();
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

	if (args_info.print_invocation_given)
	{
		std::cerr << "Invocation:";
		for (int i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << '\n';
	}

	if (args_info.print_pid_given)
		std::cerr << "PID: " << getpid() << '\n';
	
	process(args_info);
	dispatch_main();
	
	// Not reached.
	return EXIT_SUCCESS;
}
