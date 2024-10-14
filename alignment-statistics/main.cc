/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <sstream> // Needed by boost::accumulators.

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#pragma GCC diagnostic pop

#include <deque>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/bits.hh>
#include <libbio/file_handling.hh>
#include <libbio/utility.hh>			// libbio::is_lt
#include <panvc3/alignment_input.hh>
#include <range/v3/all.hpp>
#include <set>
#include "cmdline.h"

namespace accs		= boost::accumulators;
namespace dispatch	= libbio::dispatch;
namespace lb		= libbio;
namespace rsv		= ranges::views;
namespace sam		= libbio::sam;


namespace {
	
	struct interval
	{
		std::size_t	pos{};
		std::size_t	end_pos{};
		
		interval() = default;
		
		interval(std::size_t const pos_, std::size_t const end_pos_):
			pos(pos_),
			end_pos(end_pos_)
		{
		}
	};
	
	
	struct interval_end_pos_cmp
	{
		using is_transparent = std::true_type;
		
		bool operator()(interval const &lhs, interval const &rhs) const { return lhs.end_pos < rhs.end_pos; }
		bool operator()(std::size_t const lhs, interval const &rhs) const { return lhs < rhs.end_pos; }
		bool operator()(interval const &lhs, std::size_t const rhs) const { return lhs.end_pos < rhs; }
	};


	struct alignment_statistics
	{
		std::uint64_t flags_not_matched{};
		std::uint64_t seq_missing{};
		std::uint64_t ref_id_mismatches{};
		std::uint64_t mate_ref_id_mismatches{};
	};
	
	
	std::size_t calculate_record_length(sam::record const &aln_record, bool const should_include_soft_clipping)
	{
		std::size_t reference_length{};
		
		for (auto const &cigar_item : aln_record.cigar)
		{
			auto const op_count(cigar_item.count());
			auto const operation(cigar_item.operation());
			auto const cc(sam::to_char(operation));
			
			switch (cc)
			{
				case 'M': // Consume both query and reference.
				case '=':
				case 'X':
					reference_length += op_count;
					break;
				
				case 'I': // Consume query.
				case 'S':
					if (should_include_soft_clipping)
						reference_length += op_count;
					break;
					
				case 'D': // Consume reference.
				case 'N':
					reference_length += op_count;
					break;
				
				case 'H':
				case 'P':
					break;
				
				default:
					libbio_fail("Unexpected CIGAR operation “", cc, "”");
					break;
			}
		}
		
		return reference_length;
	}
	
	
	class task_base : public panvc3::alignment_input_delegate
	{
	protected:
		alignment_statistics		m_statistics;
		std::vector <std::size_t>	m_filtered_ref_ids;								// Filter == keep.
		std::vector <std::size_t>	m_ref_id_eq_classes;
		char						**m_rname_arg{};								// Type from gengetopt
		std::size_t					m_reads_processed{};
		unsigned int				m_rname_given{};								// Type from gengetopt
		bool						m_rname_prefix_given{};
		bool						m_should_consider_primary_alignments_only{};
		bool						m_requires_same_contig_in_next{};
		
	public:
		task_base(gengetopt_args_info const &args_info):
			m_rname_arg(args_info.rname_arg),
			m_rname_given(args_info.rname_given),
			m_rname_prefix_given(args_info.rname_prefix_given),
			m_should_consider_primary_alignments_only(args_info.primary_only_flag),
			m_requires_same_contig_in_next(args_info.same_ref_flag)
		{
		}
		
		alignment_statistics &statistics() { return m_statistics; }
		alignment_statistics const &statistics() const { return m_statistics; }
		
		virtual ~task_base() {};
		virtual void prepare() = 0;
		virtual void finish() = 0;
		
		void handle_header(sam::header &header) override; // Needs to be called by subclasses.
		
		bool log_alignment_count_and_check(sam::record const &aln_rec);
		void run(panvc3::alignment_input &input) { prepare(); input.run(); }
	};
	
	
	void task_base::handle_header(sam::header &header)
	{
		// Header owned by panvc3::alignment_input.
		
		auto const &reference_ids(header.reference_sequences);
		m_ref_id_eq_classes.clear();
		m_ref_id_eq_classes.resize(reference_ids.size(), SIZE_MAX);
		
		// Check if any RNAME (prefixes) were given.
		for (unsigned int i(0); i < m_rname_given; ++i)
		{
			auto const *rname(m_rname_arg[i]);
			std::string_view const rname_sv(rname);
			
			for (auto const &[ref_id, ref_entry] : rsv::enumerate(reference_ids))
			{
				auto const &ref_name(ref_entry.name);
				if (ref_name == rname_sv || (m_rname_prefix_given && ref_name.starts_with(rname_sv)))
				{
					std::cerr << "Filtering by reference ‘" << ref_name << "’ (" << ref_id << ").\n";
					m_filtered_ref_ids.push_back(ref_id);
					m_ref_id_eq_classes[ref_id] = i;
				}
			}
		}
		
		std::sort(m_filtered_ref_ids.begin(), m_filtered_ref_ids.end());
	}
	
	
	bool task_base::log_alignment_count_and_check(sam::record const &aln_rec)
	{
		++m_reads_processed;
		if (0 == m_reads_processed % 10000000)
			lb::log_time(std::cerr) << "Processed " << m_reads_processed << " alignments…\n";
	
		// If POS is zero, skip the record.
		auto const flags(aln_rec.flag);
		
		if (std::to_underlying(flags & (
			sam::flag::unmapped					|
			sam::flag::failed_filter			|
			sam::flag::duplicate				|
			sam::flag::supplementary_alignment
		))) // Ignore unmapped, filtered, duplicate and supplementary.
		{
			++m_statistics.flags_not_matched;
			return false;
		}
		
		if (aln_rec.seq.empty())
		{
			++m_statistics.seq_missing;
			return false;
		}
		
		// Ignore secondary if requested.
		if (m_should_consider_primary_alignments_only && std::to_underlying(flags & sam::flag::secondary_alignment))
		{
			++m_statistics.flags_not_matched;
			return false;
		}
		
		// Check the contig name if requested.
		if (!m_filtered_ref_ids.empty())
		{
			auto const ref_id(aln_rec.rname_id);

			// Check the contig prefix for the current alignment.
			if (sam::INVALID_REFERENCE_ID == ref_id)
			{
				++m_statistics.ref_id_mismatches;
				return false;
			}
		
			if (!std::binary_search(m_filtered_ref_ids.begin(), m_filtered_ref_ids.end(), ref_id))
			{
				++m_statistics.ref_id_mismatches;
				return false;
			}
		
			// Check the contig (prefix) of the next primary alignment.
			// (Currently we do not try to determine the reference ID of all the possible alignments
			// of the next read.)
			if (m_requires_same_contig_in_next)
			{
				auto const mate_ref_id(aln_rec.rnext_id);
				if (sam::INVALID_REFERENCE_ID == mate_ref_id)
				{
					++m_statistics.mate_ref_id_mismatches;
					return false;
				}
			
				if (m_ref_id_eq_classes[ref_id] != m_ref_id_eq_classes[mate_ref_id])
				{
					++m_statistics.mate_ref_id_mismatches;
					return false;
				}
			}
		}
		
		auto const rec_ref_pos(aln_rec.pos);
		if (sam::INVALID_POSITION == rec_ref_pos)
		{
			++m_statistics.flags_not_matched;
			return false;
		}
		
		return true;
	}
	
	
	class calculate_coverage_task final : public task_base
	{
	private:
		std::size_t										m_prev_record_pos{};
		std::multiset <interval, interval_end_pos_cmp>	m_coverage_left;
		std::multiset <interval, interval_end_pos_cmp>	m_coverage_right;
		bool											m_should_include_clipping{};
		
	public:
		calculate_coverage_task(gengetopt_args_info const &args_info):
			task_base(args_info),
			m_should_include_clipping(args_info.include_clipping_flag)
		{
		}
		
		void prepare() override;
		void handle_alignment(sam::record &aln_rec) override;
		void finish() override;
	};
	
	
	void calculate_coverage_task::prepare()
	{
		lb::log_time(std::cerr) << "Calculating coverage…\n";
		std::cout << "POSITION\tCOVERAGE\n";
	}
	
	
	void calculate_coverage_task::handle_alignment(sam::record &aln_rec)
	{
		if (!log_alignment_count_and_check(aln_rec))
			return;
		
		auto const rec_ref_pos(aln_rec.pos); // Check done earlier.
		
		// The following assertion checks that rec_ref_pos is non-negative, too.
		libbio_always_assert_lte(m_prev_record_pos, rec_ref_pos); // FIXME: error message: alignments need to be sorted by position.
		auto const ref_length(calculate_record_length(aln_rec, m_should_include_clipping));
		auto const rec_ref_end(rec_ref_pos + ref_length);
		
		// Check if the position was incremented and there are intervals to report.
		if (lb::is_lt(m_prev_record_pos, rec_ref_pos))
		{
			// Move the intervals from the right to the left.
			// This is done first b.c. the left bound of every such record is equal to prev_record_pos and hence
			// they need to be counted when calculating the coverage.
			{
				auto it(m_coverage_right.begin());
				auto const end(m_coverage_right.end()); // extract() only invalidates iterators to the extracted element.
				while (it != end)
				{
					libbio_assert_lt(it->pos, rec_ref_pos);
					auto const it_(it);
					++it;
					m_coverage_left.insert(m_coverage_left.end(), m_coverage_right.extract(it_));
				}
			}

			while (!m_coverage_left.empty() && lb::is_lt(m_prev_record_pos, rec_ref_pos))
			{
				// Remove from the left.
				// Find the intervals where the end position is greater than the current position,
				// i.e. ones that cover the current position.
				auto it(m_coverage_left.upper_bound(m_prev_record_pos));
		
				// Erase the complement.
				m_coverage_left.erase(m_coverage_left.begin(), it);

				std::cout << m_prev_record_pos << '\t' << m_coverage_left.size() << '\n';

				++m_prev_record_pos;
			}
		}
		
		// Add to the right.
		if (ref_length)
			m_coverage_right.emplace(rec_ref_pos, rec_ref_end);
		
		m_prev_record_pos = rec_ref_pos; // We could also consider some of the filtered reads’ positions for prev_record_pos, but I think this will suffice.
	}
	
	
	void calculate_coverage_task::finish()
	{
		while (!m_coverage_left.empty())
		{
			auto it(m_coverage_left.upper_bound(m_prev_record_pos));
			m_coverage_left.erase(m_coverage_left.begin(), it);
			std::cout << m_prev_record_pos << '\t' << m_coverage_left.size() << '\n';
			++m_prev_record_pos;
		}
	}
	
	
	class count_alignments_by_contig_task final : public task_base
	{
	private:
		sam::reference_sequence_entry_vector const	*m_ref_names_by_id{};
		std::vector <std::size_t>					m_counts_by_ref_id;
		bool										m_should_include_clipping{};
		
	public:
		count_alignments_by_contig_task(gengetopt_args_info const &args_info):
			task_base(args_info),
			m_should_include_clipping(args_info.include_clipping_flag)
		{
		}
		
		void prepare() override {}
		void handle_header(sam::header &header) override;
		void handle_alignment(sam::record &aln_rec) override;
		void finish() override;
	};
	
	
	void count_alignments_by_contig_task::handle_header(sam::header &header)
	{
		// Header owned by panvc3::alignment_input.
		task_base::handle_header(header);
		m_ref_names_by_id = &header.reference_sequences;
		m_counts_by_ref_id.resize(m_ref_names_by_id->size(), 0);
	}
	
	
	void count_alignments_by_contig_task::handle_alignment(sam::record &aln_rec)
	{
		if (!log_alignment_count_and_check(aln_rec))
			return;
		
		auto const ref_id(aln_rec.rname_id);
		if (sam::INVALID_REFERENCE_ID == ref_id)
			return;

		++m_counts_by_ref_id[ref_id];
	}
	
	
	void count_alignments_by_contig_task::finish()
	{
		for (auto const &[idx, name] : rsv::enumerate(*m_ref_names_by_id))
			std::cout << name << '\t' << m_counts_by_ref_id[idx] << '\n';
	}
	
	
	struct mapq_histogram_task final : public task_base
	{
	private:
		constexpr static std::size_t const	MAPQ_LIMIT{256};
		std::vector <std::size_t>			m_histogram;
		
	public:
		mapq_histogram_task(gengetopt_args_info const &args_info):
			task_base(args_info),
			m_histogram(256, 0)
		{
		}
		
		void prepare() override;
		void handle_alignment(sam::record &aln_rec) override;
		void finish() override;
	};
	
	
	void mapq_histogram_task::prepare()
	{
		lb::log_time(std::cerr) << "Calculating the histogram of the MAPQ values…\n";
		std::cout << "VALUE\tCOUNT\n";
	}
	
	
	void mapq_histogram_task::handle_alignment(sam::record &aln_rec)
	{
		if (!log_alignment_count_and_check(aln_rec))
			return;
		
		auto const mapq(aln_rec.mapq);
		typedef decltype(mapq) mapq_type;
		static_assert(std::numeric_limits <mapq_type>::max() < MAPQ_LIMIT);
		++m_histogram[mapq];
	}
	
	
	void mapq_histogram_task::finish()
	{
		for (auto const &[val, count] : rsv::enumerate(m_histogram))
		{
			if (0 == count)
				continue;

			std::cout << val << '\t' << count << '\n';
		}
	}
	
	
	class mapq_box_plot_task final : public task_base
	{
	public:
		typedef accs::accumulator_set <double, accs::stats <
			accs::tag::count,
			accs::tag::min,
			accs::tag::max,
			accs::tag::extended_p_square_quantile
		>> accumulator_type;
		
	private:
		constexpr static std::array const		probabilities{0.1, 0.25, 0.50, 0.75, 0.9};
		accumulator_type						m_acc;
		std::uint64_t							m_bin_width{};
		std::uint64_t							m_current_bin{};
		
	public:
		mapq_box_plot_task(gengetopt_args_info const &args_info):
			task_base(args_info),
			m_acc(accs::extended_p_square_probabilities = probabilities),
			m_bin_width(args_info.bin_width_arg)
		{
			if (args_info.bin_width_arg <= 0)
			{
				std::cerr << "ERROR: Bin width must be positive.\n";
				std::exit(EXIT_FAILURE);
			}
		}
		
		void prepare() override;
		void handle_alignment(sam::record &aln_rec) override;
		void finish() override { output_accumulator(); }
		
		void output_accumulator();
	};
	
	
	void mapq_box_plot_task::prepare()
	{
		lb::log_time(std::cerr) << "Calculating the box plot parameters from the MAPQ values…\n";
		std::cerr << "NOTE: Currently the data should be sorted and filtered by chromosome.\n";
		std::cout << "BIN\tMIN\tq10\tq25\tMED\tq75\tq90\tMAX\tCOUNT\n";
	}
	
	
	void mapq_box_plot_task::output_accumulator()
	{
		std::cout
			<< m_current_bin << '\t'
			<< accs::min(m_acc) << '\t'
			<< accs::quantile(m_acc, accs::quantile_probability = 0.1) << '\t'
			<< accs::quantile(m_acc, accs::quantile_probability = 0.25) << '\t'
			<< accs::quantile(m_acc, accs::quantile_probability = 0.5) << '\t'
			<< accs::quantile(m_acc, accs::quantile_probability = 0.75) << '\t'
			<< accs::quantile(m_acc, accs::quantile_probability = 0.9) << '\t'
			<< accs::max(m_acc) << '\t'
			<< accs::count(m_acc) << '\n';
	}
	
	
	void mapq_box_plot_task::handle_alignment(sam::record &aln_rec)
	{
		if (!log_alignment_count_and_check(aln_rec))
			return;
		
		auto const rec_ref_pos(aln_rec.pos); // Check done earlier.
		auto const bin(rec_ref_pos / m_bin_width);
		if (bin != m_current_bin)
		{
			output_accumulator();
			m_current_bin = bin;
			m_acc = accumulator_type(accs::extended_p_square_probabilities = probabilities);
		}

		auto const mapq(aln_rec.mapq);
		if (mapq < 255) // Ignore invalid values.
			m_acc(mapq);
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	if (args_info.same_ref_flag)
	{
		if (!args_info.rname_given)
		{
			std::cerr << "ERROR: --same-ref requires --rname.\n";
			return EXIT_FAILURE;
		}
		
		if (!args_info.primary_only_flag)
		{
			std::cerr << "ERROR: --same-ref currently requires --primary-only.\n";
			return EXIT_FAILURE;
		}
	}
	
	dispatch::thread_pool thread_pool;
	thread_pool.set_min_workers(3); // Needed for the BAM reader.
	
	if (args_info.threads_arg < 0)
	{
		std::cerr << "ERROR: Number of threads must be non-negative.\n";
		return EXIT_FAILURE;
	}
	else if (args_info.threads_arg < 4)
	{
		std::cerr << "INFO: Using four threads.\n";
	}
	else
	{
		thread_pool.set_max_workers(args_info.threads_arg - 1);
	}
	
	dispatch::parallel_queue parallel_queue(thread_pool);
	
	dispatch::group group;
	auto &main_queue(dispatch::main_queue());
	
	std::unique_ptr <task_base> task;
	if (args_info.coverage_given)
		task = std::make_unique <calculate_coverage_task>(args_info);
	else if (args_info.count_alignments_given)
		task = std::make_unique <count_alignments_by_contig_task>(args_info);
	else if (args_info.mapq_histogram_given)
		task = std::make_unique <mapq_histogram_task>(args_info);
	else if (args_info.mapq_box_plot_given)
		task = std::make_unique <mapq_box_plot_task>(args_info);
	else
	{
		std::cerr << "ERROR: No mode given." << std::endl;
		return EXIT_FAILURE;
	}
	
	auto aln_input(panvc3::alignment_input::open_path_or_stdin(
		args_info.alignments_arg,
		parallel_queue,
		main_queue,
		group,
		*task
	));
	
	main_queue.group_async(group, [&task, &aln_input]{
		task->run(aln_input); // Does not block.
	});
	
	group.notify(main_queue, [&task, &main_queue]{
		task->finish();
		lb::log_time(std::cerr) << "Done.\n";
		
		auto const &statistics(task->statistics());
		std::cerr << "Flags not matched:       " << statistics.flags_not_matched << '\n';
		std::cerr << "Sequence missing:        " << statistics.seq_missing << '\n';
		std::cerr << "Ref. ID mismatches:      " << statistics.ref_id_mismatches << '\n';
		std::cerr << "Mate ref. ID mismatches: " << statistics.mate_ref_id_mismatches << '\n';
		
		main_queue.stop();
	});
	
	main_queue.run();
	return EXIT_SUCCESS;
}
