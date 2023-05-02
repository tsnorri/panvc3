/*
 * Copyright (c) 2022-2023 Tuukka Norri
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
#include <range/v3/all.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include "cmdline.h"

namespace accs	= boost::accumulators;
namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace rsv	= ranges::views;


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
		std::size_t flags_not_matched{};
		std::size_t ref_id_mismatches{};
		std::size_t mate_ref_id_mismatches{};
	};
	
	
	template <typename t_aln_record>
	std::size_t calculate_record_length(t_aln_record const &aln_record, bool const should_include_soft_clipping)
	{
		std::size_t reference_length{};
		
		for (auto const &cigar_item : aln_record.cigar_sequence())
		{
			using seqan3::get;
			
			auto const op_count(get <0>(cigar_item));
			auto const operation(get <1>(cigar_item));
			auto const cc(operation.to_char());
			
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
	
	
	constexpr static auto sam_reader_fields()
	{
		return seqan3::fields <
			seqan3::field::ref_id,
			seqan3::field::ref_offset,
			seqan3::field::mapq,
			seqan3::field::mate,
			seqan3::field::flag,
			seqan3::field::cigar
		>{};
	}
	
	
	static auto open_alignment_input_file(fs::path const &path)
	{
		return seqan3::sam_file_input(path, sam_reader_fields());
	}
	
	
	template <typename t_format>
	static auto open_alignment_input_stream(std::istream &stream, t_format &&format)
	{
		return seqan3::sam_file_input(stream, std::forward <t_format>(format), sam_reader_fields());
	}


	template <typename t_aln_input, typename t_cb>
	void process_alignments(t_aln_input &&aln_input, gengetopt_args_info const &args_info, alignment_statistics &stats, t_cb &&cb)
	{
		typedef typename std::remove_cvref_t <t_aln_input>::traits_type	alignment_input_traits_type;
		typedef typename alignment_input_traits_type::ref_ids			reference_ids_type;
		
		bool const should_consider_primary_alignments_only(args_info.primary_only_flag);
		bool const requires_same_contig_in_next(args_info.same_ref_flag);
		
		reference_ids_type const &reference_ids(aln_input.header().ref_ids());
		std::vector <std::size_t> filtered_ref_ids; // Filter == keep.
		std::vector <std::size_t> ref_id_eq_classes(reference_ids.size(), SIZE_MAX);

		// Check if any RNAME (prefixes) were given.
		for (unsigned int i(0); i < args_info.rname_given; ++i)
		{
			auto const *rname(args_info.rname_arg[i]);
			std::string_view const rname_sv(rname);
			
			for (auto const &[ref_id, ref_name] : rsv::enumerate(reference_ids))
			{
				if (ref_name == rname_sv || (args_info.rname_prefix_given && ref_name.starts_with(rname_sv)))
				{
					std::cerr << "Filtering by reference '" << ref_name << "' (" << ref_id << ").\n";
					filtered_ref_ids.push_back(ref_id);
					ref_id_eq_classes[ref_id] = i;
				}
			}
		}
		std::sort(filtered_ref_ids.begin(), filtered_ref_ids.end());

		std::size_t reads_processed{};

		for (auto const &aln_rec : aln_input)
		{
			++reads_processed;
			if (0 == reads_processed % 10000000)
				lb::log_time(std::cerr) << "Processed " << reads_processed << " alignments…\n";
			
			// If POS is zero, skip the record.
			auto const flags(aln_rec.flag());
			
			if (lb::to_underlying(flags & (
				seqan3::sam_flag::unmapped					|
				seqan3::sam_flag::failed_filter				|
				seqan3::sam_flag::duplicate					|
				seqan3::sam_flag::supplementary_alignment
			))) // Ignore unmapped, filtered, duplicate and supplementary.
			{
				++stats.flags_not_matched;
				continue;
			}
			
			// Ignore secondary if requested.
			if (should_consider_primary_alignments_only && lb::to_underlying(flags & seqan3::sam_flag::secondary_alignment))
			{
				++stats.flags_not_matched;
				continue;
			}
			
			// Check the contig name if requested.
			if (!filtered_ref_ids.empty())
			{
				auto const ref_id(aln_rec.reference_id());

				// Check the contig prefix for the current alignment.
				if (!ref_id.has_value())
				{
					++stats.ref_id_mismatches;
					continue;
				}
				
				if (!std::binary_search(filtered_ref_ids.begin(), filtered_ref_ids.end(), *ref_id))
				{
					++stats.ref_id_mismatches;
					continue;
				}
				
				// Check the contig (prefix) of the next primary alignment.
				// (Currently we do not try to determine the reference ID of all the possible alignments
				// of the next read.)
				if (requires_same_contig_in_next)
				{
					auto const mate_ref_id(aln_rec.mate_reference_id());
					if (!mate_ref_id.has_value())
					{
						++stats.mate_ref_id_mismatches;
						continue;
					}
					
					if (ref_id_eq_classes[*ref_id] != ref_id_eq_classes[*mate_ref_id])
					{
						++stats.mate_ref_id_mismatches;
						continue;
					}
				}
			}
			
			auto const rec_ref_pos_(aln_rec.reference_position());
			if (!rec_ref_pos_.has_value())
			{
				++stats.flags_not_matched;
				continue;
			}

			cb(aln_rec);
		}
	}
	
	
	template <typename t_aln_input>
	void calculate_coverage(
		t_aln_input &aln_input,
		alignment_statistics &stats,
		gengetopt_args_info const &args_info
	)
	{
		bool const should_include_clipping(args_info.include_clipping_flag);

		std::size_t prev_record_pos{};
		std::multiset <interval, interval_end_pos_cmp> coverage_left;
		std::multiset <interval, interval_end_pos_cmp> coverage_right;
		
		lb::log_time(std::cerr) << "Calculating coverage…\n";
		std::cout << "POSITION\tCOVERAGE\n";
		
		process_alignments(aln_input, args_info, stats,
			[
				should_include_clipping,
				&prev_record_pos,
				&coverage_left,
				&coverage_right
			](auto const &aln_rec){
				auto const rec_ref_pos_(aln_rec.reference_position());
				
				// The following assertion checks that rec_ref_pos is non-negative, too.
				libbio_always_assert_lte(prev_record_pos, *rec_ref_pos_); // FIXME: error message: alignments need to be sorted by position.
				std::size_t const rec_ref_pos(*rec_ref_pos_);
				auto const ref_length(calculate_record_length(aln_rec, should_include_clipping));
				auto const rec_ref_end(rec_ref_pos + ref_length);
				
				// Check if the position was incremented and there are intervals to report.
				if (prev_record_pos < rec_ref_pos)
				{
					// Move the intervals from the right to the left.
					// This is done first b.c. the left bound of every such record is equal to prev_record_pos and hence
					// they need to be counted when calculating the coverage.
					{
						auto it(coverage_right.begin());
						auto const end(coverage_right.end()); // extract() only invalidates iterators to the extracted element.
						while (it != end)
						{
							libbio_assert_lt(it->pos, rec_ref_pos);
							auto const it_(it);
							++it;
							coverage_left.insert(coverage_left.end(), coverage_right.extract(it_));
						}
					}

					while (!coverage_left.empty() && prev_record_pos < rec_ref_pos)
					{
						// Remove from the left.
						// Find the intervals where the end position is greater than the current position,
						// i.e. ones that cover the current position.
						auto it(coverage_left.upper_bound(prev_record_pos));
					
						// Erase the complement.
						coverage_left.erase(coverage_left.begin(), it);

						std::cout << prev_record_pos << '\t' << coverage_left.size() << '\n';

						++prev_record_pos;
					}
				}
				
				// Add to the right.
				if (ref_length)
					coverage_right.emplace(rec_ref_pos, rec_ref_end);
				
				prev_record_pos = rec_ref_pos; // We could also consider some of the filtered reads’ positions for prev_record_pos, but I think this will suffice.
			}
		);
		
		while (!coverage_left.empty())
		{
			auto it(coverage_left.upper_bound(prev_record_pos));
			coverage_left.erase(coverage_left.begin(), it);
			std::cout << prev_record_pos << '\t' << coverage_left.size() << '\n';
			++prev_record_pos;
		}
	}


	template <typename t_aln_input>
	void count_alignments_by_contig(
		t_aln_input &aln_input,
		alignment_statistics &stats,
		gengetopt_args_info const &args_info
	)
	{
		lb::log_time(std::cerr) << "Counting alignments…\n";
		std::cout << "CONTIG\tCOUNT\n";

		auto const &ref_names_by_id(aln_input.header().ref_ids());
		std::vector <std::size_t> counts_by_ref_id(ref_names_by_id.size(), 0);

		process_alignments(aln_input, args_info, stats,
			[
				&counts_by_ref_id
			](auto const &aln_rec){
				auto const &ref_id(aln_rec.reference_id());
				if (!ref_id.has_value())
					return;

				auto const ref_id_(*ref_id);
				++counts_by_ref_id[ref_id_];
			}
		);

		for (auto const &[idx, name] : rsv::enumerate(ref_names_by_id))
			std::cout << name << '\t' << counts_by_ref_id[idx] << '\n';
	}


	template <typename t_aln_input>
	void mapq_histogram(
		t_aln_input &aln_input,
		alignment_statistics &stats,
		gengetopt_args_info const &args_info
	)
	{
		lb::log_time(std::cerr) << "Calculating the histogram of the MAPQ values…\n";
		std::cout << "VALUE\tCOUNT\n";

		constexpr std::size_t const MAPQ_LIMIT(256);
		std::vector <std::size_t> histogram(256, 0);

		process_alignments(aln_input, args_info, stats,
			[
				&histogram
			](auto const &aln_rec){
				auto const mapq(aln_rec.mapping_quality());
				typedef decltype(mapq) mapq_type;
				static_assert(std::numeric_limits <mapq_type>::max() < MAPQ_LIMIT);
				++histogram[mapq];
			}
		);

		for (auto const &[val, count] : rsv::enumerate(histogram))
		{
			if (0 == count)
				continue;

			std::cout << val << '\t' << count << '\n';
		}
	}


	template <typename t_aln_input>
	void mapq_box_plot(
		t_aln_input &aln_input,
		alignment_statistics &stats,
		gengetopt_args_info const &args_info
	)
	{
		if (args_info.bin_width_arg <= 0)
		{
			std::cerr << "ERROR: Bin width must be positive.\n";
			std::exit(EXIT_FAILURE);
		}
		std::uint64_t const bin_width(args_info.bin_width_arg);
		std::uint64_t current_bin{};

		lb::log_time(std::cerr) << "Calculating the box plot parameters from the MAPQ values…\n";
		std::cerr << "NOTE: Currently the data should be sorted and filtered by chromosome.\n";
		std::cout << "BIN\tMIN\tq10\tq25\tMED\tq75\tq90\tMAX\tCOUNT\n";

		typedef accs::accumulator_set <double, accs::stats <
			accs::tag::count,
			accs::tag::min,
			accs::tag::max,
			accs::tag::extended_p_square_quantile
		>> accumulator_type;
		std::array const probabilities{0.1, 0.25, 0.50, 0.75, 0.9};

		accumulator_type acc(accs::extended_p_square_probabilities = probabilities);

		auto print_cb([&](){
			std::cout
				<< current_bin << '\t'
				<< accs::min(acc) << '\t'
				<< accs::quantile(acc, accs::quantile_probability = 0.1) << '\t'
				<< accs::quantile(acc, accs::quantile_probability = 0.25) << '\t'
				<< accs::quantile(acc, accs::quantile_probability = 0.5) << '\t'
				<< accs::quantile(acc, accs::quantile_probability = 0.75) << '\t'
				<< accs::quantile(acc, accs::quantile_probability = 0.9) << '\t'
				<< accs::max(acc) << '\t'
				<< accs::count(acc) << '\n';
		});

		process_alignments(aln_input, args_info, stats,
			[
				bin_width,
				&current_bin,
				&acc,
				&probabilities,
				&print_cb
			](auto const &aln_rec){
				auto const rec_ref_pos(*aln_rec.reference_position()); // Check done earlier.
				auto const bin(rec_ref_pos / bin_width);
				if (bin != current_bin)
				{
					print_cb();
					current_bin = bin;
					acc = accumulator_type(accs::extended_p_square_probabilities = probabilities);
				}

				auto const mapq(aln_rec.mapping_quality());
				if (mapq < 255) // Ignore invalid values.
					acc(mapq);
			}
		);

		print_cb();
	}
	
	
	template <typename t_cb>
	void read_input(gengetopt_args_info const &args_info, t_cb &&cb)
	{
		// Open the SAM input. We expect the alignements to have been sorted by the leftmost co-ordinate.
		if (args_info.alignments_arg)
		{
			if (args_info.bam_input_flag)
				std::cerr << "WARNING: --bam-input has no effect when reading from file.\n";

			fs::path const alignments_path(args_info.alignments_arg);
			auto aln_input(open_alignment_input_file(alignments_path));
			cb(aln_input);
		}
		else
		{
			if (args_info.bam_input_flag)
			{
				auto aln_input(open_alignment_input_stream(std::cin, seqan3::format_bam{}));
				cb(aln_input);
			}
			else
			{
				auto aln_input(open_alignment_input_stream(std::cin, seqan3::format_sam{}));
				cb(aln_input);
			}
		}
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
			std::cerr << "ERROR: --same-ref requires --rname." << std::endl;
			return EXIT_FAILURE;
		}
		
		if (!args_info.primary_only_flag)
		{
			std::cerr << "ERROR: --same-ref currently requires --primary-only." << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	{
		alignment_statistics stats;
		if (args_info.coverage_given)
			read_input(args_info, [&stats, &args_info](auto &aln_input){ calculate_coverage(aln_input, stats, args_info); });
		else if (args_info.count_alignments_given)
			read_input(args_info, [&stats, &args_info](auto &aln_input){ count_alignments_by_contig(aln_input, stats, args_info); });
		else if (args_info.mapq_histogram_given)
			read_input(args_info, [&stats, &args_info](auto &aln_input){ mapq_histogram(aln_input, stats, args_info); });
		else if (args_info.mapq_box_plot_given)
			read_input(args_info, [&stats, &args_info](auto &aln_input){ mapq_box_plot(aln_input, stats, args_info); });
		else
		{
			std::cerr << "ERROR: No mode given." << std::endl;
			return EXIT_FAILURE;
		}

		std::cerr << "Flags not matched:       " << stats.flags_not_matched << '\n';
		std::cerr << "Ref. ID mismatches:      " << stats.ref_id_mismatches << '\n';
		std::cerr << "Mate ref. ID mismatches: " << stats.mate_ref_id_mismatches << '\n';
	}

	lb::log_time(std::cerr) << "Done.\n";
	
	return EXIT_SUCCESS;
}
