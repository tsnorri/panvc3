/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <deque>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/bits.hh>
#include <libbio/file_handling.hh>
#include <range/v3/all.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include "cmdline.h"

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
	
	
	template <typename t_aln_record>
	std::size_t calculate_record_length(t_aln_record const &aln_record, bool const should_include_soft_clipping)
	{
		std::size_t reference_length{};
		
		for (auto const &cigar_item : aln_record.cigar_sequence())
		{
			// For some reason the correct seqan3::get does not get used when accessing cigar_item
			// with it. The type coversion operator does work, though.
			// FIXME: check if this is still true.
			std::uint32_t const op_count(cigar_item);
			seqan3::cigar::operation const operation(cigar_item);
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
	
	
	template <typename t_aln_input>
	void calculate_coverage_(
		t_aln_input &aln_input,
		gengetopt_args_info const &args_info
	)
	{
		typedef typename t_aln_input::traits_type				alignment_input_traits_type;
		typedef typename alignment_input_traits_type::ref_ids	reference_ids_type;
		
		auto const contig_prefix(args_info.contig_prefix_arg);
		bool const should_include_clipping(args_info.include_clipping_flag);
		bool const should_consider_primary_alignments_only(args_info.primary_only_flag);
		bool const requires_same_config_prefix_in_next(args_info.same_ref_flag);
		
		reference_ids_type const &reference_ids(aln_input.header().ref_ids());
		
		std::size_t reads_processed{};
		std::size_t flags_not_matched{};
		std::size_t ref_id_mismatches{};
		std::size_t mate_ref_id_mismatches{};
		
		std::multiset <interval, interval_end_pos_cmp> coverage_left;
		std::multiset <interval, interval_end_pos_cmp> coverage_right;
		
		lb::log_time(std::cerr) << "Calculating coverage…\n";
		
		std::size_t prev_record_pos{};
		for (auto const &aln_rec : aln_input)
		{
			++reads_processed;
			if (0 == reads_processed % 10000000)
				lb::log_time(std::cerr) << "Processed " << reads_processed << " reads…\n";
			
			// If POS is zero, skip the record.
			auto const flags(aln_rec.flag());
			
			if (lb::to_underlying(flags & (
				seqan3::sam_flag::unmapped					|
				seqan3::sam_flag::failed_filter				|
				seqan3::sam_flag::duplicate					|
				seqan3::sam_flag::supplementary_alignment
			))) // Ignore unmapped, filtered, duplicate and supplementary.
			{
				++flags_not_matched;
				continue;
			}
			
			// Ignore secondary if requested.
			if (should_consider_primary_alignments_only && lb::to_underlying(flags & seqan3::sam_flag::secondary_alignment))
			{
				++flags_not_matched;
				continue;
			}
			
			// Check the contig name if requested.
			auto const ref_id(aln_rec.reference_id());
			if (contig_prefix)
			{
				// Check the contig prefix for the current alignment.
				if (!ref_id.has_value())
				{
					++ref_id_mismatches;
					continue;
				}
				
				if (!reference_ids[*ref_id].starts_with(contig_prefix))
				{
					++ref_id_mismatches;
					continue;
				}
				
				// Check the contig prefix of the next primary alignment.
				// (Currently we do not try to determine the reference ID of all the possible alignments
				// of the next read.)
				if (requires_same_config_prefix_in_next)
				{
					auto const mate_ref_id(aln_rec.mate_reference_id());
					if (!mate_ref_id.has_value())
					{
						++mate_ref_id_mismatches;
						continue;
					}
					
					if (!reference_ids[*mate_ref_id].starts_with(contig_prefix))
					{
						++mate_ref_id_mismatches;
						continue;
					}
				}
			}
			
			auto const rec_ref_pos_(aln_rec.reference_position());
			if (!rec_ref_pos_.has_value())
			{
				++flags_not_matched;
				continue;
			}
			
			// The following assertion checks that rec_ref_pos is non-negative, too.
			libbio_always_assert_lte(prev_record_pos, *rec_ref_pos_); // FIXME: error message: alignments need to be sorted by position.
			std::size_t const rec_ref_pos(*rec_ref_pos_);
			auto const ref_length(calculate_record_length(aln_rec, should_include_clipping));
			auto const rec_ref_end(rec_ref_pos + ref_length);
			
			// Check if the position was incremented and there are intervals to report.
			if (prev_record_pos < rec_ref_pos)
			{
				// Move the intervals from the right to the left.
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
				
				if (!coverage_left.empty())
				{
					std::cout << prev_record_pos << '\t' << coverage_left.size() << '\n';

					// Remove from the left.
					// Find the intervals where the end position is greater than the current position,
					// i.e. ones that cover the current position.
					auto it(coverage_left.upper_bound(rec_ref_pos));
				
					// Erase the rest.
					coverage_left.erase(coverage_left.begin(), it);
				}
			}
			
			// Add to the right.
			if (ref_length)
				coverage_right.emplace(rec_ref_pos, rec_ref_end);
			
			prev_record_pos = rec_ref_pos; // We could also consider some of the filtered reads’ positions for m_prev_record_pos, but I think this will suffice.
		}
		
		std::cout << prev_record_pos << '\t' << coverage_left.size() << '\n';
		lb::log_time(std::cerr) << "Done.\n";
	}
	
	
	void calculate_coverage(gengetopt_args_info const &args_info)
	{
		// Open the SAM input. We expect the alignements to have been sorted by the leftmost co-ordinate.
		if (args_info.alignments_arg)
		{
			fs::path const alignments_path(args_info.alignments_arg);
			auto aln_input(open_alignment_input_file(alignments_path));
			calculate_coverage_(aln_input, args_info);
		}
		else
		{
			auto aln_input(open_alignment_input_stream(std::cin, seqan3::format_sam{}));
			calculate_coverage_(aln_input, args_info);
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
		if (!args_info.contig_prefix_arg)
		{
			std::cerr << "ERROR: --same-ref requires --contig-prefix." << std::endl;
			return EXIT_FAILURE;
		}
		
		if (!args_info.primary_only_flag)
		{
			std::cerr << "ERROR: --same-ref currently requires --primary-only." << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	if (args_info.coverage_given)
		calculate_coverage(args_info);
	else
	{
		std::cerr << "ERROR: No mode given." << std::endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
