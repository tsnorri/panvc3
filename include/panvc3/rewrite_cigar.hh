/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_REWRITE_CIGAR_HH
#define PANVC3_REWRITE_CIGAR_HH

#include <libbio/assert.hh>
#include <panvc3/cigar.hh>
#include <panvc3/msa_index.hh>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <vector>


namespace panvc3 {
	
	template <typename t_cigar, typename t_query_seq, typename t_dst_seq>
	std::size_t rewrite_cigar(
		std::size_t src_pos,						// Position in the source reference.
		t_cigar const &cigar_seq,
		msa_index::sequence_entry const &src_entry,
		msa_index::sequence_entry const &dst_entry,
		t_query_seq const &query_seq,				// Typically std::vector <seqan3::dna5>
		t_dst_seq const &dst_seq,
		cigar_buffer &destination
	)
	{
		// Algorithm
		// ---------
		// Iterate over the CIGAR operations as follows.
		// 1. Check for an operation that does nothing, i.e. H or P and keep as-is (whole length).
		// 2. Check for an operation that consumes only query, i.e. I or S and keep as-is (whole length).
		// 3. Handle the remaining operations that consume reference as follows, one character at a time.
		//    1. Calculate the number of nucleotides between the previous and current MSA positions of
		//       the target sequence. In the source, this is always zero (exclusive; since the reference
		//       positions cannot be skipped in the alignments) but additional nucleotides in the target
		//       correspond to deletions in the alignments.
		//    2. Check the character in the target sequence at the current MSA position. This is done by
		//       calculating the MSA position with select on the bit vector built from the source sequence
		//       and reading the value in the bit vector of the target sequence.
		//       a) Sp. the value is 1 (gap). If the alignment operation is D (or N), the operation is
		//          removed. If the operations is M, =, or X, it is changed to I.
		//       b) Sp. the value is 0 (non-gap). If the alignment operation is D (or N), it is kept as-is.
		//          If the operation is M, =, or X, the query character is compared to the target sequence
		//          character and the operation is changed to = or X accordingly.
		
		using seqan3::operator""_cigar_operation;
		
		std::size_t query_pos{};															// Position in the query (read).
		auto const aln_pos(src_entry.gap_positions_select0_support(1 + src_pos));			// Convert to an aligned position.
		auto prev_excess(dst_entry.gap_positions_rank0_support(aln_pos));					// May be zero.
		
		// prev_excess has the number of non-gaps in the target sequence up to and excluding the current position.
		// Hence it is also the zero-based index of the first non-gap character beginning from the current position.
		auto const retval(prev_excess);
		
		destination.clear();
		for (auto const &cigar_item : cigar_seq)
		{
			using seqan3::get;
			
			auto const op_count(get <0>(cigar_item));
			auto const operation(get <1>(cigar_item));
			auto const operation_(operation.to_char());
			
			switch (operation_)
			{
				case 'I':	// Insertion, consumes query.
				case 'S':	// Soft clipping, consumes query.
					destination.push_back(operation, op_count);
					query_pos += op_count;
					break;
					
				case 'H':	// Hard clipping, consumes nothing.
				case 'P':	// Padding (silent deletion from padded reference), consumes nothing.
					destination.push_back(operation, op_count);
					break;
				
				case 'M':	// Match or mismatch, consumes both.
				case '=':	// Match, consumes both.
				case 'X':	// Mismatch, consumes both.
				{
					// Process one character at a time.
					for (std::uint32_t i(0); i < op_count; ++i)
					{
						// Add a deletion if needed.
						auto const aln_pos(src_entry.gap_positions_select0_support(1 + src_pos));	// Convert to an aligned position.
						auto const excess(dst_entry.gap_positions_rank0_support(aln_pos));
						if (prev_excess < excess) // Add a deletion if there was a non-gap character between the current position and the previous one.
							destination.push_back('D'_cigar_operation, excess - prev_excess);
						
						if (1 == dst_entry.gap_positions[aln_pos])
						{
							// Gap; change to I.
							destination.push_back('I'_cigar_operation, 1);
							prev_excess = excess;
						}
						else
						{
							// Non-gap; excess is equal to the unaligned position of the dst. character.
							auto const query_cc(query_seq[query_pos]);
							auto const dst_cc(dst_seq[excess]);
							
							if (query_cc == dst_cc)
								destination.push_back('='_cigar_operation, 1);
							else
								destination.push_back('X'_cigar_operation, 1);
							
							prev_excess = excess + 1;
						}
						
						++query_pos;
						++src_pos;
					}
					break;
				}
				
				case 'D':	// Deletion, consumes reference.
				case 'N':	// Skipped region, consumes reference. (In SAMv1, this is only relevant in mRNA-to-genome alignments.)
				{
					// Process one character at a time.
					for (std::uint32_t i(0); i < op_count; ++i)
					{
						// Add a deletion if needed.
						auto const aln_pos(src_entry.gap_positions_select0_support(1 + src_pos));	// Convert to an aligned position.
						auto const excess(dst_entry.gap_positions_rank0_support(aln_pos));
						if (prev_excess < excess) // Add a deletion if there was a non-gap character between the current position and the previous one.
							destination.push_back('D'_cigar_operation, excess - prev_excess);
						
						if (0 == dst_entry.gap_positions[aln_pos])
						{
							// Destination has non-gap. (Gaps are ignored.)
							destination.push_back('D'_cigar_operation, 1);
							prev_excess = excess + 1;
						}
						else
						{
							prev_excess = excess;
						}
						
						++src_pos;
					}
					break;
				}
				
				default:
					libbio_fail("Unexpected CIGAR operation “", operation_, "”");
					break;
			}
		}
		
		destination.finish();
		return retval;
	}
}

#endif
