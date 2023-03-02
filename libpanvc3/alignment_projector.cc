/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/align.hh>
#include <panvc3/alignment_projector.hh>
#include <panvc3/rewrite_cigar.hh>
#include <panvc3/utility.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/iterator/insert_iterators.hpp>
#include <range/v3/view/transform.hpp>

namespace rsv	= ranges::views;


namespace panvc3 {
	
	void alignment_projector::reset()
	{
		// indel_run_checker::reset() called in project_alignment().
		m_cigar_realigned.clear();
		m_rewrite_buffer.clear();
		m_realign_buffer.clear();
		m_realigned_reference_ranges.clear();
		m_realigned_query_ranges.clear();
	}


	std::size_t alignment_projector::project_alignment(
		std::size_t const src_pos,
		msa_index::sequence_entry const &src_seq_entry,
		msa_index::sequence_entry const &dst_seq_entry,
		std::vector <char> const &ref_seq,
		std::vector <sequence_alphabet> const &query_seq,
		std::vector <seqan3::cigar> const &cigar_seq,
		std::vector <quality_alphabet> const &base_qualities,
		std::int32_t const gap_opening_cost,
		std::int32_t const gap_extension_cost
	)
	{
		auto const dst_pos(rewrite_cigar(
			src_pos,
			cigar_seq,
			src_seq_entry,
			dst_seq_entry,
			query_seq | rsv::transform([](auto const dna_cc){ return dna_cc.to_char(); }),
			ref_seq,
			m_rewrite_buffer
		));
		
		// Realign where needed, i.e. if there are runs of adjacent insertions and deletions.
		auto const &cigar_rewritten(m_rewrite_buffer.operations());
		m_indel_run_checker.reset(cigar_rewritten, dst_pos);
		m_cigar_realigned.clear();
		
		auto cigar_begin(cigar_rewritten.begin());
		while (m_indel_run_checker.find_next_range_for_realigning())
		{
			// Copy the non-realigned part to cigar_realigned.
			auto const realn_range(m_indel_run_checker.cigar_realigned_range());
			ranges::copy(
				ranges::subrange(cigar_begin, realn_range.first),
				ranges::back_inserter(m_cigar_realigned)
			);
			cigar_begin = realn_range.second;
			
			// Store the realigned ranges.
			auto const ref_pos(m_indel_run_checker.reference_position());
			auto const ref_range(m_indel_run_checker.reference_range()); // Has segment-relative position.
			auto const query_range(m_indel_run_checker.query_range());
			m_realigned_reference_ranges.emplace_back(ref_pos, ref_range.length);
			m_realigned_query_ranges.emplace_back(query_range);
			
			// Realign the found range.
			if (base_qualities.empty())
			{
				auto ref_part(ref_seq | ref_range.slice() | rsv::transform([](auto const cc){
					sequence_alphabet scc;
					scc.assign_char(cc);
					return scc;
				}));
				auto query_part(query_seq | query_range.slice());
				align_global <false>(
					ref_part,
					query_part,
					gap_opening_cost,
					gap_extension_cost,
					m_realign_buffer
				);
			}
			else
			{
				typedef seqan3::qualified <sequence_alphabet, quality_alphabet>	qualified_alphabet;
				auto ref_part(
					ref_seq
					| ref_range.slice()
					| rsv::transform([](auto const cc) -> qualified_alphabet {
						sequence_alphabet scc;
						scc.assign_char(cc);
						qualified_alphabet retval;
						retval = scc;
						retval = max_letter <quality_alphabet>();	// Assign maximum qualities to the reference.
						return retval;								// (We could take the variant likelihoods into account, though.)
					})
				);
				auto query_part(
					rsv::zip(
						query_seq,
						base_qualities
					)
					| query_range.slice()
					| rsv::transform([](auto const tup) -> qualified_alphabet {
						qualified_alphabet retval;
						retval = std::get <0>(tup);
						retval = std::get <1>(tup);
						return retval;
					})
				);
			
				align_global <true, sequence_alphabet, quality_alphabet>(
					ref_part,
					query_part,
					gap_opening_cost,
					gap_extension_cost,
					m_realign_buffer
				);
			}
		
			// Copy the new operations.
			auto const &realigned_part(m_realign_buffer.operations());
			m_cigar_realigned.insert(m_cigar_realigned.end(), realigned_part.begin(), realigned_part.end());
		}
	
		// Copy to cigar_realigned.
		ranges::copy(
			ranges::subrange(cigar_begin, cigar_rewritten.end()),
			ranges::back_inserter(m_cigar_realigned)
		);

		// Collapse the operations.
		collapse_cigar_operations(m_cigar_realigned);
		
		return dst_pos;
	}
}
