/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_ALIGNMENT_PROJECTOR_HH
#define PANVC3_ALIGNMENT_PROJECTOR_HH

#include <ostream>
#include <panvc3/cigar.hh>
#include <panvc3/indel_run_checker.hh>
#include <panvc3/msa_index.hh>
#include <panvc3/range.hh>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <vector>


namespace panvc3 {
	
	class alignment_projector
	{
	public:
		typedef seqan3::dna5		sequence_alphabet;
		typedef seqan3::phred42		quality_alphabet;
		
	protected:
		indel_run_checker			m_indel_run_checker;
		cigar_vector				m_cigar_realigned;
		cigar_buffer				m_rewrite_buffer;
		cigar_buffer				m_realign_buffer;
		range_vector				m_realigned_reference_ranges;
		range_vector				m_realigned_query_ranges;		// OK b.c. we project one alignment at a time, i.e. there will not be ranges from multiple alignments.
	
	public:
		void reset();
		
		std::size_t project_alignment(
			std::size_t const src_pos,
			msa_index::sequence_entry const &src_seq_entry,
			msa_index::sequence_entry const &dst_seq_entry,
			std::vector <char> const &ref_seq,
			std::vector <sequence_alphabet> const &query_seq,
			std::vector <seqan3::cigar> const &cigar_seq,
			std::vector <quality_alphabet> const &base_qualities,	// Pass empty vector to not use base qualities.
			std::int32_t const gap_opening_cost,
			std::int32_t const gap_extension_cost
		);
		
		cigar_vector const &alignment() const { return m_cigar_realigned; }
		range_vector const &realigned_reference_ranges() const { return m_realigned_reference_ranges; }
		range_vector const &realigned_query_ranges() const { return m_realigned_query_ranges; }
	};
}

#endif
