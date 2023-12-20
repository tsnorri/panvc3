/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_ALIGNMENT_PROJECTOR_HH
#define PANVC3_ALIGNMENT_PROJECTOR_HH

#include <ostream>
#include <panvc3/align.hh>
#include <panvc3/cigar.hh>
#include <panvc3/cigar_adapter.hh>
#include <panvc3/indel_run_checker.hh>
#include <panvc3/msa_index.hh>
#include <panvc3/range.hh>
#include <panvc3/rewrite_cigar.hh>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <vector>


namespace panvc3::alignment_projection::detail {
	
	typedef seqan3::dna5		sequence_alphabet;
	typedef seqan3::phred42		quality_alphabet;
	typedef seqan3::qualified <
		sequence_alphabet,
		quality_alphabet
	>							qualified_alphabet;
	
	
	template <bool t_should_add_qualities, typename t_range>
	auto transformed_reference_range(t_range &&range)
	{
		return range | ranges::views::transform([](auto const cc) {
			sequence_alphabet scc;
			scc.assign_char(cc);
			if constexpr (t_should_add_qualities)
			{
				qualified_alphabet qc;
				qc = scc;
				qc = max_letter <quality_alphabet>();	// Assign maximum qualities to the reference.
				return qc;								// (We could take the variant likelihoods into account, though.)
			}
			else
			{
				return scc;
			}
		});
	}
}


namespace panvc3 {

	struct alignment_projector_base
	{
		virtual ~alignment_projector_base() {}
	};


	struct alignment_projector_delegate
	{
		virtual ~alignment_projector_delegate() {}
		virtual void alignment_projector_begin_realignment(alignment_projector_base const &) = 0;
		virtual void alignment_projector_end_realignment(alignment_projector_base const &) = 0;
	};
	
	
	template <typename t_adapter>
	class alignment_projector_tpl final : public alignment_projector_base
	{
	public:
		typedef t_adapter::vector_type				cigar_vector;
		typedef t_adapter::sequence_alphabet		sequence_alphabet;
		typedef t_adapter::quality_alphabet			quality_alphabet;
		typedef cigar_buffer_tpl <t_adapter>		cigar_buffer;
		typedef indel_run_checker_tpl <t_adapter>	indel_run_checker_type;
		
	private:
		indel_run_checker_type			m_indel_run_checker;
		cigar_vector					m_cigar_realigned;
		cigar_buffer					m_rewrite_buffer;
		cigar_buffer					m_realign_buffer;
		range_vector					m_realigned_reference_ranges;
		range_vector					m_realigned_query_ranges;		// OK b.c. we project one alignment at a time, i.e. there will not be ranges from multiple alignments.
		alignment_projector_delegate	*m_delegate{};
		
	public:
		alignment_projector_tpl() = default;

		alignment_projector_tpl(alignment_projector_delegate &delegate):
			m_delegate(&delegate)
		{
		}

		void reset();
		
		std::size_t project_alignment(
			std::size_t const src_pos,
			msa_index::sequence_entry const &src_seq_entry,
			msa_index::sequence_entry const &dst_seq_entry,
			std::vector <char> const &ref_seq,
			std::vector <sequence_alphabet> const &query_seq,
			cigar_vector const &cigar_seq,
			std::vector <quality_alphabet> const &base_qualities,	// Pass empty vector to not use base qualities.
			std::int32_t const gap_opening_cost,
			std::int32_t const gap_extension_cost
		);
		
		cigar_vector const &alignment() const { return m_cigar_realigned; }
		range_vector const &realigned_reference_ranges() const { return m_realigned_reference_ranges; }
		range_vector const &realigned_query_ranges() const { return m_realigned_query_ranges; }
		
		indel_run_checker_type const &indel_run_checker() const { return m_indel_run_checker; }
	};
	
	
	struct alignment_projector_adapter_libbio : public cigar_adapter_libbio
	{
		typedef char												sequence_alphabet;
		typedef char												quality_alphabet;
		
		typedef alignment_projection::detail::sequence_alphabet		sequence_alphabet_;
		typedef alignment_projection::detail::quality_alphabet		quality_alphabet_;
		typedef alignment_projection::detail::qualified_alphabet	qualified_alphabet_;
		
		typedef std::vector <sequence_alphabet>						sequence_vector;
		typedef std::vector <quality_alphabet>						quality_vector;
		
		template <typename t_value>
		decltype(auto) to_characters(std::vector <t_value> const &vec) const { return vec; }
		
		template <typename t_query_seq>
		auto transformed_query_range(t_query_seq &&query_seq) const;
		
		template <typename t_query_seq, typename t_quality_seq>
		auto transformed_query_range_with_base_qualities(t_query_seq &&query_seq, t_quality_seq &&base_qualities) const;
	};
	
	
	struct alignment_projector_adapter_seqan3 : public cigar_adapter_seqan3
	{
		typedef alignment_projection::detail::sequence_alphabet		sequence_alphabet;
		typedef alignment_projection::detail::quality_alphabet		quality_alphabet;
		typedef alignment_projection::detail::qualified_alphabet	qualified_alphabet;
		
		typedef std::vector <sequence_alphabet>						sequence_vector;
		typedef std::vector <quality_alphabet>						quality_vector;
		
		auto to_characters(sequence_vector const &vec) { return vec | ranges::views::transform([](auto const dna_cc){ return dna_cc.to_char(); }); }
		
		template <typename t_query_seq>
		decltype(auto) transformed_query_range(t_query_seq &&query_seq) const { return std::forward <t_query_seq>(query_seq); }
		
		template <typename t_query_seq, typename t_quality_seq>
		auto transformed_query_range_with_base_qualities(t_query_seq &&query_seq, t_quality_seq &&base_qualities) const;
	};
	
	
	typedef alignment_projector_tpl <alignment_projector_adapter_libbio>	alignment_projector_libbio;
	typedef alignment_projector_tpl <alignment_projector_adapter_seqan3>	alignment_projector_seqan3;
	
	
	template <typename t_adapter>
	void alignment_projector_tpl <t_adapter>::reset()
	{
		// indel_run_checker::reset() called in project_alignment().
		m_cigar_realigned.clear();
		m_rewrite_buffer.clear();
		m_realign_buffer.clear();
		m_realigned_reference_ranges.clear();
		m_realigned_query_ranges.clear();
	}
	
	
	template <typename t_adapter>
	std::size_t alignment_projector_tpl <t_adapter>::project_alignment(
		std::size_t const src_pos,
		msa_index::sequence_entry const &src_seq_entry,
		msa_index::sequence_entry const &dst_seq_entry,
		std::vector <char> const &ref_seq,
		std::vector <sequence_alphabet> const &query_seq,
		cigar_vector const &cigar_seq,
		std::vector <quality_alphabet> const &base_qualities,	// Pass empty vector to not use base qualities.
		std::int32_t const gap_opening_cost,
		std::int32_t const gap_extension_cost
	)
	{
		namespace ap = alignment_projection::detail;
		
		t_adapter adapter{};
		
		auto const dst_pos(rewrite_cigar <t_adapter>(src_pos, cigar_seq, src_seq_entry, dst_seq_entry, adapter.to_characters(query_seq), ref_seq, m_rewrite_buffer));
		
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
			if (m_delegate)
				m_delegate->alignment_projector_begin_realignment(*this);
			
			if (base_qualities.empty())
			{
				align_global <false>(
					ap::transformed_reference_range <false>(ref_seq | ref_range.slice()),
					adapter.transformed_query_range(query_seq | query_range.slice()),
					gap_opening_cost,
					gap_extension_cost,
					m_realign_buffer
				);
			}
			else
			{
				align_global <true, ap::sequence_alphabet, ap::quality_alphabet>(
					ap::transformed_reference_range <true>(ref_seq | ref_range.slice()),
					adapter.transformed_query_range_with_base_qualities(
						query_seq | query_range.slice(),
						base_qualities | query_range.slice()
					),
					gap_opening_cost,
					gap_extension_cost,
					m_realign_buffer
				);
			}
			
			if (m_delegate)
				m_delegate->alignment_projector_end_realignment(*this);
			
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
	
	
	template <typename t_query_seq>
	auto alignment_projector_adapter_libbio::transformed_query_range(t_query_seq &&query_seq) const
	{
		return query_seq | ranges::views::transform([](auto const cc){
			sequence_alphabet_ scc;
			scc.assign_char(cc);
			return scc;
		});
	}
	
	
	template <typename t_query_seq, typename t_quality_seq>
	auto alignment_projector_adapter_libbio::transformed_query_range_with_base_qualities(t_query_seq &&query_seq, t_quality_seq &&base_qualities) const
	{
		sequence_alphabet_ scc;
		quality_alphabet_ qc;
		
		return (
			ranges::views::zip(
				query_seq,
				base_qualities
			)
			| ranges::views::transform([](auto const tup) -> qualified_alphabet_ {
				sequence_alphabet_ scc;
				quality_alphabet_ qc;
				qualified_alphabet_ retval;
				scc.assign_char(std::get <0>(tup));
				qc.assign_char(std::get <1>(tup));
				retval = scc;
				retval = qc;
				return retval;
			})
		);
	}
	
	
	template <typename t_query_seq, typename t_quality_seq>
	auto alignment_projector_adapter_seqan3::transformed_query_range_with_base_qualities(t_query_seq &&query_seq, t_quality_seq &&base_qualities) const
	{
		return (
			ranges::views::zip(
				query_seq,
				base_qualities
			)
			| ranges::views::transform([](auto const tup) -> qualified_alphabet {
				qualified_alphabet retval;
				retval = std::get <0>(tup);
				retval = std::get <1>(tup);
				return retval;
			})
		);
	}
}

#endif
