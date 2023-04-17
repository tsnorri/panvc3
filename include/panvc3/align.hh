/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_ALIGN_HH
#define PANVC3_ALIGN_HH

#include <cmath>
#include <panvc3/cigar.hh>
#include <range/v3/view/zip.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>


namespace panvc3 {
	
	// Formulae from Malde, K. The effect of sequence quality on sequence alignment.
	// https://academic.oup.com/bioinformatics/article/24/7/897/296249
	template <typename t_qual_alphabet>
	struct base_quality_scorer
	{
		typedef double	score_type;
		
		score_type score(bool const matches, t_qual_alphabet const lhs_qual, t_qual_alphabet const rhs_qual) const
		{
			auto const lhs_e(1.0 / std::pow(10.0, lhs_qual.to_phred() / 10.0));
			auto const rhs_e(1.0 / std::pow(10.0, rhs_qual.to_phred() / 10.0));
			auto const combined_e(lhs_e + rhs_e - lhs_e / 3.0 * rhs_e * 4.0);
			if (matches)
				return 2.0 + std::log2(1.0 - combined_e);
			else
				return 2.0 - std::log2(3.0) + std::log2(combined_e);
		}
	};
	
	
	template <typename t_seq_alphabet, typename t_qual_alphabet, typename t_base = base_quality_scorer <t_qual_alphabet>>
	struct base_quality_scoring_scheme : public t_base
	{
		typedef typename t_base::score_type							score_type;
		typedef seqan3::qualified <t_seq_alphabet, t_qual_alphabet>	alphabet_type;
		
		score_type score(alphabet_type const lhs, alphabet_type const rhs) const
		{
			auto const [lhs_cc, lhs_qual] = lhs;
			auto const [rhs_cc, rhs_qual] = rhs;
			return t_base::score(lhs_cc == rhs_cc, lhs_qual, rhs_qual);
		}
	};
	
	template <>
	struct base_quality_scoring_scheme <void, void> {};
	
	
	template <
		bool t_should_use_base_probabilities,
		typename t_sequence_alphabet = void,
		typename t_quality_alphabet = void,
		typename t_seq_1,
		typename t_seq_2
	>
	auto align_global(
		t_seq_1 &&seq1,
		t_seq_2 &&seq2,
		std::int32_t const gap_opening_penalty,
		std::int32_t const gap_extension_penalty,
		cigar_buffer &dst
	)
	{
		using seqan3::operator""_cigar_operation;
		
		dst.clear();
		
		auto const config{
			seqan3::align_cfg::method_global{}
			| seqan3::align_cfg::scoring_scheme{std::conditional_t <
				t_should_use_base_probabilities,
				base_quality_scoring_scheme <t_sequence_alphabet, t_quality_alphabet>,
				seqan3::nucleotide_scoring_scheme <>
			>{}}
			| seqan3::align_cfg::gap_cost_affine{
				seqan3::align_cfg::open_score{gap_opening_penalty},
				seqan3::align_cfg::extension_score{gap_extension_penalty}
		}};
		
		auto res_range(seqan3::align_pairwise(std::tie(seq1, seq2), config));
		auto const &res(*res_range.begin());
		auto const &alignment(res.alignment()); // std::tuple of seqan3::gap_decorators.
		
		for (auto const &ccs : ranges::views::zip(std::get <0>(alignment), std::get <1>(alignment)))
		{
			auto const [ll, rr] = ccs;
			auto const lcc(ll.to_char());
			auto const rcc(rr.to_char());
			
			if ('-' == lcc)
				dst.push_back('I'_cigar_operation, 1);
			else if ('-' == rcc)
				dst.push_back('D'_cigar_operation, 1);
			else if (lcc == rcc)
				dst.push_back('='_cigar_operation, 1);
			else
				dst.push_back('X'_cigar_operation, 1);
		}
		
		dst.finish();
		return res.score();
	}
}

#endif
