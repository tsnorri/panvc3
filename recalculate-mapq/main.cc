/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <condition_variable>
#include <iostream>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/utility.hh> // lb::make_array
#include <mutex>
#include <panvc3/align.hh>
#include <panvc3/sam_tag.hh>
#include <panvc3/utility.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/algorithm/upper_bound.hpp>
#include <range/v3/range/operations.hpp>
#include <range/v3/range/primitives.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/reverse.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <syncstream>
#include <thread>
#include "cmdline.h"

namespace chrono	= std::chrono;
namespace fs		= std::filesystem;
namespace lb		= libbio;
namespace ios		= boost::iostreams;
namespace rsv		= ranges::views;


using seqan3::operator""_tag;


namespace {
	typedef std::uint32_t						chromosome_id_type;
	typedef std::uint32_t						sequence_length_type;
	typedef panvc3::seqan3_sam_tag_type			sam_tag_type;
	
	typedef seqan3::sam_tag_type_t <"AS"_tag>	alignment_score_tag_type;
	typedef double								alignment_score_type;
	typedef std::vector <alignment_score_type>	alignment_score_vector;
	typedef std::uint8_t						mapping_quality_type;
	
	constexpr static inline auto const SEQUENCE_LENGTH_MAX{std::numeric_limits <sequence_length_type>::max()};
	constexpr static inline auto const ALIGNMENT_SCORE_MIN{-std::numeric_limits <alignment_score_type>::max()};
	constexpr static inline mapping_quality_type const MAPQ_NO_NEXT_RECORD{255};
	
	
	template <typename>
	struct is_optional : public std::false_type {};
	
	template <typename t_type>
	struct is_optional <std::optional <t_type>> : public std::true_type {};
	
	template <typename t_type>
	constexpr static inline const bool is_optional_v{is_optional <t_type>::value};
	
	
	template <typename ... t_fns> struct overloaded : t_fns... { using t_fns::operator()...; };
	
	
	struct field_value_out_of_bounds {};
	struct tag_type_mismatch { sam_tag_type	tag{}; };
	struct tag_value_out_of_bounds { sam_tag_type tag{}; };
	
	struct input_error : public std::exception
	{
		typedef std::variant <
			field_value_out_of_bounds,
			tag_type_mismatch,
			tag_value_out_of_bounds
		> error_info_type;
			
		error_info_type	error_info{};
		
		template <typename t_type>
		input_error(t_type &&error_info_):
			error_info(std::forward <t_type>(error_info_))
		{
		}
	};
	
	
	template <typename t_type, typename ... t_args>
	input_error make_input_error(t_args && ... args)
	{
		return input_error{t_type{std::forward <t_args>(args)...}};
	}


	struct sam_tag_specification
	{
		panvc3::seqan3_sam_tag_type ref_n_positions_tag{};
		panvc3::seqan3_sam_tag_type original_reference_tag{};
		panvc3::seqan3_sam_tag_type original_position_tag{};
		panvc3::seqan3_sam_tag_type original_rnext_tag{};
		panvc3::seqan3_sam_tag_type original_pnext_tag{};
		panvc3::seqan3_sam_tag_type original_alignment_score_tag{};
		panvc3::seqan3_sam_tag_type new_alignment_score_tag{};
	};


	struct alignment_scoring
	{
		typedef std::int32_t	score_type;

		score_type	min_mismatch_penalty{};
		score_type	max_mismatch_penalty{};
		score_type	n_penalty{};
		score_type	gap_opening_penalty{};
		score_type	gap_extension_penalty{};
	};
	
	
	template <typename t_aln_rec>
	struct alignment_scorer
	{
		virtual ~alignment_scorer() {}
		virtual alignment_score_type operator()(t_aln_rec &rec, sam_tag_specification const &tag_spec) const = 0;
	};
	
	
	template <typename t_aln_rec>
	struct as_tag_alignment_scorer final : public alignment_scorer <t_aln_rec>
	{
		alignment_score_type operator()(t_aln_rec &rec, sam_tag_specification const &tag_spec) const override;
	};
	
	
	template <typename t_aln_rec>
	class cigar_alignment_scorer final : public alignment_scorer <t_aln_rec>
	{
	public:
		typedef std::int32_t	score_type;

	private:
		constexpr static inline std::size_t const		QUALITY_MAX{256};
		typedef std::array <score_type, QUALITY_MAX>	quality_lookup_array;
		typedef std::vector <std::uint8_t>				ref_n_position_vector;
		
	private:
		quality_lookup_array		m_mismatch_penalties;
		quality_lookup_array		m_n_penalties;
		alignment_scoring			m_scoring{};
		panvc3::seqan3_sam_tag_type	m_ref_n_positions_tag{};
		
	public:
		cigar_alignment_scorer(alignment_scoring &&scoring, panvc3::seqan3_sam_tag_type const ref_n_positions_tag):
			m_scoring(std::move(scoring)),
			m_ref_n_positions_tag(ref_n_positions_tag)
		{
			fill_penalties(m_mismatch_penalties, m_scoring.min_mismatch_penalty, m_scoring.max_mismatch_penalty);
			fill_penalties(m_n_penalties, m_scoring.n_penalty, m_scoring.n_penalty);
		}
		
		alignment_score_type operator()(t_aln_rec &rec, sam_tag_specification const &tag_spec) const override;
		
	private:
		void fill_penalties(quality_lookup_array &penalties, score_type const min, score_type const max);
		alignment_score_type calculate_score(t_aln_rec const &rec, ref_n_position_vector const &ref_n_positions) const;
	};
	
	
	template <typename t_aln_rec>
	alignment_score_type as_tag_alignment_scorer <t_aln_rec>::operator()(t_aln_rec &aln_rec, sam_tag_specification const &) const
	{
		// Retrieve the alignment score from the AS tag.
		auto const &tags(aln_rec.tags());
		auto const it(tags.find("AS"_tag));
		if (tags.end() == it)
			return ALIGNMENT_SCORE_MIN;

		using seqan3::get;
		return get <alignment_score_tag_type>(it->second);
	}


	template <typename t_aln_rec>
	void cigar_alignment_scorer <t_aln_rec>::fill_penalties(quality_lookup_array &penalties, score_type const min, score_type const max)
	{
		for (std::size_t i(0); i < QUALITY_MAX; ++i)
			penalties[i] = std::min(i, std::size_t(40)) / double(40) * (max - min) + min; // Cut-off at 40 like in Bowtie 2.
	}
	
	
	template <typename t_aln_rec>
	alignment_score_type cigar_alignment_scorer <t_aln_rec>::calculate_score(
		t_aln_rec const &aln_rec,
		ref_n_position_vector const &ref_n_positions
	) const
	{
		// Calculate the alignment score from the CIGAR and the base qualities.
		using seqan3::get;
		using seqan3::operator""_cigar_operation;
		using seqan3::operator""_dna5;
		
		auto const &query_seq(aln_rec.sequence());
		auto const &cigar_seq(aln_rec.cigar_sequence());
		auto const &base_qualities(aln_rec.base_qualities());
		std::size_t query_pos{};
		std::size_t ref_pos{};
		char prev_op{};
		alignment_score_type retval{};
		for (auto const &item : cigar_seq)
		{
			auto const count(get <0>(item));
			auto const op(get <1>(item));
			auto const op_(op.to_char());
			
			switch (op_)
			{
				case 'H':	// Hard clipping, consumes nothing.
				case 'P':	// Padding (silent deletion from padded reference), consumes nothing.
					// Do not set prev_op.
					break;
					
				case 'N':	// Skipped region, consumes reference. (In SAMv1, this is only relevant in mRNA-to-genome alignments.)
					// Do not set prev_op.
					ref_pos += count;
					break;
					
				case 'I':	// Insertion, consumes query.
					if (! ('I' == prev_op || 'D' == prev_op))
						retval -= m_scoring.gap_opening_penalty;
					
					retval -= count * m_scoring.gap_extension_penalty;
					query_pos += count;
					prev_op = op_;
					break;
				
				case 'D':	// Deletion, consumes reference.
					if (! ('I' == prev_op || 'D' == prev_op))
						retval -= m_scoring.gap_opening_penalty;
					
					retval -= count * m_scoring.gap_extension_penalty;
					ref_pos += count;
					prev_op = op_;
					break;
				
				case 'S':	// Soft clipping, consumes query.
					query_pos += count;
					prev_op = op_;
					break;

				case '=':	// Match, consumes both.
					for (std::size_t i(0); i < count; ++i)
					{
						if ('N'_dna5 == query_seq[query_pos + i]) // Match implies that reference also has N.
							retval -= m_n_penalties[base_qualities[query_pos + i].to_rank()];
					}
					query_pos += count;
					ref_pos += count;
					prev_op = op_;
					break;
				
				case 'X':	// Mismatch, consumes both.
					for (std::size_t i(0); i < count; ++i)
					{
						auto const qp(query_pos + i);
						auto const rp(ref_pos + i);
						auto const qual(base_qualities[qp].to_rank());
						if ('N'_dna5 == query_seq[qp] || ((0x1 << (rp % 8)) & ref_n_positions[rp / 8]))
							retval -= m_n_penalties[qual];
						else
							retval -= m_mismatch_penalties[qual];
					}
					query_pos += count;
					ref_pos += count;
					prev_op = op_;
					break;
				
				case 'M':	// Match or mismatch, consumes both.
				default:
					libbio_fail("Unexpected CIGAR operation “", op_, "”");
					break;
			}
		}
		
		return retval;
	}
	
	
	template <typename t_aln_rec>
	alignment_score_type cigar_alignment_scorer <t_aln_rec>::operator()(t_aln_rec &aln_rec, sam_tag_specification const &tag_spec) const
	{
		using seqan3::get;

		auto &tags(aln_rec.tags());
		auto const it(tags.find(tag_spec.ref_n_positions_tag));
		if (tags.end() == it)
		{
			std::cerr << "WARNING: Reference N positions not set in record '" << aln_rec.id() << "'; skipping.\n";
			return ALIGNMENT_SCORE_MIN;
		}

		auto const new_score(calculate_score(aln_rec, get <ref_n_position_vector>(it->second)));
		
		// Retrieve the alignment score from the AS tag.
		{
			auto const it(tags.find("AS"_tag));
			if (tags.end() == it)
			{
				auto &as_tag(tags.template get <"AS"_tag>());
				as_tag = new_score;
			}
			else
			{
				tags[tag_spec.original_alignment_score_tag] = get <alignment_score_tag_type>(it->second);
				tags[tag_spec.new_alignment_score_tag] = float(new_score);
				it->second = score_type(new_score);
			}
		}
		
		return new_score;
	}
	
	
	struct mapq_score_calculator
	{
		virtual ~mapq_score_calculator() {}
		virtual mapping_quality_type calculate_mapq(		// of one segment, remember to sum for a pair.
			sequence_length_type const read_length,
			sequence_length_type const other_read_length,	// pass zero if not paired.
			alignment_score_type const score,				// AS 
			alignment_score_type const next_score			// AS/XS, pass ALIGNMENT_SCORE_MIN if no other alignments.
		) const = 0;
	};

	
	// These need to be outside bowtie2_v2_score_calculator so that the definitions of operator() are complete before using in a static assertion.
	struct score_entry
	{
		alignment_score_type	normalised_score_threshold{};
		mapping_quality_type	mapping_quality{};
		
		struct project_key
		{
			constexpr alignment_score_type operator()(score_entry const &entry) const { return entry.normalised_score_threshold; }
		};
	};


	struct score_entry_2
	{
		alignment_score_type	diff_next_threshold{};
		alignment_score_type	normalised_score_threshold{};
		mapping_quality_type	mapping_quality{};
		
		struct project_key
		{
			constexpr auto operator()(score_entry_2 const &entry) const
			{
				return std::make_pair(entry.diff_next_threshold, entry.normalised_score_threshold);
			}
		};
	};


	// Calculates the mapping quality, given the alignment score of an aligned segment (or a pair thereof),
	// the next candidate score, and the segment length(s).
	// Used these as references for calculating MAPQ:
	// Bowtie 2's source
	// https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/531-bowtie-an-ultrafast-memory-efficient-open-source-short-read-aligner/page34?postcount=510#post229040
	// http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html#bt2py
	struct bowtie2_v2_score_calculator final : public mapq_score_calculator
	{
		typedef panvc3::cmp_proj <score_entry,		score_entry::project_key>	score_entry_cmp;
		typedef panvc3::cmp_proj <score_entry_2,	score_entry_2::project_key>	score_entry_2_cmp;

		constexpr static std::array const unique_alignment_scores{
			score_entry{0.0,	0},
			score_entry{0.3,	3},
			score_entry{0.4,	8},
			score_entry{0.5,	23},
			score_entry{0.6,	24},
			score_entry{0.7,	40},
			score_entry{0.8,	42}
		};

		static_assert(
			std::is_sorted(
				unique_alignment_scores.begin(),
				unique_alignment_scores.end(),
				score_entry_cmp{}
			)
		);

		constexpr static std::array const non_unique_alignment_scores{
			score_entry_2{0.0,	0.0,	2},	// For matching positive first component less than 0.1.
			score_entry_2{0.0,	0.67,	6},
			score_entry_2{0.1,	0.0,	0},
			score_entry_2{0.1,	0.67,	7},
			score_entry_2{0.1,	0.88,	12},
			score_entry_2{0.1,	1.0,	30},
			score_entry_2{0.2,	0.0,	0},
			score_entry_2{0.2,	0.67,	11},
			score_entry_2{0.2,	0.88,	17},
			score_entry_2{0.2,	1.0,	31},
			score_entry_2{0.3,	0.0,	3},
			score_entry_2{0.3,	0.67,	15},
			score_entry_2{0.3,	0.88,	18},
			score_entry_2{0.3,	1.0,	32},
			score_entry_2{0.4,	0.0,	4},
			score_entry_2{0.4,	0.68,	14},
			score_entry_2{0.4,	0.84,	21},
			score_entry_2{0.4,	1.0,	34},
			score_entry_2{0.5,	0.0,	5},
			score_entry_2{0.5,	0.68,	16},
			score_entry_2{0.5,	0.84,	25},
			score_entry_2{0.5,	1.0,	35},
			score_entry_2{0.6,	0.0,	22},
			score_entry_2{0.6,	1.0,	36},
			score_entry_2{0.7,	0.0,	26},
			score_entry_2{0.7,	1.0,	37},
			score_entry_2{0.8,	0.0,	27},
			score_entry_2{0.8,	1.0,	38},
			score_entry_2{0.9,	0.0,	33},
			score_entry_2{0.9,	1.0,	39}
		};

		static_assert(
			std::is_sorted(
				non_unique_alignment_scores.begin(),
				non_unique_alignment_scores.end(),
				score_entry_2_cmp{}
			)
		);

		inline alignment_score_type calculate_read_min_score(sequence_length_type const read_length) const;
		inline alignment_score_type calculate_read_max_score(sequence_length_type const read_length) const;

		mapping_quality_type calculate_mapq(				// of one segment, remember to sum for a pair.
			sequence_length_type const read_length,
			sequence_length_type const other_read_length,   // Pass zero if not paired.
			alignment_score_type const score,  			    // AS 
			alignment_score_type const next_score			// AS/XS, pass ALIGNMENT_SCORE_MIN if no other alignments.
		) const override;
	};


	// Bowtie 2's scoring function (version 2) calculates the minimum score from the sum of applying
	//
	//   f(x) = max(I, min(X, C + L * X)) [simple_func.h]
	//
	// where I = -DBL_MAX (and hence unimportant), C = -0.6, L = -0.6 (call to scoreMin.init in
	// bt2_search.cpp, DEFAULT_MIN_CONST and DEFAULT_MIN_LINEAR set to -0.6 in scoring.h), and X
	// is the read length (and hence the right side of the minimum always applies). See also
	// "Valid alignments meet or exceed the minimum score threshold" in MANUAL.markdown.
	auto bowtie2_v2_score_calculator::calculate_read_min_score(sequence_length_type const read_length) const -> alignment_score_type
	{
		if (!read_length)
			return 0;
		return -0.6 + (-0.6 * read_length);
	}

	 
	// In end-to-end mode, the maximum score is always zero.
	// To reaffirm, see Scoring::perfectScore(), Scoring::match() in scoring.h.
	// Scoring passed to new_mapq() in bt2_search.cpp:3050 and 4171, is multiseed_sc in either case.
	// Said instance is initialised on line 4916 (I think). It seems to use the default values,
	// DEFAULT_MATCH_BONUS_TYPE and DEFAULT_MATCH_BONUS, which are COST_MODEL_CONSTANT and zero
	// respectively.
	//
	// We would otherwise multiply read length with a precalculated value plus 0.5 (see scoring.h:304),
	// but since the constant const model with zero constant is used, monotone is set to true and hence
	// Scoring::perfectScore() always returns zero.
	auto bowtie2_v2_score_calculator::calculate_read_max_score(sequence_length_type const read_length) const -> alignment_score_type
	{
		return 0;
	}


	mapping_quality_type bowtie2_v2_score_calculator::calculate_mapq(	// of one segment, remember to sum for a pair.
		sequence_length_type const read_length,
		sequence_length_type const other_read_length,					// Pass zero if not paired.
		alignment_score_type const score,								// AS
		alignment_score_type const next_score_							// AS/XS, pass ALIGNMENT_SCORE_MIN if no other alignments.
	) const
	{
		auto const min_score(calculate_read_min_score(read_length) + calculate_read_min_score(other_read_length)); // Score of a barely valid match
		auto const max_score(calculate_read_max_score(read_length) + calculate_read_max_score(other_read_length)); // Score of a perfect match.
		auto const score_range(std::max(alignment_score_type(1), max_score - min_score));
		
		// Check if the score is too low.
		if (score < min_score)
			return 0;

		// Try to cope with the situation that the next alignment is not good enough.
		auto const next_score(min_score < next_score_ ? next_score_ : min_score);
		
		libbio_assert_lte(score, max_score);
		auto const normalised_score(score - min_score);
		auto const normalised_score_quotient(normalised_score / score_range);
		static_assert(std::is_floating_point_v <alignment_score_type>); // For avoiding UB below.
		libbio_assert_lte(next_score, score);
		auto const diff_next(score - next_score); // bestDiff in Bowtie 2; not undefined behaviour even if next_score == ALIGNMENT_SCORE_MIN.
		auto const diff_next_quotient(diff_next / score_range); // bestdiff / diff in Bowtie 2.
		libbio_assert_lte(0, normalised_score_quotient);
		libbio_assert_lte(normalised_score_quotient, 1.0);

		// Replaced branches with table lookup.
		// See BowtieMapq2::mapq() in Bowtie 2’s unique.h.
		if (ALIGNMENT_SCORE_MIN == next_score) // next_score < min_score in the sample code.
		{
			auto const it(
				std::upper_bound(
					unique_alignment_scores.begin(),
					unique_alignment_scores.end(),
					normalised_score_quotient,
					score_entry_cmp{}
				)
			);
			return (it - 1)->mapping_quality;
		}
		else if (diff_next)
		{
			try
			{
				libbio_assert_lte(0, diff_next_quotient);
				libbio_assert_lte(diff_next_quotient, 1.0);
				
				// The following approach is safer than assuming that the difference of the scores in terms of diff_next_threshold is always 0.1.
				auto const begin(non_unique_alignment_scores.begin());
				auto const it(ranges::upper_bound(
					non_unique_alignment_scores,
					diff_next_quotient,
					ranges::less{},
					[](auto const &entry){ return entry.diff_next_threshold; }
				));
				libbio_assert_neq(begin, it);
				auto const diff_next_threshold((it - 1)->diff_next_threshold);
				auto const it_(std::upper_bound(
					begin,
					it,
					std::make_pair(diff_next_threshold, normalised_score_quotient),
					score_entry_2_cmp{}
				));
				libbio_assert_neq(begin, it_);
				return (it_ - 1)->mapping_quality;
			}
			catch (...)
			{
				std::cerr
					<< "diff_next_quotient: " << diff_next_quotient
					<< " score: " << score
					<< " next_score: " << next_score
					<< " max_score: " << max_score
					<< " min_score: " << min_score
					<< " score_range: " << score_range
					<< std::endl;
				throw;
			}
		}
		else
		{
			if (normalised_score_quotient >= 0.67) return 1;
			return 0;
		}
	}


	template <bool t_nullptr_eq_result, typename t_sequence>
	bool sequences_eq(t_sequence const *lhs, t_sequence const *rhs)
	{
		if (lhs == rhs) // Always true if both are nullptr.
			return true;

		if (!lhs)
			return t_nullptr_eq_result;

		if (!rhs)
			return t_nullptr_eq_result;

		return *lhs == *rhs;
	}
	
	
	struct position
	{
		struct mate_tag {};
		
		chromosome_id_type		chr{};
		sequence_length_type	pos{};
		
		position() = default;
		
		constexpr position(chromosome_id_type const chr_, sequence_length_type const pos_):
			chr(chr_),
			pos(pos_)
		{
		}
		
		template <typename t_reference_id, typename t_position>
		constexpr static inline position from_fields(
			std::optional <t_reference_id> const ref_id,
			std::optional <t_position> const position
		);
		
		template <typename t_aln_record>
		constexpr static inline position from_record_with_tags(
			t_aln_record const &rec,
			sam_tag_type const reference_tag,
			sam_tag_type const position_tag
		);
		
		template <typename t_aln_record>
		constexpr explicit position(t_aln_record const &rec):
			position(from_fields(rec.reference_id(), rec.reference_position()))
		{
		}
		
		template <typename t_aln_record>
		constexpr position(t_aln_record const &rec, mate_tag const):
			position(from_fields(rec.mate_reference_id(), rec.mate_position()))
		{
		}
		
		template <typename t_aln_record>
		constexpr position(
			t_aln_record const &rec,
			sam_tag_type const reference_tag,
			sam_tag_type const position_tag
		):
			position(from_record_with_tags(rec, reference_tag, position_tag))
		{
		}
		
		constexpr auto to_tuple() const { return std::make_tuple(chr, pos); }
		constexpr bool operator<(position const &other) const { return to_tuple() < other.to_tuple(); }
		constexpr bool operator==(position const &other) const { return to_tuple() == other.to_tuple(); }
	};
	
	
	constexpr static inline const position INVALID_POSITION{UINT32_MAX, SEQUENCE_LENGTH_MAX};


	std::ostream &operator<<(std::ostream &os, position const &pos)
	{
		os << pos.chr << ':' << pos.pos;
		return os;
	}
	
	
	template <typename t_reference_id, typename t_position>
	constexpr auto position::from_fields(
		std::optional <t_reference_id> const ref_id,
		std::optional <t_position> const pos
	) -> position
	{
		if (!ref_id)
			return INVALID_POSITION;
		
		if (!pos)
			return INVALID_POSITION;
		
		if (*ref_id < 0) throw make_input_error <field_value_out_of_bounds>();
		if (*pos < 0) throw make_input_error <field_value_out_of_bounds>();
		
		return {chromosome_id_type(*ref_id), sequence_length_type(*pos)};
	}
	
	
	template <typename t_type, typename t_tags, typename t_tag>
	inline std::optional <t_type> find_tag_value(t_tags const &tags, t_tag const tag)
	{
		auto const it(tags.find(tag));
		if (tags.end() == it)
			return {};
		
		auto const *val(std::get_if <t_type>(&it->second));
		if (val)
			return {*val};
		
		throw make_input_error <tag_type_mismatch>(tag);
	}
	
	
	template <typename t_aln_record>
	constexpr auto position::from_record_with_tags(
		t_aln_record const &aln_rec,
		sam_tag_type const reference_tag,
		sam_tag_type const position_tag
	) -> position
	{
		typedef std::remove_cvref_t <decltype(aln_rec.reference_id())>			reference_type;
		typedef std::remove_cvref_t <decltype(aln_rec.reference_position())>	position_type;
		
		// If these fail, handle the situation in find_tag_value.
		static_assert(is_optional_v <reference_type>);
		static_assert(is_optional_v <position_type>);
		
		typedef typename reference_type::value_type	reference_value_type;
		typedef typename position_type::value_type	position_value_type;
		
		auto const &tags(aln_rec.tags());
		auto const ref_id(find_tag_value <reference_value_type>(tags, reference_tag));
		auto const pos(find_tag_value <position_value_type>(tags, position_tag));
		if (! (ref_id && pos))
			return INVALID_POSITION;
		
		if (*ref_id < 0) throw make_input_error <tag_value_out_of_bounds>(reference_tag);
		if (*pos < 0) throw make_input_error <tag_value_out_of_bounds>(position_tag);
		
		return {chromosome_id_type(*ref_id), sequence_length_type(*pos)};
	}
	
	
	struct position_pair
	{
		struct normalised_tag {};
		
		position lhs{};
		position rhs{};
		
		// Again get the generated constructors.
		position_pair() = default;
		
		constexpr position_pair(position const &lhs_, position const &rhs_):
			lhs(lhs_),
			rhs(rhs_)
		{
		}
		
		constexpr position_pair(position const &lhs_, position const &rhs_, normalised_tag const):
			position_pair(lhs_, rhs_)
		{
			this->normalise();
		}
		
		constexpr void normalise() { using std::swap; if (! (lhs < rhs)) swap(lhs, rhs); }
		constexpr auto to_tuple() const { return std::make_tuple(lhs, rhs); }
		constexpr bool operator<(position_pair const &other) const { return to_tuple() < other.to_tuple(); }
		constexpr bool operator==(position_pair const &other) const { return to_tuple() == other.to_tuple(); }
		constexpr bool has_mate() const { return INVALID_POSITION != rhs; }
	};


	std::ostream &operator<<(std::ostream &os, position_pair const &pp)
	{
		os << "(lhs: " << pp.lhs << ", rhs: " << pp.rhs << ")";
		return os;
	}

		
	template <typename t_aln_record>
	struct scored_record
	{
		t_aln_record			*record{};
		alignment_score_type	alignment_score{};
		alignment_score_type	pairwise_score{};
		sequence_length_type	mate_length{};
	};


	template <typename t_aln_record>
	std::ostream &operator<<(std::ostream &os, scored_record <t_aln_record> const &sr)
	{
		os << "(rec: " << sr.record;
		if (sr.record)
			os << " id: " << sr.record->id();
		os << " AS: " << sr.alignment_score << " PS: " << sr.pairwise_score << " mate_length: " << sr.mate_length << ")";
		return os;
	}


	template <typename t_sequence>
	struct paired_segment_score
	{
		typedef t_sequence sequence_type;

		position_pair			positions{};
		sequence_type			*sequence{};
		alignment_score_type	score{};
		alignment_score_type	other_score{};
		bool					has_mate{};
		
		alignment_score_type total_score() const { return score + other_score; }
		alignment_score_type max_score() const { return (has_mate ? std::max(score, other_score) : score); }
		
		struct project_positions
		{
			constexpr auto operator()(paired_segment_score const &desc) const { return desc.positions; }
		};
	};
	
	
	template <typename t_aln_record>
	class mapq_scorer
	{
	public:
		typedef alignment_scorer <t_aln_record>	alignment_scorer_type;
		
	private:
		typedef std::remove_cvref_t <
			decltype(std::declval <t_aln_record>().sequence())
		>										sequence_type;
		typedef std::vector <t_aln_record>		alignment_vector;

		struct statistics
		{
			std::uint64_t	total_alignments{};
			std::uint64_t	unpaired_alignments{};
			std::uint64_t	reads_with_and_without_mate{};
			std::uint64_t	mate_not_found{};
		};
		
		struct segment_description
		{
			struct position			position{};
			alignment_score_type	score{};
			sequence_length_type	length{};
			
			constexpr auto to_position_score_tuple() const { return std::make_tuple(position, score); }
			constexpr bool operator<(segment_description const &other) const { return to_position_score_tuple() < other.to_position_score_tuple(); }
			
			struct project_position
			{
				constexpr auto operator()(segment_description const &desc) const { return desc.position; }
			};
		};
		
		typedef panvc3::cmp_proj <segment_description, typename segment_description::project_position> cmp_segment_description_position;
		typedef scored_record <t_aln_record>			scored_record_type;
		typedef paired_segment_score <sequence_type>	paired_segment_score_type;
		typedef panvc3::cmp_proj <paired_segment_score_type, typename paired_segment_score_type::project_positions> cmp_paired_segment_score_positions;
		
	private:
		alignment_scorer_type					*m_aln_scorer{};
		mapq_score_calculator					*m_scorer{};
		std::vector <segment_description>		m_segment_descriptions_by_original_position;
		std::vector <paired_segment_score_type>	m_paired_segment_scores_by_projected_position;
		std::vector <scored_record_type>		m_scored_records;
		sam_tag_specification					m_sam_tags{};
		struct statistics						m_statistics{};
		
	private:
		void add_paired_segment_score(paired_segment_score_type const &pss);
		
	public:
		mapq_scorer(
			alignment_scorer_type &aln_scorer,
			mapq_score_calculator &scorer,
			sam_tag_specification const &sam_tags
		):
			m_aln_scorer(&aln_scorer),
			m_scorer(&scorer),
			m_sam_tags(sam_tags)
		{
		}
		
		template <typename t_aln_output>
		void process_alignment_group(
			alignment_vector &alignments,
			t_aln_output &&aln_output
		);
			
		struct statistics const &statistics() const { return m_statistics; }
	};


	template <typename t_sequence>
	std::ostream &operator<<(std::ostream &os, paired_segment_score <t_sequence> const &pss)
	{
		os << "(positions: " << pss.positions << " sequence: " << pss.sequence << " score: " << pss.score << " other_score: " << pss.other_score << " has_mate: " << pss.has_mate << ")";
		return os;
	}
	
	
	template <typename t_aln_record>
	void mapq_scorer <t_aln_record>::add_paired_segment_score(paired_segment_score_type const &pss)
	{
		auto it(std::lower_bound(
			m_paired_segment_scores_by_projected_position.begin(),
			m_paired_segment_scores_by_projected_position.end(),
			pss.positions,
			cmp_paired_segment_score_positions{}
		));

		for (; m_paired_segment_scores_by_projected_position.end() != it; ++it)
		{
			if (it->positions != pss.positions)
				break;

			if (it->sequence == pss.sequence)
			{
				if (it->total_score() < pss.total_score())
					*it = pss;
				return;
			}
		}

		m_paired_segment_scores_by_projected_position.emplace(it, pss);
	}
	
	
	template <typename t_aln_record, typename t_error_info>
	void handle_input_error(t_aln_record const &aln_rec, t_error_info const &error_info)
	{
		std::osyncstream cerr(std::cerr);

		auto print_tag([&cerr](auto const tag){
			std::array <char, 2> buffer;
			panvc3::from_tag(tag, buffer);
			ranges::copy(buffer, std::ostream_iterator <char>(cerr));
		});
		
		std::visit(overloaded{
			[&](tag_type_mismatch const &mm){
				cerr << "ERROR: Record ‘" << aln_rec.id() << "’ had unexpected type for tag ‘";
				print_tag(mm.tag);
				cerr << "’.\n";
			},
			[&](tag_value_out_of_bounds const &oob){
				cerr << "ERROR: Value of tag ‘";
				print_tag(oob.tag);
				cerr << "’ out of bounds in record with ID ‘" << aln_rec.id() << "’.\n";
			},
			[&](field_value_out_of_bounds const &){
				cerr << "ERROR: Unexpected field value in record with ID ‘" << aln_rec.id() << "’.\n";
			},
			[&](auto const &){
				cerr << "ERROR: Unexpected error in record with ID ‘" << aln_rec.id() << "’.\n";
			}
		}, error_info);
		cerr << std::flush;
		std::exit(EXIT_FAILURE);
	}
	
	
	// Precondition: alignments contains (all) the alignments with the same identifier but not necessarily the same sequence.
	template <typename t_aln_record>
	template <typename t_aln_output>
	void mapq_scorer <t_aln_record>::process_alignment_group(
		alignment_vector &alignments,
		t_aln_output &&aln_output
	)
	{
		// Bowtie 2 calculates the mapping qualities by considering the best and the second best alignment for each
		// segment in case of single-ended reads (and possibly unpaired alignments?). In case of paired alignments,
		// each pair is given the sum of the alignment scores. (See AlnSinkWrap::selectByScore() and BowtieMapq2::mapq().)
		//
		// Hence, in order to re-calculate the mapping qualities, we need to score each projected, (normalised) position
		// (pair) and  use that information to find the next best alignments s.t. at least the projected position of one
		// segment differs.
		//
		// In SAM, each segment has its own mapping quality even though in Bowtie2, the score is determined from the alignment
		// score of the segment pair. I am not sure if we can assume that the mapping quality of the pair is symmetric,
		// so we calculate the score for each segment. (For instance, there could be segments a and b aligned to position p1,
		// and for both the next segment is located at position p2, where only segment c is aligned; from c’s point of view,
		// the next segment could be at position p1 but now if a and b have different alignment scores, they could also have
		// different mapping qualities.)

		// NOTE: The algorithm works with single-ended and paired-ended reads but not with three or more segments per template.
		// (Not sure if such sequencing technology exists, though.)
		
		if (alignments.empty())
			return;
		
		m_statistics.total_alignments += alignments.size();
		
		// Sort s.t. the records may be searched with the original values of RNEXT and PNEXT.
		// First cache the (original) positions and scores s.t. sorting does not take O(n log^2 n) time.
		m_segment_descriptions_by_original_position.clear();
		m_segment_descriptions_by_original_position.reserve(alignments.size());
		m_scored_records.clear();
		m_scored_records.resize(alignments.size(), scored_record_type{nullptr, ALIGNMENT_SCORE_MIN, ALIGNMENT_SCORE_MIN, 0});
		std::uint8_t seen_record_types{}; // 0x3 for both.
		for (auto &&[aln_rec, scored_rec] : rsv::zip(alignments, m_scored_records))
		{
			try
			{
				bool const has_mate(aln_rec.mate_reference_id() && aln_rec.mate_position());
				seen_record_types |= 0x1 << has_mate;
				m_statistics.unpaired_alignments += (!has_mate);
				
				auto const score((*m_aln_scorer)(aln_rec, m_sam_tags));
				
				scored_rec.record = &aln_rec;
				scored_rec.alignment_score = score;
				
				m_segment_descriptions_by_original_position.emplace_back(
					position{aln_rec, m_sam_tags.original_reference_tag, m_sam_tags.original_position_tag},
					score,
					aln_rec.sequence().size()
				);
			}
			catch (input_error const &exc)
			{
				handle_input_error(aln_rec, exc.error_info);
			}
		}
		
		// Check if there are both alignments with and without a mate (obscure?).
		// We currently handle this differently than Bowtie 2. In other words, we could keep such runs but instead of comparing the scores
		// of all records, we should partition by pairedness. Also the formula in calculate_mapq() should be updated s.t. the scores would
		// be scaled to [0, 1] before calculating the difference.
		if (0x3 == seen_record_types)
		{
			++m_statistics.reads_with_and_without_mate;
			std::osyncstream(std::cerr) << "WARNING: Read ‘" << alignments.front().id() << "’ has both paired and unpaired alignment records; skipping.\n";
			return;
		}
		
		m_segment_descriptions_by_original_position.emplace_back(INVALID_POSITION, 0, 0); // Add a sentinel for extra safety.
		std::sort(m_segment_descriptions_by_original_position.begin(), m_segment_descriptions_by_original_position.end());
		
		// Determine the sum of the scores for each pair.
		// m_scored_records is needed so that we don’t need to re-calculate the score after partitioning the records.
		m_paired_segment_scores_by_projected_position.clear();
		m_paired_segment_scores_by_projected_position.reserve(alignments.size());
		for (auto &sr : m_scored_records)
		{
			try
			{
				auto &aln_rec(*sr.record);
				position_pair const projected_pos{position{aln_rec}, position{aln_rec, typename position::mate_tag{}}, typename position_pair::normalised_tag{}};
				bool const has_mate(projected_pos.has_mate());
				paired_segment_score_type pss{projected_pos, has_mate ? nullptr : &aln_rec.sequence(), sr.alignment_score, 0, false};
				
				// Min. score only allowed for invalid positions.
				libbio_assert_msg((INVALID_POSITION != pss.positions.lhs) ^ (ALIGNMENT_SCORE_MIN == pss.score), "pss: ", pss, " sr: ", sr);
				
				// Determine the mate’s position only if the alignment has a valid position.
				sequence_length_type mate_length{};
				if (INVALID_POSITION != pss.positions.lhs && has_mate)
				{
					position const mate_original_pos{aln_rec, m_sam_tags.original_rnext_tag, m_sam_tags.original_pnext_tag};
					auto const it(std::upper_bound( // Couldn’t get ranges::upper_bound() to work.
						m_segment_descriptions_by_original_position.begin(),
						m_segment_descriptions_by_original_position.end(),
						mate_original_pos,
						cmp_segment_description_position{}
					));
					if (it == m_segment_descriptions_by_original_position.begin())
						++m_statistics.mate_not_found;
					else
					{
						auto const it_(it - 1);
						if (mate_original_pos == it_->position)
						{
							// Store the best score at the mate position if there is one, also store the length of the mate.
							pss.other_score = it_->score;
							pss.has_mate = true;
							mate_length = it_->length;
						}
						else
						{
							++m_statistics.mate_not_found;
						}
					}
				}
				
				sr.pairwise_score = pss.score + pss.other_score;
				sr.mate_length = mate_length;
				
				// Insert if needed.
				add_paired_segment_score(pss);
			}
			catch (input_error const &exc)
			{
				handle_input_error(*sr.record, exc.error_info);
			}
		}
		// Sort by score.
		libbio_assert(!m_paired_segment_scores_by_projected_position.empty());
		ranges::sort(m_paired_segment_scores_by_projected_position, ranges::less{}, [](auto const &pss){ return pss.total_score(); });
		
		// Calculate the mapping qualities.
		for (auto const &sr : m_scored_records)
		{
			libbio_assert(sr.record);
			auto &aln_rec(*sr.record);
			
			try
			{
				auto const &seq(aln_rec.sequence());
				position_pair const pos{position{aln_rec}, position{aln_rec, position::mate_tag{}}, position_pair::normalised_tag{}};
				bool const has_mate(pos.has_mate());

				// Find the alignments with better scores.
				auto const begin(m_paired_segment_scores_by_projected_position.begin());
				auto const end(m_paired_segment_scores_by_projected_position.end());
				auto it(ranges::upper_bound(
					m_paired_segment_scores_by_projected_position,
					sr.pairwise_score,
					ranges::less{},
					[](auto const &pss){ return pss.total_score(); }
				));

				// Linear search only applies if 0x3 == seen_record_types in both cases below (std::find_if() and the for loop).
				bool const is_best_score{
					0x3 == seen_record_types
					? end != std::find_if(it, end, [&seq](auto const &other){ return sequences_eq <true>(&seq, other.sequence); })
					: it == end // If all records are either paired or unpaired, it suffices (initially) to check if there were greater scores.
				};
				for (auto const &other : rsv::reverse(ranges::subrange(begin, it)))
				{
					libbio_assert_neq(INVALID_POSITION, other.positions.lhs);
					
					// Check if there is a segment description that refers to the sequence of the current alignment.
					if (!sequences_eq <true>(it->sequence, other.sequence)) [[unlikely]] // implies 0x1 & seen_record_types.
						continue;
					
					// The sequences match.
					if (other.positions == pos) // The pairs (positions, sequence) are unique, so this must be the entry for aln_rec.
						continue;
					
					// At least one of the positions in the pair does not match (but the scores may match).
					// Bowtie 2 distinguishes mate 1 and 2, which we are not able to do. Hence we choose the better score for the next alignment.
					aln_rec.mapping_quality() = m_scorer->calculate_mapq(seq.size(), sr.mate_length, sr.pairwise_score, has_mate ? other.total_score() : other.max_score());
					goto output_record;
				}

				if (is_best_score)
					aln_rec.mapping_quality() = m_scorer->calculate_mapq(seq.size(), sr.mate_length, sr.pairwise_score, ALIGNMENT_SCORE_MIN);
				else // There was a better score, hence the current record had the worst score.
					aln_rec.mapping_quality() = MAPQ_NO_NEXT_RECORD;

			output_record:
				aln_rec.header_ptr() = &aln_output.header();
				aln_output.push_back(aln_rec);
			}
			catch (...)
			{
				std::osyncstream(std::cerr) << "*** Read ID ‘" << aln_rec.id() << "’\n" << std::flush;
				throw;
			}
		}
	}


	template <typename t_aln_input, typename t_aln_output, typename t_record = typename std::remove_cvref_t <t_aln_input>::record_type>
	void process_(
		t_aln_input &&aln_input,
		t_aln_output &&aln_output,
		mapq_scorer <t_record> &scorer,
		std::uint16_t const status_output_interval
	)
	{
		typedef std::remove_cvref_t <t_aln_input>	input_type;
		typedef t_record							record_type;
		typedef chrono::steady_clock				clock_type;
		
		std::vector <record_type> rec_buffer;
		std::uint64_t rec_idx{};

		auto const start_time(clock_type::now());
		panvc3::timer status_output_timer;
		auto status_output_thread{[&]() -> std::jthread {
			if (!status_output_interval)
				return {};

			return std::jthread{[&](){
				while (status_output_timer.wait_for(chrono::minutes(status_output_interval)))
				{
					auto const pp(clock_type::now());
					auto const running_time{pp - start_time};

					std::osyncstream cerr(std::cerr);
					lb::log_time(cerr) << "Time spent processing: ";
					panvc3::log_duration(cerr, running_time);
					cerr << "; processed " << (rec_idx - 1) << " records";

					if (rec_idx)
					{
						double usecs_per_record(chrono::duration_cast <chrono::microseconds>(running_time).count());
						usecs_per_record /= rec_idx;
						cerr << " (in " << usecs_per_record << " µs / record)";
					}

					cerr << ".\n" << std::flush;
				}
			}};
		}()};

		for (auto &aln_rec : aln_input)
		{
			if (rec_idx && 0 == rec_idx % 10'000'000)
			{
				std::osyncstream cerr(std::cerr);
				lb::log_time(cerr) << "Processed " << rec_idx << " alignments…\n" << std::flush;
			}
			++rec_idx;

			if (rec_buffer.empty())
			{
				rec_buffer.emplace_back(std::move(aln_rec));
				continue;
			}
			
			auto const &id(aln_rec.id());
			auto const &eq_class_id(rec_buffer.front().id());
			
			if (id != eq_class_id)
			{
				scorer.process_alignment_group(rec_buffer, aln_output);
				rec_buffer.clear();
			}
			
			rec_buffer.emplace_back(std::move(aln_rec));
		}
		
		if (!rec_buffer.empty())
			scorer.process_alignment_group(rec_buffer, aln_output);

		status_output_timer.stop();
		if (status_output_thread.joinable())
			status_output_thread.join();

		{
			std::osyncstream cerr(std::cerr);
			lb::log_time(cerr) << "Done.\n";
			auto const &statistics(scorer.statistics());
			cerr << "\tTotal alignments: " << statistics.total_alignments << '\n';
			cerr << "\tUnpaired alignments: " << statistics.unpaired_alignments << '\n';
			cerr << "\tRecords with mate missing: " << statistics.mate_not_found << '\n';
			cerr << "\tReads with and without a mate: " << statistics.reads_with_and_without_mate << '\n';
			cerr << std::flush;
		}
	}


	template <typename t_header, typename t_header_>
	void copy_program_info(t_header const &input_header, t_header_ &output_header)
	{
		output_header.program_infos.resize(input_header.program_infos.size());
		// FIXME: Not a particularly good idea; filed a feature request to make program_info_t not depend on the template parameter.
		for (auto &&[src, dst] : rsv::zip(input_header.program_infos, output_header.program_infos))
			dst = *reinterpret_cast <t_header_::program_info_t const *>(&src);
	}


	template <typename t_header, typename t_header_>
	void append_program_info(t_header const &input_header, t_header_ &output_header, int const argc, char const * const * const argv)
	{
		typedef typename t_header_::program_info_t program_info_type;
		//copy_program_info(input_header, output_header);
		panvc3::append_sam_program_info(
			"panvc3.recalculate-mapq.",
			"PanVC 3 recalculate_mapq",
			argc,
			argv,
			CMDLINE_PARSER_VERSION,
			output_header.program_infos
		);
	}
	
	
	void process(gengetopt_args_info const &args_info, int const argc, char const * const * const argv)
	{
		// Status output interval
		if (args_info.status_output_interval_arg < 0)
		{
			std::cerr << "ERROR: Status output interval must be non-negative.\n";
			std::exit(EXIT_FAILURE);
		}

		// Mismatch penalties.
		if (args_info.min_mismatch_penalty_arg < 0)		{ std::cerr << "ERROR: Minimum mismatch penalty must be non-negative.\n";	std::exit(EXIT_FAILURE); }
		if (args_info.max_mismatch_penalty_arg < 0)		{ std::cerr << "ERROR: Maximum mismatch penalty must be non-negative.\n";	std::exit(EXIT_FAILURE); }
		if (args_info.n_penalty_arg < 0)				{ std::cerr << "ERROR: N character penalty must be non-negative.\n";		std::exit(EXIT_FAILURE); }
		if (args_info.gap_opening_penalty_arg < 0)		{ std::cerr << "ERROR: Gap opening penalty must be non-negative.\n";		std::exit(EXIT_FAILURE); }
		if (args_info.gap_extension_penalty_arg < 0)	{ std::cerr << "ERROR: Gap extension penalty must be non-negative.\n";		std::exit(EXIT_FAILURE); }

		// Open the SAM input and output.
		typedef seqan3::sam_file_input <>								input_type;
		typedef typename input_type::record_type						record_type;
		auto aln_input{[&](){
			auto const make_input_type([&]<typename t_fmt>(t_fmt &&fmt){
				if (args_info.alignments_arg)
					return input_type(fs::path{args_info.alignments_arg});
				else
					return input_type(std::cin, std::forward <t_fmt>(fmt));
			});
			
			if (args_info.bam_input_flag)
				return make_input_type(seqan3::format_bam{});
			else
				return make_input_type(seqan3::format_sam{});
		}()};

		auto const &aln_input_header(aln_input.header());
		auto const &input_ref_ids(aln_input.header().ref_ids()); // ref_ids() not const.
		typedef std::remove_cvref_t <decltype(input_ref_ids)>			ref_ids_type;
		ref_ids_type output_ref_ids;
		typedef seqan3::sam_file_output <
			typename input_type::selected_field_ids,
			seqan3::type_list <seqan3::format_sam, seqan3::format_bam>,
			ref_ids_type
		>																output_type;
		typedef mapq_scorer <record_type>								mapq_scorer_type;
		
		auto aln_output{[&](){
			// Make sure that aln_output has some header information.
			auto const make_output_type([&]<typename t_fmt>(t_fmt &&fmt){
				ref_ids_type empty_ref_ids;
				if (args_info.output_path_arg)
				{
					return output_type(
						fs::path{args_info.output_path_arg},
						output_ref_ids,
						rsv::empty <std::size_t>()	// Reference lengths; the constructor expects a forward range.
					);
				}
				else
				{
					return output_type(
						std::cout,
						output_ref_ids,
						rsv::empty <std::size_t>(),	// Reference lengths; the constructor expects a forward range.
						std::forward <t_fmt>(fmt)
					);
				}
			});
			
			if (args_info.output_bam_flag)
				return make_output_type(seqan3::format_bam{});
			else
				return make_output_type(seqan3::format_sam{});
		}()};
		
		// Now that we have valid header information in aln_output, we can access it and copy the headers from aln_input.
		{
			output_ref_ids = input_ref_ids;
			
			auto &aln_output_header(aln_output.header());
			aln_output_header.sorting = aln_input_header.sorting;
			aln_output_header.subsorting = aln_input_header.subsorting;
			aln_output_header.grouping = aln_input_header.subsorting;
			aln_output_header.program_infos = aln_input_header.program_infos;
			aln_output_header.comments = aln_input_header.comments;
			aln_output_header.ref_id_info = aln_input_header.ref_id_info;
			aln_output_header.ref_dict = aln_input_header.ref_dict;
			aln_output_header.read_groups = aln_input_header.read_groups;
			
			append_program_info(aln_input_header, aln_output_header, argc, argv);
		}
		
		auto const make_sam_tag([](char const *tag, bool const should_allow_any = false) -> sam_tag_type {
			if (!tag)
				return 0;
			
			std::array <char, 2> buffer{tag[0], tag[1]};
			return panvc3::to_tag(buffer);
		});
		
		sam_tag_specification const sam_tags{
			.ref_n_positions_tag{make_sam_tag(args_info.ref_n_positions_tag_arg)},
			.original_reference_tag{make_sam_tag(args_info.original_rname_tag_arg)},
			.original_position_tag{make_sam_tag(args_info.original_pos_tag_arg)},
			.original_rnext_tag{make_sam_tag(args_info.original_rnext_tag_arg)},
			.original_pnext_tag{make_sam_tag(args_info.original_pnext_tag_arg)},
			.original_alignment_score_tag{make_sam_tag(args_info.original_alignment_score_tag_arg)},
			.new_alignment_score_tag{make_sam_tag(args_info.new_alignment_score_tag_arg)}
		};
		
		if (args_info.print_reference_names_flag)
		{
			auto const &ref_ids(aln_input.header().ref_ids());
			std::cerr << "Reference IDs:\n";
			for (auto const &[idx, name] : rsv::enumerate(ref_ids))
				std::cerr << idx << '\t' << name << '\n';
		}

		lb::log_time(std::cerr) << "Processing the alignments…\n";
		bowtie2_v2_score_calculator score_calculator;
		if (args_info.rescore_alignments_given)
		{
			// FIXME: Check upper bound (INT32_MAX).
			cigar_alignment_scorer <record_type> aln_scorer(alignment_scoring{
				.min_mismatch_penalty	= alignment_scoring::score_type(args_info.min_mismatch_penalty_arg),
				.max_mismatch_penalty	= alignment_scoring::score_type(args_info.max_mismatch_penalty_arg),
				.n_penalty				= alignment_scoring::score_type(args_info.n_penalty_arg),
				.gap_opening_penalty	= alignment_scoring::score_type(args_info.gap_opening_penalty_arg),
				.gap_extension_penalty	= alignment_scoring::score_type(args_info.gap_extension_penalty_arg)
			}, sam_tags.ref_n_positions_tag);
			mapq_scorer_type scorer(aln_scorer, score_calculator, sam_tags);
			process_(aln_input, aln_output, scorer, args_info.status_output_interval_arg);
		}
		else
		{
			as_tag_alignment_scorer <record_type> aln_scorer;
			mapq_scorer_type scorer(aln_scorer, score_calculator, sam_tags);
			process_(aln_input, aln_output, scorer, args_info.status_output_interval_arg);
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	//lb::setup_allocated_memory_logging();
	
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
	
	process(args_info, argc, argv);
	
	return EXIT_SUCCESS;
}
