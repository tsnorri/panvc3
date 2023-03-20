/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <iostream>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/utility.hh> // lb::make_array
#include <panvc3/sam_tag.hh>
#include <panvc3/utility.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/algorithm/upper_bound.hpp>
#include <range/v3/range/operations.hpp>
#include <range/v3/range/primitives.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/reverse.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include "cmdline.h"

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace ios	= boost::iostreams;
namespace rsv	= ranges::views;


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
		panvc3::seqan3_sam_tag_type original_reference_tag{};
		panvc3::seqan3_sam_tag_type original_position_tag{};
		panvc3::seqan3_sam_tag_type original_rnext_tag{};
		panvc3::seqan3_sam_tag_type original_pnext_tag{};
	};


	struct alignment_scorer
	{
		virtual ~alignment_scorer() {}
		virtual mapping_quality_type calculate_mapq(		// of one segment, remember to sum for a pair.
			sequence_length_type const read_length,
			sequence_length_type const other_read_length,	// pass zero if not paired.
			alignment_score_type const score,				// AS 
			alignment_score_type const next_score			// AS/XS, pass ALIGNMENT_SCORE_MIN if no other alignments.
		) const = 0;
	};

	
	// These need to be outside bowtie2_v2_scorer so that the definitions of operator() are complete before using in a static assertion.
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
	struct bowtie2_v2_scorer final : public alignment_scorer
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
	auto bowtie2_v2_scorer::calculate_read_min_score(sequence_length_type const read_length) const -> alignment_score_type
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
	auto bowtie2_v2_scorer::calculate_read_max_score(sequence_length_type const read_length) const -> alignment_score_type
	{
		return 0;
	}


	mapping_quality_type bowtie2_v2_scorer::calculate_mapq(	// of one segment, remember to sum for a pair.
		sequence_length_type const read_length,
		sequence_length_type const other_read_length,		// Pass zero if not paired.
		alignment_score_type const score,					// AS
		alignment_score_type const next_score				// AS/XS, pass ALIGNMENT_SCORE_MIN if no other alignments.
	) const
	{
		auto const min_score(calculate_read_min_score(read_length) + calculate_read_min_score(other_read_length)); // Score of a barely valid match
		auto const max_score(calculate_read_max_score(read_length) + calculate_read_max_score(other_read_length)); // Score of a perfect match.
		auto const score_range(std::max(alignment_score_type(1), max_score - min_score));

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
	};
	
	
	template <typename t_aln_record>
	class mapq_scorer
	{
	private:
		typedef std::remove_cvref_t <
			decltype(std::declval <t_aln_record>().sequence())
		>										sequence_type;
		typedef std::vector <t_aln_record>		alignment_vector;

		struct statistics
		{
			std::uint64_t	unpaired_alignments{};
			std::uint64_t	reads_with_and_without_mate{};
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
		
		struct paired_segment_score
		{
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
		
		typedef panvc3::cmp_proj <paired_segment_score, typename paired_segment_score::project_positions> cmp_paired_segment_score_positions;
		
		struct scored_record
		{
			t_aln_record			*record{};
			alignment_score_type	score{};
			sequence_length_type	mate_length{};
		};
		
	private:
		alignment_scorer					*m_scorer{};
		std::vector <segment_description>	m_segment_descriptions_by_original_position;
		std::vector <paired_segment_score>	m_paired_segment_scores_by_projected_position;
		std::vector <scored_record>			m_scored_records;
		sam_tag_specification				m_sam_tags{};
		struct statistics					m_statistics{};
		
	private:
		alignment_score_type alignment_score(t_aln_record const &aln_rec) const;
		void add_paired_segment_score(paired_segment_score const &pss);
		
	public:
		mapq_scorer(alignment_scorer &scorer, sam_tag_specification const &sam_tags):
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
	
	
	template <typename t_aln_record>
	alignment_score_type mapq_scorer <t_aln_record>::alignment_score(t_aln_record const &aln_rec) const
	{
		auto const &tags(aln_rec.tags());
		auto const it(tags.find("AS"_tag));
		if (tags.end() == it)
			return ALIGNMENT_SCORE_MIN;

		return std::get <alignment_score_tag_type>(it->second);
	}


	template <typename t_aln_record>
	void mapq_scorer <t_aln_record>::add_paired_segment_score(paired_segment_score const &pss)
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
		auto print_tag([](auto const tag){
			std::array <char, 2> buffer;
			panvc3::from_tag(tag, buffer);
			ranges::copy(buffer, std::ostream_iterator <char>(std::cerr));
		});
		
		std::visit(overloaded{
			[&](tag_type_mismatch const &mm){
				std::cerr << "ERROR: Record ‘" << aln_rec.id() << "’ had unexpected type for tag ‘";
				print_tag(mm.tag);
				std::cerr << "’.\n";
			},
			[&](tag_value_out_of_bounds const &oob){
				std::cerr << "ERROR: Value of tag ‘";
				print_tag(oob.tag);
				std::cerr << "’ out of bounds in record with ID ‘" << aln_rec.id() << "’.\n";
			},
			[&](field_value_out_of_bounds const &){
				std::cerr << "ERROR: Unexpected field value in record with ID ‘" << aln_rec.id() << "’.\n";
			},
			[&](auto const &){
				std::cerr << "ERROR: Unexpected error in record with ID ‘" << aln_rec.id() << "’.\n";
			}
		}, error_info);
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
		
		// Sort s.t. the records may be searched with the original values of RNEXT and PNEXT.
		// First cache the (original) positions and scores s.t. sorting does not take O(n log^2 n) time.
		m_segment_descriptions_by_original_position.clear();
		m_segment_descriptions_by_original_position.reserve(alignments.size());
		for (auto const &aln_rec : alignments)
		{
			try
			{
				m_segment_descriptions_by_original_position.emplace_back(
					position{aln_rec, m_sam_tags.original_reference_tag, m_sam_tags.original_position_tag},
					alignment_score(aln_rec),
					aln_rec.sequence().size()
				);
			}
			catch (input_error const &exc)
			{
				handle_input_error(aln_rec, exc.error_info);
			}
		}
		m_segment_descriptions_by_original_position.emplace_back(INVALID_POSITION, 0, 0); // Add a sentinel to make the mate searching below slightly simpler.
		std::sort(m_segment_descriptions_by_original_position.begin(), m_segment_descriptions_by_original_position.end());
		
		// Determine the sum of the scores for each pair.
		// m_scored_records is needed so that we don’t need to re-calculate the score after partitioning the records.
		m_paired_segment_scores_by_projected_position.clear();
		m_paired_segment_scores_by_projected_position.reserve(alignments.size());
		m_scored_records.clear();
		m_scored_records.resize(alignments.size(), scored_record{nullptr, ALIGNMENT_SCORE_MIN, 0});
		std::uint8_t seen_record_types{}; // 0x3 for both.
		for (auto &&[aln_rec, scored_rec] : rsv::zip(alignments, m_scored_records))	
		{
			try
			{
				bool const has_mate(aln_rec.mate_reference_id() && aln_rec.mate_position());
				seen_record_types |= 0x1 << has_mate;
				m_statistics.unpaired_alignments += (!has_mate);
			
				position_pair const projected_pos{position{aln_rec}, position{aln_rec, typename position::mate_tag{}}, typename position_pair::normalised_tag{}};
				paired_segment_score pss{projected_pos, has_mate ? nullptr : &aln_rec.sequence(), alignment_score(aln_rec), 0, false};
			
				// Min. score only allowed for invalid positions.
				libbio_assert((INVALID_POSITION != pss.positions.lhs) ^ (ALIGNMENT_SCORE_MIN == pss.score));
			
				// Determine the mate’s position only if the alignment has a valid position.
				sequence_length_type mate_length{};
				if (INVALID_POSITION != pss.positions.lhs)
				{
					position const mate_original_pos{aln_rec, m_sam_tags.original_rnext_tag, m_sam_tags.original_pnext_tag};
					auto const it(std::upper_bound( // Couldn’t get ranges::upper_bound() to work.
						m_segment_descriptions_by_original_position.begin(),
						m_segment_descriptions_by_original_position.end(),
						mate_original_pos,
						cmp_segment_description_position{}
					));
					if (it != m_segment_descriptions_by_original_position.begin())
					{
						auto const it_(it - 1);
						if (mate_original_pos == it_->position)
						{
							// Store the best score at the mate position if there is one, also store the length of the mate.
							pss.other_score = it_->score;
							pss.has_mate = true;
							mate_length = it_->length;
						}
					}
				}
				
				scored_rec = scored_record(&aln_rec, pss.score + pss.other_score, mate_length);
				
				// Insert if needed.
				add_paired_segment_score(pss);
			}
			catch (input_error const &exc)
			{
				handle_input_error(aln_rec, exc.error_info);
			}
		}
		// Sort by score.
		libbio_assert(!m_paired_segment_scores_by_projected_position.empty());
		ranges::sort(m_paired_segment_scores_by_projected_position, ranges::less{}, [](auto const &pss){ return pss.total_score(); });
		
		// Check if there are both alignments with and without a mate (obscure?).
		if (0x3 == seen_record_types)
			++m_statistics.reads_with_and_without_mate;
		
		// Calculate the mapping qualities.
		for (auto const &sr : m_scored_records)
		{
			libbio_assert(sr.record);
			auto &aln_rec(*sr.record);
			
			try
			{
				auto const &seq(aln_rec.sequence());
				auto const seq_len(seq.size());
				position_pair const pos{position{aln_rec}, position{aln_rec, position::mate_tag{}}, position_pair::normalised_tag{}};

				// Find the alignments with better scores.
				auto const begin(m_paired_segment_scores_by_projected_position.begin());
				auto const end(m_paired_segment_scores_by_projected_position.end());
				auto it(ranges::upper_bound(
					m_paired_segment_scores_by_projected_position,
					sr.score,
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
					if (0x3 == seen_record_types && !sequences_eq <true>(it->sequence, other.sequence)) [[unlikely]]
						continue;
					
					// The sequences match.
					if (other.positions == pos) // The pairs (positions, sequence) are unique, so this must be the entry for aln_rec.
						continue;
					
					// At least one of the positions in the pair does not match (but the scores may match).
					// Bowtie 2 distinguishes mate 1 and 2, which we are not able to do. Hence we choose the better score for the next alignment.
					aln_rec.mapping_quality() = m_scorer->calculate_mapq(seq_len, sr.mate_length, sr.score, seq_len ? other.total_score() : other.max_score());
					goto output_record;
				}

				if (is_best_score)
					aln_rec.mapping_quality() = m_scorer->calculate_mapq(seq_len, sr.mate_length, sr.score, ALIGNMENT_SCORE_MIN);
				else // There was a better score, hence the current record had the worst score.
					aln_rec.mapping_quality() = MAPQ_NO_NEXT_RECORD;

			output_record:	
				aln_output.push_back(aln_rec);
			}
			catch (...)
			{
				std::cerr << "*** Read ID ‘" << aln_rec.id() << "’\n" << std::flush;
				throw;
			}
		}
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	void process_(t_aln_input &&aln_input, t_aln_output &&aln_output, sam_tag_specification const &sam_tags)
	{
		typedef std::remove_cvref_t <t_aln_input>	input_type;
		typedef typename input_type::record_type	record_type;
		
		bowtie2_v2_scorer aln_scorer;
		mapq_scorer <record_type> scorer(aln_scorer, sam_tags);
		
		std::vector <record_type> rec_buffer;
		std::uint64_t rec_idx{};
		for (auto &aln_rec : aln_input)
		{
			++rec_idx;
			if (0 == (rec_idx % 10'000'000 - 1))
				lb::log_time(std::cerr) << "Processed " << (1 + rec_idx) << " alignments…\n";
			
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

		lb::log_time(std::cerr) << "Done.\n";
		auto const &statistics(scorer.statistics());
		std::cerr << "\tUnpaired alignments: " << statistics.unpaired_alignments << '\n';
		std::cerr << "\tReads with and without a mate: " << statistics.reads_with_and_without_mate << '\n';
	}
	
	
	void process(gengetopt_args_info const &args_info)
	{
		// Open the SAM input and output.
		typedef seqan3::sam_file_input <>								input_type;
		typedef seqan3::sam_file_output <
			typename input_type::selected_field_ids,
			seqan3::type_list <seqan3::format_sam, seqan3::format_bam>
		>																output_type;
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
		
		auto aln_output{[&](){
			// Make sure that aln_output has some header information.
			auto const make_output_type([&]<typename t_fmt>(t_fmt &&fmt){
				if (args_info.output_path_arg)
					return output_type(fs::path{args_info.output_path_arg});
				else
					return output_type(std::cout, std::forward <t_fmt>(fmt));
			});
			
			if (args_info.output_bam_flag)
				return make_output_type(seqan3::format_bam{});
			else
				return make_output_type(seqan3::format_sam{});
		}()};
		
		auto const make_sam_tag([](char const *tag, bool const should_allow_any = false) -> sam_tag_type {
			if (!tag)
				return 0;
			
			std::array <char, 2> buffer{tag[0], tag[1]};
			return panvc3::to_tag(buffer);
		});
		
		sam_tag_specification const sam_tags{
			.original_reference_tag{make_sam_tag(args_info.original_rname_tag_arg)},
			.original_position_tag{make_sam_tag(args_info.original_pos_tag_arg)},
			.original_rnext_tag{make_sam_tag(args_info.original_rnext_tag_arg)},
			.original_pnext_tag{make_sam_tag(args_info.original_pnext_tag_arg)}
		};
		
		if (args_info.print_reference_names_flag)
		{
			auto const &ref_ids(aln_input.header().ref_ids());
			std::cerr << "Reference IDs:\n";
			for (auto const &[idx, name] : rsv::enumerate(ref_ids))
				std::cerr << idx << '\t' << name << '\n';
		}
		
		lb::log_time(std::cerr) << "Processing the alignments…\n";
		process_(aln_input, aln_output, sam_tags);
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
	
	process(args_info);
	
	return EXIT_SUCCESS;
}
