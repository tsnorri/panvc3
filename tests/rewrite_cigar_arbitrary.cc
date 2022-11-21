/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>								// std::iter_swap
#include <boost/pfr/core.hpp>						// boost::pfr::structure_to_tuple
#include <catch2/catch.hpp>
#include <libbio/assert.hh>
#include <panvc3/fmap.hh>
#include <panvc3/rewrite_cigar.hh>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/iterator/stream_iterators.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>						// rc::prop
#include "rapidcheck_additions.hh"					// rc::gen::inClosedRange


namespace gen	= Catch::Generators;
namespace lb	= libbio;
namespace rsv	= ranges::views;

using seqan3::operator""_cigar_operation;


namespace {
	
	template <typename t_type>
	using pair = std::tuple <t_type, t_type>;
	
	template <typename t_type>
	using triplet = std::tuple <t_type, t_type, t_type>;
	
	
	template <typename t_range>
	auto to_vector(t_range &&rng)
	{
		return std::vector(ranges::begin(rng), ranges::end(rng));
	}
	
	
	// Helper for getting the tuple type instead of an instance.
	template <typename t_type>
	using structure_to_tuple_t = decltype(boost::pfr::structure_to_tuple(std::declval <t_type>()));
	
	
	// We could just use SeqAn’s dna4 alphabet but this way debugging will be easier,
	// since each enum value is a readable character.
	struct nucleotide
	{
		enum class value_type : char
		{
			A = 'A',
			C = 'C',
			G = 'G',
			T = 'T'
		};
		
		constexpr static inline auto values{std::to_array({
			value_type::A,
			value_type::C,
			value_type::G,
			value_type::T
		})};
		
		value_type value{value_type::A};
		
		constexpr static value_type to_nucleotide_value(char const cc)
		{
			switch (cc)
			{
				case 'A': return value_type::A;
				case 'C': return value_type::C;
				case 'G': return value_type::G;
				case 'T': return value_type::T;
				
				default:
					libbio_fail("Unexpected value");
			}
		}
		
		constexpr nucleotide() = default;
		
		constexpr /* implicit */ nucleotide(value_type const value_):
			value(value_)
		{
		}
		
		constexpr explicit nucleotide(char const value_):
			value(to_nucleotide_value(value_))
		{
		}
		
		// For getting a different nucleotide.
		constexpr nucleotide next_nucleotide() const
		{
			switch (value)
			{
				case value_type::A: return nucleotide(value_type::C);
				case value_type::C: return nucleotide(value_type::G);
				case value_type::G: return nucleotide(value_type::T);
				case value_type::T: return nucleotide(value_type::A);
			}
			
			// Silence a compiler error.
			libbio_fail("Should not be reached");
			return nucleotide(value_type::A);
		}
		
		constexpr char to_char() const { return std::to_underlying(value); }
		constexpr /* implicit */ operator char() const { return to_char(); }
	};
	
	static_assert(1 == sizeof(nucleotide));
	
	
	// For constructing MSAs.
	struct gapped_nucleotide
	{
		enum class value_type : char
		{
			A	= 'A',
			C	= 'C',
			G	= 'G',
			T	= 'T',
			gap	= '-'
		};
		
		constexpr static inline auto values{std::to_array({
			value_type::A,
			value_type::C,
			value_type::G,
			value_type::T,
			value_type::gap
		})};
		
		value_type value{value_type::gap};
		
		constexpr static value_type to_nucleotide_value(char const cc)
		{
			switch (cc)
			{
				case 'A': return value_type::A;
				case 'C': return value_type::C;
				case 'G': return value_type::G;
				case 'T': return value_type::T;
				case '-': return value_type::gap;
				
				default:
					libbio_fail("Unexpected value");
			}
		}
		
		constexpr gapped_nucleotide() = default;
		
		constexpr /* implicit */ gapped_nucleotide(value_type const value_):
			value(value_)
		{
		}
		
		constexpr explicit gapped_nucleotide(char const value_):
			value(to_nucleotide_value(value_))
		{
		}
		
		constexpr /* implicit */ gapped_nucleotide(nucleotide const nt):
			value(static_cast <value_type>(nt.value))
		{
		}
		
		constexpr char to_char() const { return std::to_underlying(value); }
		constexpr /* implicit */ operator char() const { return to_char(); }
	};
	
	static_assert(1 == sizeof(gapped_nucleotide));
	
	
	constexpr gapped_nucleotide operator ""_gapped_nucleotide(char const cc) { return gapped_nucleotide(cc); }
	
	
	// For constructing sequences of any enum.
	template <typename t_enum, t_enum t_max, t_enum t_min = static_cast <t_enum>(0)>
	struct enum_wrapper
	{
		typedef t_enum								value_type;
		typedef std::underlying_type_t <value_type>	underlying_type;
		
		constexpr static inline underlying_type const min{std::to_underlying(t_min)};
		constexpr static inline underlying_type const max{std::to_underlying(t_max)};
		
		value_type value{t_min};
		
		constexpr static value_type to_value_type(underlying_type const val)
		{
			// Works with unsigned underlying_types.
			if (! (min <= val && val <= max))
				throw std::invalid_argument("Unexpected value");
			
			return static_cast <t_enum>(val);
		}
		
		constexpr enum_wrapper() = default;
		
		constexpr /* implicit */ enum_wrapper(value_type const value_):
			value(value_)
		{
		}
		
		constexpr explicit enum_wrapper(underlying_type const value_):
			value(to_value_type(value_))
		{
		}
		
		constexpr underlying_type to_underlying() const { return std::to_underlying(value); }
		constexpr /* implicit */ operator underlying_type() const { return to_underlying(); }
		constexpr bool operator==(value_type const value_) const { return value == value_; }
	};
	
	
	// Sequence and segment output.
	enum class alignment_output : std::uint8_t
	{ 
		unaligned,
		align_sequences,
		align_all
	};
	
	
	// Part of the MSA that does not contain the aligned segment.
	struct msa_non_segment_section
	{
		enum class operation_type : std::uint8_t
		{
			match = 0x0,
			mismatch,
			insertion,	// in the source sequence
			deletion
		};
		
		typedef enum_wrapper <operation_type, operation_type::deletion>	operation;
		
		
		std::vector <pair <gapped_nucleotide>>	sequence;
		std::size_t								src_position{};
		std::size_t								dst_position{};
		
		msa_non_segment_section() = default;
		
		msa_non_segment_section(
			std::vector <nucleotide> const &sequence_,
			std::vector <operation> const &operations_
		):
			sequence(sequence_.size())
		{
			libbio_assert_eq(sequence_.size(), operations_.size());
			for (auto &&[cc, op, dst] : rsv::zip(sequence_, operations_, sequence))
			{
				switch (op.value)
				{
					case operation_type::match:
						dst = {cc, cc};
						++src_position;
						++dst_position;
						break;
					
					case operation_type::mismatch:
					{
						dst = {cc, cc.next_nucleotide()};
						++src_position;
						++dst_position;
						break;
					}
					
					case operation_type::insertion:
					{
						dst = {cc, '-'_gapped_nucleotide};
						++src_position;
						break;
					}
					
					case operation_type::deletion:
					{
						dst = {'-'_gapped_nucleotide, cc};
						++dst_position;
						break;
					}
				}
			}
		}
		
		bool empty() const { return sequence.empty(); }
		
		template <std::size_t t_idx, alignment_output t_alignment_output>
		auto sequence_at_index() const
		{
			auto rng(sequence | rsv::transform([](auto const pair){ return std::get <t_idx>(pair); }));
			if constexpr (t_alignment_output == alignment_output::unaligned)
				return rng | rsv::filter([](auto const cc){ return '-'_gapped_nucleotide != cc; });
			else
				return rng;
		}
	};
	
	
	// Part of the MSA that contains the “boundary” of the aligned segment, i.e. clipping and padding.
	// (SAM requires that the operations have a certain order but we do not.)
	struct msa_boundary_section
	{
		enum class operation_type : std::uint8_t
		{
			soft_clipping = 0x0,
			hard_clipping,
			padding
		};
		
		typedef enum_wrapper <operation_type, operation_type::padding>	operation;
		
		
		std::vector <nucleotide>	segment;
		panvc3::cigar_vector		cigar;
		
		static std::size_t expected_sequence_length(std::vector <operation> const &ops)
		{
			return std::count(ops.begin(), ops.end(), operation(operation_type::soft_clipping));
		}
		
		
		msa_boundary_section() = default;
		
		msa_boundary_section(
			std::vector <nucleotide> const &segment_,
			std::vector <operation> const &operations_
		):
			segment(segment_)
		{
			// Since soft clipping is the only consuming operation here,
			// the caller should use rc::gen::mapcat to get a sequence of the appropriate length.
			panvc3::cigar_buffer cigar_buffer;
			for (auto const op : operations_)
			{
				switch (op.value)
				{
					case operation_type::soft_clipping:
						cigar_buffer.push_back('S'_cigar_operation);
						break;
						
					case operation_type::hard_clipping:
						cigar_buffer.push_back('H'_cigar_operation);
						break;
						
					case operation_type::padding:
						cigar_buffer.push_back('P'_cigar_operation);
						break;
				}
			}
			
			cigar_buffer.finish();
			cigar = std::move(cigar_buffer.operations());
		}
		
		bool empty() const { return segment.empty(); }
	};
	
	
	// Part of the MSA that contains the aligned segment.
	struct msa_segment_section
	{
		enum class operation_type : std::uint8_t
		{
			// Matches and mismatches
			segment_matches_both = 0x0,			// AAA	(segment character, src character, dst character)
			segment_matches_src_mismatches_dst,	// AAB
			segment_mismatches_src_matches_dst,	// ABA
			segment_mismatches_both,			// ABC
			// Insertions
			segment_matches_src_insertion,		// AA-
			segment_mismatches_src_insertion,	// AB-
			segment_insertion,					// A--
			// Deletions
			dst_deletion,						// --A
			segment_deletion_match,				// -AA
			segment_deletion_mismatch,			// -AB
			// No-ops
			src_deletion,						// -A-
			padding
		};
		
		typedef enum_wrapper <operation_type, operation_type::padding>	operation;
		
		constexpr static auto const match_and_mismatch_operations{std::to_array <operation>({
			operation_type::segment_matches_both,				// AAA	(segment character, src character, dst character)
			operation_type::segment_matches_src_mismatches_dst,	// AAB
			operation_type::segment_mismatches_src_matches_dst,	// ABA
			operation_type::segment_mismatches_both				// ABC
		})};
		
		
		// MSA.
		std::vector <triplet <gapped_nucleotide>>	sequence;
		panvc3::cigar_vector						input_cigar;
		panvc3::cigar_vector						expected_cigar;
		std::vector <operation>						operations; // For debugging.
		
		
		// Only padding does not consume anything.
		static std::size_t expected_sequence_length(std::vector <operation> const &ops)
		{
			return std::count_if(ops.begin(), ops.end(), [](auto const op){ return op.value != operation_type::padding; });
		}
		
		static std::tuple <triplet <gapped_nucleotide>, seqan3::cigar::operation, seqan3::cigar::operation>
		msa_operation(nucleotide const cc, operation const op)
		{
			switch (op.value)
			{
				// Matches
				case operation_type::segment_matches_both:																	// AAA	(segment character, src character, dst character)
					return {{cc, cc, cc}, '='_cigar_operation, '='_cigar_operation};
				
				case operation_type::segment_matches_src_mismatches_dst:													// AAB
					return {{cc, cc, cc.next_nucleotide()}, '='_cigar_operation, 'X'_cigar_operation};
				
				case operation_type::segment_mismatches_src_matches_dst:													// ABA
					return {{cc, cc.next_nucleotide(), cc}, 'X'_cigar_operation, '='_cigar_operation};
				
				case operation_type::segment_mismatches_both:																// ABC
				{
					auto const ncc(cc.next_nucleotide());
					auto const ncc_(ncc.next_nucleotide());
					return {{cc, ncc, ncc_}, 'X'_cigar_operation, 'X'_cigar_operation};
				}
				
				// Insertions
				case operation_type::segment_matches_src_insertion:															// AA-
					return {{cc, cc, '-'_gapped_nucleotide}, '='_cigar_operation, 'I'_cigar_operation};
				
				case operation_type::segment_mismatches_src_insertion:														// AB-
					return {{cc, cc.next_nucleotide(), '-'_gapped_nucleotide}, 'X'_cigar_operation, 'I'_cigar_operation};
				
				case operation_type::segment_insertion:																		// A--
					return {{cc, '-'_gapped_nucleotide, '-'_gapped_nucleotide}, 'I'_cigar_operation, 'I'_cigar_operation};
				
				// Deletions
				case operation_type::segment_deletion_match:																// -AA
					return {{'-'_gapped_nucleotide, cc, cc}, 'D'_cigar_operation, 'D'_cigar_operation};
				
				case operation_type::segment_deletion_mismatch:																// -AB
					return {{'-'_gapped_nucleotide, cc, cc.next_nucleotide()}, 'D'_cigar_operation, 'D'_cigar_operation};
				
				// No-ops, handled by the caller separately.
				case operation_type::src_deletion:																			// -A-
				case operation_type::dst_deletion:																			// --A
				case operation_type::padding:
					libbio_fail("Unexpected operation");
			}
			
			// Not reached but keep the compiler happy by returning something.
			libbio_fail("Should not be reached");
			return {{'-'_gapped_nucleotide, '-'_gapped_nucleotide, '-'_gapped_nucleotide}, 'P'_cigar_operation, 'P'_cigar_operation};
		}
		
		msa_segment_section() = default;
		
		msa_segment_section(
			std::vector <nucleotide> const &segment_,
			std::vector <operation> const &operations_
		):
			sequence(segment_.size()),
			operations(operations_)
		{
			// The caller should use rc::gen::mapcat to get a sequence of the appropriate length.
			panvc3::cigar_buffer input_cigar_buffer;
			panvc3::cigar_buffer expected_cigar_buffer;
			std::size_t pos{};
			for (auto const op : operations_)
			{
				switch (op.value)
				{
					case operation_type::segment_matches_both:					// AAA	(segment character, src character, dst character)
					case operation_type::segment_matches_src_mismatches_dst:	// AAB
					case operation_type::segment_mismatches_src_matches_dst:	// ABA
					case operation_type::segment_mismatches_both:				// ABC
					case operation_type::segment_matches_src_insertion:			// AA-
					case operation_type::segment_mismatches_src_insertion:		// AB-
					case operation_type::segment_insertion:						// A--
					case operation_type::segment_deletion_match:				// -AA
					case operation_type::segment_deletion_mismatch:				// -AB
					{
						auto const cc(segment_[pos]);
						auto const res(msa_operation(cc, op));
						auto const &[characters, input_op, expected_op] = res;
						sequence[pos] = characters;
						input_cigar_buffer.push_back(input_op);
						expected_cigar_buffer.push_back(expected_op);
						++pos;
						break;
					}
				
					// Others
					case operation_type::dst_deletion:							// --A
					{
						auto const cc(segment_[pos]);
						sequence[pos] = {'-'_gapped_nucleotide, '-'_gapped_nucleotide, cc};
						expected_cigar_buffer.push_back('D'_cigar_operation);
						++pos;
						break;
					}
				
					case operation_type::src_deletion:							// -A-
					{
						auto const cc(segment_[pos]);
						sequence[pos] = {'-'_gapped_nucleotide, cc, '-'_gapped_nucleotide};
						input_cigar_buffer.push_back('D'_cigar_operation);
						++pos;
						break;
					}
				
					case operation_type::padding:
						input_cigar_buffer.push_back('P'_cigar_operation);
						expected_cigar_buffer.push_back('P'_cigar_operation);
						break;
				}
			}
			
			libbio_assert_eq(pos, sequence.size());
			
			input_cigar_buffer.finish();
			expected_cigar_buffer.finish();
			input_cigar = std::move(input_cigar_buffer.operations());
			expected_cigar = std::move(expected_cigar_buffer.operations());
		}
		
		auto segment_unaligned() const
		{
			return sequence
			| rsv::transform([](auto const tup){ return std::get <0>(tup); })
			| rsv::filter([](auto const cc){ return '-'_gapped_nucleotide != cc; });
		}
		
		auto alignment_without_segment() const
		{
			return sequence
			| rsv::transform([](auto const tup) -> pair <gapped_nucleotide> { return {std::get <1>(tup), std::get <2>(tup)}; })
			| rsv::filter([](auto const pp){ return std::make_tuple('-'_gapped_nucleotide, '-'_gapped_nucleotide) != pp; });
		}
		
		// If t_alignment_output == alignment_output::align_all, t_idx == 0 causes the aligned segment to be returned.
		template <std::size_t t_idx, alignment_output t_alignment_output>
		auto sequence_at_index() const
		{
			if constexpr (t_alignment_output == alignment_output::align_all)
				return sequence | rsv::transform([](auto const tup){ return std::get <t_idx>(tup); });
			else
			{
				auto rng(alignment_without_segment() | rsv::transform([](auto const pair){ return std::get <t_idx>(pair); }));
				if constexpr (t_alignment_output == alignment_output::align_sequences)
					return rng;
				else
					return rng | rsv::filter([](auto const cc){ return '-'_gapped_nucleotide != cc; });
			}
		}
		
		bool empty() const { return sequence.empty(); }
	};
	
	
	// Build a complete MSA of two sequences and one aligned segment.
	// NOTE: Needs to be SimpleAggregate for Boost.PFR to work.
	struct msa_builder
	{
		// We store these mostly for debugging purposes.
		msa_segment_section		segment_section;
		msa_boundary_section	left_boundary;
		msa_boundary_section	right_boundary;
		msa_non_segment_section	left_non_segment;
		msa_non_segment_section	right_non_segment;
		
		auto segment_unaligned() const
		{
			return rsv::concat(
				left_boundary.segment,
				segment_section.segment_unaligned(),
				right_boundary.segment
			);
		}
		
		template <std::size_t t_idx, alignment_output t_alignment_output>
		auto sequence_at_index() const
		{
			return rsv::concat(
				left_non_segment.sequence_at_index <t_idx, t_alignment_output>(),
				segment_section.sequence_at_index <t_idx, t_alignment_output>(),
				right_non_segment.sequence_at_index <t_idx, t_alignment_output>()
			);
		}
		
		auto input_cigar() const
		{
			return rsv::concat(
				left_boundary.cigar,
				segment_section.input_cigar,
				right_boundary.cigar
			);
		}
		
		auto expected_cigar() const
		{
			return rsv::concat(
				left_boundary.cigar,
				segment_section.expected_cigar,
				right_boundary.cigar
			);
		}
	};
	
	
	// For RapidCheck.
	void showValue(seqan3::cigar const cigar_item, std::ostream &os)
	{
		using seqan3::get;
		
		auto const op_count(get <0>(cigar_item));
		auto const operation(get <1>(cigar_item));
		auto const operation_(operation.to_char());
		
		os << op_count << operation_;
	}
	
	
	void showValue(std::vector <seqan3::cigar> const &cigar_vec, std::ostream &os)
	{
		for (auto const cc : cigar_vec)
			showValue(cc, os);
	}
	
	
	void showValue(msa_builder const &builder, std::ostream &os)
	{
		auto left_non_segment_src_seq(builder.left_non_segment.sequence_at_index <0, alignment_output::align_all>());
		auto left_non_segment_dst_seq(builder.left_non_segment.sequence_at_index <1, alignment_output::align_all>());
		auto right_non_segment_src_seq(builder.right_non_segment.sequence_at_index <0, alignment_output::align_all>());
		auto right_non_segment_dst_seq(builder.right_non_segment.sequence_at_index <1, alignment_output::align_all>());
		auto left_boundary_seq(builder.left_boundary.segment);
		auto right_boundary_seq(builder.right_boundary.segment);
		auto segment_read_seq(builder.segment_section.sequence_at_index <0, alignment_output::align_all>());
		auto segment_src_seq(builder.segment_section.sequence_at_index <1, alignment_output::align_all>());
		auto segment_dst_seq(builder.segment_section.sequence_at_index <2, alignment_output::align_all>());
		
		auto it{std::ostream_iterator <char>(os)};
		auto const left_non_segment_size(ranges::size(left_non_segment_src_seq));
		auto const right_non_segment_size(ranges::size(right_non_segment_src_seq));
		auto const left_boundary_size(ranges::size(left_boundary_seq));
		auto const right_boundary_size(ranges::size(right_boundary_seq));
		auto const seg_read_size(ranges::size(segment_read_seq));
		auto const seg_src_size(ranges::size(segment_src_seq));
		auto const seg_dst_size(ranges::size(segment_dst_seq));
		
		os << "Alignment:  ";
		std::fill_n(it, left_non_segment_size, '_');
		ranges::copy(left_boundary_seq, it);
		ranges::copy(segment_read_seq, it);
		ranges::copy(right_boundary_seq, it);
		os << " (lb: " << left_boundary_size << ", sr: " << seg_read_size << ", rb: " << right_boundary_size << ")\n";
		
		os << "Source:      ";
		ranges::copy(left_non_segment_src_seq, it);
		std::fill_n(it, left_boundary_size, '-');
		ranges::copy(segment_src_seq, it);
		std::fill_n(it, right_boundary_size, '-');
		ranges::copy(right_non_segment_src_seq, it);
		os << " (lns: " << left_non_segment_size << ", lb: " << left_boundary_size << ", src: " << seg_src_size << ", rb: " << right_boundary_size << ", rns: " << right_non_segment_size << ")\n";
		
		os << "Destination: ";
		ranges::copy(left_non_segment_dst_seq, it);
		std::fill_n(it, left_boundary_size, '-');
		ranges::copy(segment_dst_seq, it);
		std::fill_n(it, right_boundary_size, '-');
		ranges::copy(right_non_segment_dst_seq, it);
		os << " (lns: " << left_non_segment_size << ", lb: " << left_boundary_size << ", dst: " << seg_dst_size << ", rb: " << right_boundary_size << ", rns: " << right_non_segment_size << ")\n";
		
		os << "Input CIGAR:    ";
		for (auto const cc : builder.input_cigar()) showValue(cc, os);
		os << '\n';
		os << "Expected CIGAR: ";
		for (auto const cc : builder.expected_cigar()) showValue(cc, os);
		os << '\n';
		os << "Expected dst position: " << left_non_segment_size << '\n';
		
		os << "Segment section operations: ";
		ranges::copy(
			builder.segment_section.operations
			| rsv::transform([](auto const cc){ return std::uint32_t(cc); }),
			ranges::make_ostream_joiner(os, ", ")
		);
	}
}


namespace rc {
	
	template <>
	struct Arbitrary <nucleotide>
	{
		static Gen <nucleotide> arbitrary()
		{
			return gen::construct <nucleotide>(gen::elementOf(nucleotide::values));
		}
	};
	
	
	template <>
	struct Arbitrary <gapped_nucleotide>
	{
		static Gen <gapped_nucleotide> arbitrary()
		{
			return gen::construct <gapped_nucleotide>(gen::elementOf(gapped_nucleotide::values));
		}
	};
	
	
	template <typename t_enum, t_enum t_max, t_enum t_min>
	struct Arbitrary <enum_wrapper <t_enum, t_max, t_min>>
	{
		typedef enum_wrapper <t_enum, t_max, t_min>	value_type;
		
		static Gen <value_type> arbitrary()
		{
			return gen::construct <value_type>(gen::inClosedRange(value_type::min, value_type::max));
		}
	};
	
	
	template <>
	struct Arbitrary <msa_non_segment_section>
	{
		static Gen <msa_non_segment_section> arbitrary()
		{
			typedef msa_non_segment_section::operation operation;
			
			return gen::withSize([](int const size){
				// The vectors are expected to have equal lengths.
				return gen::construct <msa_non_segment_section>(
					gen::container <std::vector <nucleotide>>(size, gen::arbitrary <nucleotide>()),
					gen::container <std::vector <operation>>(size, gen::arbitrary <operation>())
				);
			});
		}
	};
}


namespace {
	
	// Post-processing for msa_segment_section.
	template <typename t_section>
	struct fix_operations_for_arbitrary_msa_section
	{
		template <typename t_cb>
		static auto fix(std::vector <typename t_section::operation> const &operations, t_cb &&cb) -> rc::Gen <t_section>
		{
			return cb(operations);
		}
	};
	
	
	template <>
	struct fix_operations_for_arbitrary_msa_section <msa_segment_section>
	{
		template <typename t_cb>
		static auto fix(std::vector <msa_segment_section::operation> const &operations, t_cb &&cb) -> rc::Gen <msa_segment_section>
		{
			// Deletions may not occur on the boundary of msa_segment_section,
			// so we add matches if needed.
			auto const has_deletions_on_boundary([](auto const &rng){
				for (auto const op : rng)
				{
					switch (op.value)
					{
						case msa_segment_section::operation_type::segment_matches_both:					// AAA	(segment character, src character, dst character)
						case msa_segment_section::operation_type::segment_matches_src_mismatches_dst:	// AAB
						case msa_segment_section::operation_type::segment_mismatches_src_matches_dst:	// ABA
						case msa_segment_section::operation_type::segment_mismatches_both:				// ABC
							return false;
						
						case msa_segment_section::operation_type::dst_deletion:							// --A
						case msa_segment_section::operation_type::segment_deletion_match:				// -AA
						case msa_segment_section::operation_type::segment_deletion_mismatch:			// -AB
							return true;
							
						default:
							break;
					}
				}
				
				return false;
			});
			
			auto const res(has_deletions_on_boundary(operations) | (has_deletions_on_boundary(rsv::reverse(operations)) << 0x1));
			if (!res)
				return cb(operations);
			
			auto const new_operation([](){
				return rc::gen::elementOf(msa_segment_section::match_and_mismatch_operations);
			});
			
			switch (res)
			{
				case 0x1: // Left end has deletions.
					return rc::gen::mapcat(
						rc::gen::tuple(rc::gen::just(operations), new_operation()),
						[cb = std::forward <t_cb>(cb)](auto const &tup){
							auto operations(std::get <0>(tup)); // Copy.
							operations.insert(operations.begin(), std::get <1>(tup));
							return cb(std::move(operations));
						}
					);
				
				case 0x2: // Right end has deletions.
					return rc::gen::mapcat(
						rc::gen::tuple(rc::gen::just(operations), new_operation()),
						[cb = std::forward <t_cb>(cb)](auto const &tup){
							auto operations(std::get <0>(tup)); // Copy.
							operations.push_back(std::get <1>(tup));
							return cb(std::move(operations));
						}
					);
				
				case 0x3: // Both ends have deletions.
					return rc::gen::mapcat(
						rc::gen::tuple(rc::gen::just(operations), new_operation(), new_operation()),
						[cb = std::forward <t_cb>(cb)](auto const &tup) {
							auto operations(std::get <0>(tup)); // Copy.
							operations.insert(operations.begin(), std::get <1>(tup));
							operations.push_back(std::get <2>(tup));
							return cb(std::move(operations));
						}
					);
				
				default:
					libbio_fail("Unexpected result");
					return cb(operations); // Not reached.
			}
		}
	};
	
	
	// Construction helper, passes the operations and the sequence to t_msa’s constructor.
	template <typename t_msa>
	rc::Gen <t_msa> arbitrary_msa_section()
	{
		typedef typename t_msa::operation operation;
		
		return rc::gen::withSize([](int const size){
			return rc::gen::mapcat(
				rc::gen::container <std::vector <operation>>(size, rc::gen::arbitrary <operation>()),
				[size](auto const &operations){
					
					typedef fix_operations_for_arbitrary_msa_section <t_msa> fix_type;
					return fix_type::fix(operations, []<typename t_operations>(t_operations &&operations){
						// Build a nucleotide sequence whose length depends on the operations.
						auto const seq_length(t_msa::expected_sequence_length(operations));
						return rc::gen::construct <t_msa>(
							rc::gen::container <std::vector <nucleotide>>(
								seq_length,
								rc::gen::arbitrary <nucleotide>()
							),
							rc::gen::just(std::forward <t_operations>(operations))
						);
					});
				}
			);
		});
	}
}


namespace rc {
	
	template <>
	struct Arbitrary <msa_boundary_section>
	{
		static Gen <msa_boundary_section> arbitrary() { return arbitrary_msa_section <msa_boundary_section>(); }
	};
	
	
	template <>
	struct Arbitrary <msa_segment_section>
	{
		static Gen <msa_segment_section> arbitrary()
		{
			// Try to shrink by consuming = operations from both ends.
			return gen::shrink(
				arbitrary_msa_section <msa_segment_section>(),
				[](msa_segment_section &&section){
					using seqan3::get;
					
					auto &input_cigar(section.input_cigar);
					auto &expected_cigar(section.expected_cigar);
					if (input_cigar.empty() || expected_cigar.empty())
						return Seq <msa_segment_section>();
					
					auto const contract([](auto &input_item, auto &expected_item, auto &&erase_fn){
						auto const is_eq_op([](auto const &cigar_item){
							return '='_cigar_operation == get <1>(cigar_item);
						});
						
						if (is_eq_op(input_item) && is_eq_op(expected_item))
						{
							auto input_amt(get <0>(input_item));
							auto expected_amt(get <0>(expected_item));
							auto const amt(std::min(input_amt, expected_amt));
							
							// Keep one operation.
							if (1 < amt)
							{
								// The values returned by get <0> are proxies, hence we cannot use operator -=.
								auto const amt_(amt - 1);
								input_amt = input_amt - amt_;
								expected_amt = expected_amt - amt_;
								erase_fn(amt_);
								
								return true;
							}
						}
						
						return false;
					});
					
					auto &seq(section.sequence);
					auto &ops(section.operations);
					auto const r1(contract(input_cigar.front(), expected_cigar.front(), [&seq, &ops](auto const size){
						seq.erase(seq.begin(), seq.begin() + size);
						ops.erase(ops.begin(), ops.begin() + size);
					}));
					auto const r2(contract(input_cigar.back(), expected_cigar.back(), [&seq, &ops](auto const size){
						seq.erase(seq.end() - size, seq.end());
						ops.erase(ops.end() - size, ops.end());
					}));
					
					// Check if anything was contracted.
					if (r1 || r2)
						return seq::just(std::move(section));
					
					return Seq <msa_segment_section>();
				}
			);
		}
	};
	
	
	template <>
	struct Arbitrary <msa_builder>
	{
		static Gen <msa_builder> arbitrary()
		{
			// Scale the sizes of the sections. We achieve this by getting four factors, scaling them with
			// predefined multipliers (below) and calculating the size of the msa_segment_section by subtracting
			// the sizes of the other sections from the total size.
			
			// RapidCheck has generators for integers, so we need to convert to floating point manually.
			typedef std::uint16_t factor_type;	// Enough resolution hopefully.
			constexpr factor_type const factor_max{std::numeric_limits <factor_type>::max()};
			
			// We still need the current size.
			return gen::withSize([](int size) {
				return gen::mapcat(
					gen::tuple(
						gen::inClosedRange <factor_type>(0, factor_max),
						gen::inClosedRange <factor_type>(0, factor_max),
						gen::inClosedRange <factor_type>(0, factor_max),
						gen::inClosedRange <factor_type>(0, factor_max)
					),
					[]<typename... t_args>(std::tuple <t_args...> const factors) {
						// This can be (perhaps) removed when we have P2141R0, perhaps in C++26.
						typedef structure_to_tuple_t <msa_builder> section_tuple_type;
						constexpr auto const section_count(std::tuple_size_v <section_tuple_type>);
						constexpr auto const factor_count(sizeof...(t_args));
						static_assert(section_count - 1 == factor_count);
						
						// The multipliers for left and right boundaries and left and right non-segments, respectively.
						constexpr std::array <double, factor_count> const multipliers{1.0, 1.0, 3.0, 2.0}; // Parts of 12.
						
						// Call gen::construct <msa_builder>() with the different section types as the parameters.
						// To this end, we scale the size value individually for each section type before generating an arbitrary value.
						// Shrinking might not be needed from here.
						return std::apply(
							// For some reason passing gen::construct <msa_builder> directly as the callback causes an error.
							[](auto && ... values){ return gen::construct <msa_builder>(std::forward <decltype(values)>(values)...); },
							panvc3::map_to_tuple(
								std::make_index_sequence <section_count>{},
								// Not specifying the return type causes (at least) GCC to make a wrong deduction.
								[&factors, &multipliers]<typename t_idx>(t_idx const idx_) -> Gen <std::tuple_element_t <t_idx::value, section_tuple_type>> {
									constexpr auto const idx(t_idx::value);
									typedef std::tuple_element_t <idx, section_tuple_type> section_type;
									
									// The scaling function.
									auto const scale_fn([](auto const factor, auto const multiplier) -> double {
										return multiplier / 12.0 * (1.0 * factor / factor_max);
									});
									
									if constexpr (0 == idx)
									{
										// The segment_section gets the complement of the sum of the other scale factors as its scale factor.
										auto const weighted_factors(panvc3::map_to_array(
											std::make_index_sequence <factor_count>{},
											[&](auto const idx_){
												constexpr auto const idx(decltype(idx_)::value);
												return scale_fn(std::get <idx>(factors), std::get <idx>(multipliers));
											}
										));
											
										auto const scale_factor_(std::accumulate(weighted_factors.begin(), weighted_factors.end(), 0.0));
										libbio_assert_lte(0.0, scale_factor_);
										auto const scale_factor(1.0 - scale_factor_);
										return gen::scale(scale_factor, gen::arbitrary <section_type>());
									}
									else
									{
										// Otherwise just get the factors, scale the size, and produce the correct section type.
										constexpr auto const multiplier(std::get <idx - 1>(multipliers));
										auto const factor(std::get <idx - 1>(factors));
										return gen::scale(scale_fn(factor, multiplier), gen::arbitrary <section_type>());
									}
								}
							)
						);
					}
				);
			});
		}
	};
}


TEST_CASE("rewrite_cigar() can process an arbitrary MSA", "[rewrite_cigar]")
{
	rc::prop(
		"rewrite_cigar() works as expected",
		[](msa_builder const &builder){
			
			bool const none_empty(
				builder.segment_section.empty() &&
				builder.left_boundary.empty() &&
				builder.right_boundary.empty() &&
				builder.left_non_segment.empty() &&
				builder.right_non_segment.empty()
			);
			RC_CLASSIFY(none_empty);
			
			auto const src_pos(builder.left_non_segment.src_position);
			auto const expected_dst_pos(builder.left_non_segment.dst_position);
			auto const input_cigar(to_vector(builder.input_cigar()));
			auto const expected_cigar(to_vector(builder.expected_cigar()));
			auto const segment_unaligned(to_vector(builder.segment_unaligned()));
			auto const dst_unaligned(to_vector(builder.sequence_at_index <1, alignment_output::unaligned>()));
			
			// For some reason the ranges returned by builder are not sized.
			auto const src_aligned(to_vector(builder.sequence_at_index <0, alignment_output::align_sequences>()));
			auto const dst_aligned(to_vector(builder.sequence_at_index <1, alignment_output::align_sequences>()));
			
			panvc3::sequence_entry_pair seq_entries;
			panvc3::make_sequence_entry_pair(src_aligned, dst_aligned, seq_entries);
			
			panvc3::cigar_buffer dst_cigar_buffer;
			auto const dst_pos(panvc3::rewrite_cigar(
				src_pos,
				input_cigar,
				seq_entries.first,
				seq_entries.second,
				segment_unaligned,
				dst_unaligned,
				dst_cigar_buffer
			));
			
			auto const &actual_cigar(dst_cigar_buffer.operations());
			if (dst_pos != expected_dst_pos || !panvc3::cigar_eq(expected_cigar, actual_cigar))
			{
				std::stringstream os;
				os << "Failed.\n";
				os << "Actual dst position: " << dst_pos << '\n';
				os << "Actual CIGAR: ";
				showValue(actual_cigar, os);
				os << '\n';
				
				RC_FAIL(os.str());
			}
		}
	);
}
