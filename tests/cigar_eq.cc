/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <libbio/assert.hh>
#include <panvc3/cigar.hh>
#include <panvc3/cigar_eq.hh>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/drop_exactly.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/take_exactly.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>			// rc::prop
#include "rapidcheck_additions.hh"		// rc::gen::inClosedRange


namespace gen	= Catch::Generators;
namespace lb	= libbio;
namespace rsv	= ranges::views;

using seqan3::operator""_cigar_operation;


namespace {
	
	template <typename t_type>
	using pair = std::pair <t_type, t_type>;
	
	
	// Helper function for getting an element in the middle of a vector.
	// (Sanity checking not done.)
	template <typename t_vec>
	auto &middle_(t_vec &&vec)
	{
		auto const size(vec.size());
		libbio_assert_lt(0, size);
		return vec[size / 2];
	}
	
	template <typename t_type> t_type &middle(std::vector <t_type> &vec) { return middle_(vec); }
	template <typename t_type> t_type const &middle(std::vector <t_type> const &vec) { return middle_(vec); }
	
	
	// Helper function for building a range that conditionally holds the given value.
	// Optionally dereferences the given value on access. (Useful for e.g. std::optional.)
	template <bool t_should_deref = false, typename t_value>
	auto conditional_range(t_value &&value, bool const should_take)
	{
		if constexpr (t_should_deref)
			return rsv::single(value) | rsv::take_exactly(should_take) | rsv::transform([](auto &value){ return *value; });
		else
			return rsv::single(std::forward <t_value>(value)) | rsv::take_exactly(should_take);
	}
	
	template <bool t_should_deref = true, typename t_value>
	auto conditional_range(t_value &&value)
	{
		bool const should_take(value);
		return conditional_range <t_should_deref>(std::forward <t_value>(value), should_take);
	}
	
	
	typedef panvc3::cigar_adapter_seqan3::count_type	cigar_count_type;
	typedef std::vector <cigar_count_type>				cigar_count_vector;
	
	
	constexpr static auto const cigar_non_indel_operations(std::to_array({
		'M'_cigar_operation,
		'N'_cigar_operation,
		'S'_cigar_operation,
		'H'_cigar_operation,
		'P'_cigar_operation,
		'='_cigar_operation,
		'X'_cigar_operation
	}));
	
	
	// Helper for constructing seqan3::cigar.
	seqan3::cigar make_cigar_item(cigar_count_type const count, seqan3::cigar::operation const op)
	{
		seqan3::cigar retval;
		retval = count;
		retval = op;
		return retval;
	}
	
	
	// Represents a non-indel in CIGAR, i.e. M, N, S, H, P, =, or X.
	struct non_indel_item
	{
		seqan3::cigar	item;
		
		void set_operation(seqan3::cigar::operation const op) { item = op; }
		void set_count(cigar_count_type const count) { item = count; }
		void append_sequence_to(std::vector <seqan3::cigar> &dst) const { dst.emplace_back(item); }
	};
	
	
	// Represents an insertion followed by a deletion.
	struct indel_item
	{
		cigar_count_type	ins;
		cigar_count_type	del;
		
		void append_sequence_to(std::vector <seqan3::cigar> &dst) const
		{
			libbio_assert_lt(0, ins);
			libbio_assert_lt(0, del);
			
			dst.emplace_back(make_cigar_item(ins, 'I'_cigar_operation));
			dst.emplace_back(make_cigar_item(del, 'D'_cigar_operation));
		}
	};
	
	
	// Run of insertions and deletions.
	struct indel_run
	{
		std::vector <indel_item>	indels;
		std::uint16_t				leading_deletions{};
		std::uint16_t				trailing_insertions{};
		
		void append_sequence_to(std::vector <seqan3::cigar> &dst) const
		{
			if (leading_deletions) dst.emplace_back(make_cigar_item(leading_deletions, 'D'_cigar_operation));
			
			for (auto const &indel : indels)
				indel.append_sequence_to(dst);
			
			if (trailing_insertions) dst.emplace_back(make_cigar_item(trailing_insertions, 'I'_cigar_operation));
		}
	};
	
	
	// The middle part of a pair of CIGAR sequences; has two (possibly equivalent) runs of indels
	// and a common non-indel.
	struct cigar_middle_part
	{
		indel_run		lhs_indels;
		indel_run		rhs_indels;
		non_indel_item	non_indel;
		
		bool has_all_parts() const
		{
			return !(lhs_indels.indels.empty() || rhs_indels.indels.empty());
		}
	};
	
	
	// The tail part of a pair of CIGAR sequences; has two (possibly equivalent) runs of indels.
	struct cigar_tail_part
	{
		indel_run		lhs_indels;
		indel_run		rhs_indels;
		
		bool has_all_parts() const
		{
			return !(lhs_indels.indels.empty() || rhs_indels.indels.empty());
		}
	};
	
	
	// Complete pair of two CIGAR sequences.
	struct cigar_sequence_pair
	{
		std::vector <cigar_middle_part>	mid;
		rc::Maybe <non_indel_item>		head;
		rc::Maybe <cigar_tail_part>		tail;
		
		template <typename t_access_fn>
		auto sequence(t_access_fn &&fn) const
		{
			std::vector <seqan3::cigar> retval;
			
			if (head)
				head->append_sequence_to(retval);
			
			for (auto const &mid_part : mid)
			{
				fn(mid_part).append_sequence_to(retval);
				mid_part.non_indel.append_sequence_to(retval);
			}
			
			if (tail)
				fn(*tail).append_sequence_to(retval);
			
			return retval;
		}
		
		std::vector <seqan3::cigar> lhs() const { return sequence([](auto const &part){ return part.lhs_indels; }); }
		std::vector <seqan3::cigar> rhs() const { return sequence([](auto const &part){ return part.rhs_indels; }); }
	};
	
	
	// Set of related test inputs for cigar_eq(); has a matching sequence pair
	// and a number of non-matching sequence pairs.
	struct cigar_test_input
	{
		cigar_sequence_pair					matching_pair;
		std::vector <cigar_sequence_pair>	non_matching_pairs;
		
		// For debugging.
		std::vector <bool>					mutation_targets;
		std::vector <bool>					mutation_is_del;
	};
	
	
	void showValue(seqan3::cigar const cigar_item, std::ostream &os)
	{
		using seqan3::get;
		
		auto const op_count(get <0>(cigar_item));
		auto const operation(get <1>(cigar_item));
		auto const operation_(operation.to_char());
		
		os << op_count << operation_;
	}
	
	
	void showValue(cigar_sequence_pair const &input, std::ostream &os)
	{
		os << " lhs:             ";
		for (auto const cc : input.lhs()) showValue(cc, os);
		os << '\n';
		
		os << " rhs:             ";
		for (auto const cc : input.rhs()) showValue(cc, os);
		os << '\n';
	}
	
	
	void showValue(cigar_test_input const &input, std::ostream &os)
	{
		os << "matching_pair:\n";
		showValue(input.matching_pair, os);
		
		os << "non_matching_pairs:\n";
		for (auto const &[idx, pair]: rsv::enumerate(input.non_matching_pairs))
		{
			os << idx << '\n';
			showValue(pair, os);
		}
		
		auto it{std::ostream_iterator <int>(os)};
		
		os << "mutation_targets: ";
		ranges::copy(input.mutation_targets, it);
		os << '\n';
		
		os << "mutation_is_del:  ";
		ranges::copy(input.mutation_is_del, it);
		os << '\n';
	}
	
	
	// For placing a mutation in some part of a pair of CIGAR sequences.
	enum class mutation_position
	{
		left,
		middle,
		right
	};
}


namespace rc {
	
	template <>
	struct Arbitrary <non_indel_item>
	{
		static Gen <non_indel_item> arbitrary()
		{
			return gen::build(
				gen::construct <non_indel_item>(),
				gen::set(&non_indel_item::set_operation, gen::elementOf(cigar_non_indel_operations)),
				gen::set(&non_indel_item::set_count, gen::inClosedRange(1, 10))
			);
		}
	};
	
	
	template <>
	struct Arbitrary <cigar_sequence_pair>
	{
		static Gen <cigar_sequence_pair> arbitrary()
		{
			return gen::withSize([](int size){
				libbio_assert_lte(0, size);
				return gen::exec([size](bool const has_leading_non_indel, bool const has_trailing_indel_run) -> cigar_sequence_pair {
					
					cigar_sequence_pair retval;
					
					auto fill_indel_run([](auto ins_count, auto del_count, indel_run &ir){
						// Add the insertions and the deletions.
						while (ins_count && del_count)
						{
							auto const ins(*gen::inClosedRange(cigar_count_type(1), ins_count));
							auto const del(*gen::inClosedRange(cigar_count_type(1), del_count));
							ins_count -= ins;
							del_count -= del;
							libbio_assert_lt(0, ins);
							libbio_assert_lt(0, del);
							ir.indels.emplace_back(ins, del);
						}
						
						// Handle the remaining insertions and deletions.
						ir.leading_deletions = del_count;
						ir.trailing_insertions = ins_count;
					});
					
					// Insertion and deletion counts.
					auto const insertion_counts(*gen::container <cigar_count_vector>(size, gen::inClosedRange(1, 30)));
					auto const deletion_counts(*gen::container <cigar_count_vector>(size, gen::inClosedRange(1, 30)));
					
					// Add the leading non-indel if needed.
					if (has_leading_non_indel)
						retval.head = *gen::arbitrary <non_indel_item>();
					
					if (size)
					{
						// Add the trailing indel run if needed.
						if (has_trailing_indel_run)
						{
							retval.tail = cigar_tail_part{};
							fill_indel_run(insertion_counts.front(), deletion_counts.front(), retval.tail->lhs_indels);
							fill_indel_run(insertion_counts.front(), deletion_counts.front(), retval.tail->rhs_indels);
						}
						
						// Add the middle parts.
						retval.mid.resize(size - has_trailing_indel_run);
						if (!retval.mid.empty())
						{
							auto rng(
								rsv::zip(
									insertion_counts | rsv::drop_exactly(has_trailing_indel_run),
									deletion_counts | rsv::drop_exactly(has_trailing_indel_run),
									retval.mid
								)
							);
							
							libbio_assert_eq(ranges::size(rng), retval.mid.size());
							
							for (auto &&[ins_count, del_count, middle_part] : rng)
							{
								libbio_assert_lt(0, ins_count);
								libbio_assert_lt(0, del_count);
								fill_indel_run(ins_count, del_count, middle_part.lhs_indels);
								fill_indel_run(ins_count, del_count, middle_part.rhs_indels);
								middle_part.non_indel = *gen::arbitrary <non_indel_item>();
							}
						}
					}
					
					return retval;
				});
			});
		}
	};
	
	
	template <>
	struct Arbitrary <cigar_test_input>
	{
		static Gen <cigar_test_input> arbitrary()
		{
			return gen::exec([](cigar_sequence_pair &&seq_pair) -> cigar_test_input {
				// Introduce a change (“mutation”) to one of the CIGAR sequences.
				auto const mutate([](auto &part, mutation_position const pos, bool const should_modify_rhs, bool const mutation_is_del) -> bool {
					auto &vec(should_modify_rhs ? part.rhs_indels : part.lhs_indels);
				
					if (mutation_is_del)
					{
						switch (pos)
						{
							case mutation_position::left:
								++vec.leading_deletions;
								return true;
							
							case mutation_position::middle:
								if (vec.indels.empty())
									return false;
								++middle(vec.indels).del;
								return true;
							
							case mutation_position::right:
								if (vec.indels.empty())
									return false;
								++vec.indels.back().del;
								return true;
						}
					}
					else
					{
						switch (pos)
						{
							case mutation_position::left:
								if (vec.indels.empty())
									return false;
								++vec.indels.front().ins;
								return true;
							
							case mutation_position::middle:
								if (vec.indels.empty())
									return false;
								++middle(vec.indels).ins;
								return true;
							
							case mutation_position::right:
								++vec.trailing_insertions;
								return true;
						}
					}
				
					libbio_fail("Should not be reached");
					return false;
				});
				
				auto mutation_targets(*gen::container <std::vector <bool>>(seq_pair.mid.size() + 1, gen::arbitrary <bool>()));
				auto mutation_is_del(*gen::container <std::vector <bool>>(seq_pair.mid.size() + 1, gen::arbitrary <bool>()));
				
				std::vector <cigar_sequence_pair> non_matching_pairs;
				
				// Make changes in the CIGAR sequences.
				constexpr mutation_position const positions[]{mutation_position::left, mutation_position::middle, mutation_position::right};
				for (std::size_t i(0); i < seq_pair.mid.size(); ++i)
				{
					auto const should_modify_rhs(mutation_targets[i]);
					auto const current_mutation_is_del(mutation_is_del[i]);
				
					for (auto const pos : positions)
					{
						auto seq_pair_(seq_pair); // Copy.
						if (mutate(seq_pair_.mid[i], pos, should_modify_rhs, current_mutation_is_del))
							non_matching_pairs.emplace_back(std::move(seq_pair_));
					}
				}
			
				if (seq_pair.tail)
				{
					for (auto const pos : positions)
					{
						auto seq_pair_(seq_pair); // Copy.
						if (mutate(*seq_pair_.tail, pos, mutation_targets.back(), mutation_is_del.back()))
							non_matching_pairs.emplace_back(std::move(seq_pair_));
					}
				}
				
				return {std::move(seq_pair), std::move(non_matching_pairs), std::move(mutation_targets), std::move(mutation_is_del)};
			});
		}
	};
}


SCENARIO("cigar_eq() can handle simple non-matching sequences", "[cigar_eq]")
{
	GIVEN("a pair of non-matching CIGAR sequences")
	{
		std::vector const lhs{
			make_cigar_item(2, 'I'_cigar_operation),
			make_cigar_item(1, 'D'_cigar_operation),
			make_cigar_item(1, 'P'_cigar_operation)
		};
		
		std::vector const rhs{
			make_cigar_item(1, 'I'_cigar_operation),
			make_cigar_item(1, 'D'_cigar_operation),
			make_cigar_item(1, 'P'_cigar_operation)
		};
		
		WHEN("the sequences are compared")
		{
			THEN("the result is false")
			{
				CHECK(!panvc3::cigar_eq_seqan3(lhs, rhs));
			}
		}
	}
}


SCENARIO("cigar_eq() can handle padding", "[cigar_eq]")
{
	GIVEN("a pair of CIGAR sequences with an insertion-deletion-padding run")
	{
		std::vector const lhs{
			make_cigar_item(1, 'I'_cigar_operation),
			make_cigar_item(2, 'D'_cigar_operation),
			make_cigar_item(3, 'P'_cigar_operation),
			make_cigar_item(1, '='_cigar_operation),
		};
		
		std::vector const rhs{
			make_cigar_item(2, 'D'_cigar_operation),
			make_cigar_item(3, 'P'_cigar_operation),
			make_cigar_item(1, 'I'_cigar_operation),
			make_cigar_item(1, '='_cigar_operation),
		};
		
		WHEN("the sequences are compared")
		{
			THEN("the result is false when P ends a run")
			{
				CHECK(!panvc3::cigar_eq_seqan3 <false>(lhs, rhs));
			}
			
			THEN("the result is true when P does not end a run")
			{
				CHECK(panvc3::cigar_eq_seqan3 <true>(lhs, rhs));
			}
		}
	}
}


TEST_CASE("cigar_eq() can compare arbitrary matching sequences", "[cigar_eq]")
{
	rc::prop(
		"cigar_eq() works as expected",
		[](){
			
			cigar_test_input test_input(*rc::gen::arbitrary <cigar_test_input>());
			auto &seq_pair(test_input.matching_pair);
			
			{
				bool const sequence_pair_has_all_parts(
					!seq_pair.mid.empty() &&
					seq_pair.mid.front().has_all_parts() &&
					seq_pair.head &&
					seq_pair.tail &&
					seq_pair.tail->has_all_parts()
				);
				
				RC_CLASSIFY(sequence_pair_has_all_parts);
			}
			
			if (!panvc3::cigar_eq_seqan3(seq_pair.lhs(), seq_pair.rhs()))
				RC_FAIL("cigar_eq() returned for the matching pair");
			
			for (auto const &[idx, seq_pair] : rsv::enumerate(test_input.non_matching_pairs))
			{
				if (panvc3::cigar_eq_seqan3(seq_pair.lhs(), seq_pair.rhs()))
				{
					std::stringstream os;
					os << "cigar_eq() returned true for the non-matching pair " << idx;
					RC_FAIL(os.str());
				}
			}
			
			
			return true;
		}
	);
}
