/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <libbio/markov_chain.hh>
#include <libbio/markov_chain_rapidcheck.hh>
#include <libbio/utility/is_equal.hh>
#include <libbio/utility/is_lt.hh>
#include <panvc3/indel_run_checker.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/iterator/insert_iterators.hpp>	// ranges::back_inserter
#include <range/v3/iterator/stream_iterators.hpp>	// ranges::make_ostream_joiner()
#include <range/v3/view/tail.hpp>					// ranges::views::tail
#include <rapidcheck.h>
#include <rapidcheck/catch.h>						// rc::prop
#include "test_additions.hh"

namespace lb		= libbio;
namespace mcs		= libbio::markov_chains;
namespace rsv		= ranges::views;
namespace tests		= panvc3::tests;

using seqan3::operator""_cigar_operation;


namespace {
	
	constexpr static std::uint8_t cigar_op_count{9}; // M, I, D, N, S, H, P , =, and X.
	
	template <char t_op> // seqan3::cigar::operation is non-structural (according to GCC) and hence cannot be used as a non-type template parameter.
	struct cigar_op_specific
	{
		constexpr static auto const operation{seqan3::operator""_cigar_operation(t_op)};
		constexpr static std::array const allowed_operations{operation};
		constexpr static double const probability{1.0 / cigar_op_count};
		constexpr static double const probability_1{1.0 / (cigar_op_count - 1)};
	};
	
	typedef cigar_op_specific <'D'>	cigar_op_d;
	typedef cigar_op_specific <'I'>	cigar_op_i;
	
	struct cigar_op_indel
	{
		constexpr static std::array const allowed_operations{
			'D'_cigar_operation,
			'I'_cigar_operation
		};
		
		constexpr static double const probability{1.0 * allowed_operations.size() / cigar_op_count};
		constexpr static double const probability_1{1.0 * allowed_operations.size() / (cigar_op_count - 1)};
	};
	
	struct cigar_op_boundary
	{
		constexpr static std::array const allowed_operations{
			'M'_cigar_operation,
			'N'_cigar_operation,
			'S'_cigar_operation,
			'H'_cigar_operation,
			'P'_cigar_operation,
			'='_cigar_operation,
			'X'_cigar_operation
		};
		
		constexpr static double const probability{1.0 * allowed_operations.size() / cigar_op_count};
		constexpr static double const probability_1{1.0 * allowed_operations.size() / (cigar_op_count - 1)};
	};
	
	
	struct empty_cigar {};
	
	// I don’t know of seqan3::cigar’s constructor that would be easy to use with RapidCheck,
	// so we wrap it.
	template <typename t_op>
	struct cigar_
	{
		seqan3::cigar	value{};
		
		cigar_(seqan3::cigar const value_):
			value(value_)
		{
		}
		
		cigar_(seqan3::cigar::operation const op, panvc3::cigar_count_type const count)
		{
			value = op;
			value = count;
		}
		
		// FIXME: add assertions that the operation in value matches t_op?
		
		constexpr static auto const probability{t_op::probability};
		constexpr static auto const probability_1{t_op::probability_1};
	};
	
	typedef cigar_ <cigar_op_d>			cigar_d;
	typedef cigar_ <cigar_op_i>			cigar_i;
	typedef cigar_ <cigar_op_indel>		cigar_indel;
	typedef cigar_ <cigar_op_boundary>	cigar_boundary;
	
	// For converting from either cigar_ or empty_cigar;
	// needed by mcs::chain so that the resulting vector can have non-polymorphic objects.
	struct cigar
	{
		seqan3::cigar	value{};
		
		template <typename t_op>
		cigar(cigar_ <t_op> const cc):
			value(cc.value)
		{
		}
		
		cigar(empty_cigar const) {}
	};
	
	
	enum class run_type
	{
		no_op,
		non_indel,
		indel,
		boundary
	};
	
	struct run_base
	{
		panvc3::cigar_vector content{};
		
	protected:
		~run_base() {} // Destruction only allowed via subclasses.
	};
	
	struct no_op_run {}; // No need to inherit anything.
	
	struct indel_run : public run_base
	{
		typedef cigar_indel		cigar_type;
		constexpr static double scale{0.25};
		run_type type() const { return run_type::indel; }
	};
	
	struct boundary_run : public run_base
	{
		typedef cigar_boundary	cigar_type;
		constexpr static double scale{0.1};
		run_type type() const { return run_type::boundary; }
	};
	
	// Since we allow insertions and deletions but not next to each other,
	// we need to use a Markov chain.
	struct non_indel_run : public run_base
	{
		run_type type() const { return run_type::non_indel; }
	};
	
	// For converting from run_base.
	struct run
	{
		panvc3::cigar_vector	content{};
		run_type				type{run_type::no_op};
		
		template <typename t_run>
		run(t_run &&run_) requires std::is_base_of_v <run_base, t_run>:
			content(std::move(run_.content)),
			type(run_.type())
		{
		}
		
		run(no_op_run const) {}
	};
	
	
	std::ostream &operator<<(std::ostream &os, indel_run const &rr)
	{
		os << rr.content.size() << 'I';
		return os;
	}
	
	
	std::ostream &operator<<(std::ostream &os, boundary_run const &rr)
	{
		os << rr.content.size() << 'b';
		return os;
	}
	
	
	std::ostream &operator<<(std::ostream &os, non_indel_run const &rr)
	{
		os << rr.content.size() << 'i';
		return os;
	}
	
	
	std::ostream &operator<<(std::ostream &os, run const &rr)
	{
		os << rr.content.size();
		switch (rr.type)
		{
			case run_type::no_op:
				os << 'n';
				break;
			case run_type::non_indel:
				os << 'i';
				break;
			case run_type::indel:
				os << 'I';
				break;
			case run_type::boundary:
				os << 'b';
				break;
		}
		return os;
	}
	
	
	std::ostream &operator<<(std::ostream &os, std::vector <run> const &rv)
	{
		ranges::copy(rv, ranges::make_ostream_joiner(os, ", "));
		return os;
	}
	
	
	template <typename t_run>
	rc::Gen <t_run> arbitrary_simple_run()
	{
		return rc::gen::scale(
			t_run::scale,
			rc::gen::construct <t_run>(
				rc::gen::container <panvc3::cigar_vector>( // Produce the vector of seqan3::cigar needed by t_run.
					rc::gen::map( // fmap from generator of cigar to generator of seqan3::cigar.
						rc::gen::arbitrary <typename t_run::cigar_type>(),
						[](auto const wrapped){ return wrapped.value; }
					)
				)
			)
		);
	}
	
	template <typename t_src, typename t_dst>
	using transition = mcs::transition <t_src, t_dst, t_dst::probability>;
	
	template <typename t_src, typename t_dst>
	using transition_1 = mcs::transition <t_src, t_dst, t_dst::probability_1>;
	
	typedef mcs::chain <
		cigar,
		empty_cigar,
		mcs::transition_list <
			transition		<empty_cigar,		cigar_boundary>,
			transition		<empty_cigar,		cigar_d>,
			transition		<empty_cigar,		cigar_i>,
			transition		<cigar_boundary,	cigar_boundary>,
			transition		<cigar_boundary,	cigar_d>,
			transition		<cigar_boundary,	cigar_i>,
			transition_1	<cigar_d,			cigar_boundary>,
			transition_1	<cigar_d,			cigar_d>,
			transition_1	<cigar_i,			cigar_boundary>,
			transition_1	<cigar_i,			cigar_i>
		>
	> non_indel_markov_chain_type;
	
	
	struct test_input
	{
		typedef std::vector <run>	run_vector;
		
		run_vector					runs{};
		std::vector <std::size_t>	length_cs{}; // Cumulative sum of operation counts.
		
		test_input(run_vector &&runs_):
			runs(std::move(runs_))
		{
			length_cs.reserve(1 + runs.size());
			length_cs.push_back(0);
			
			std::size_t cs{};
			for (auto const &run : runs)
			{
				cs += run.content.size();
				length_cs.push_back(cs);
			}
		}
		
		void copy_operations(panvc3::cigar_vector &dst) const
		{
			dst.clear();
			for (auto const &run : runs)
				ranges::copy(run.content, ranges::back_inserter(dst));
		}
	};
	
	typedef mcs::chain <
		run,
		no_op_run,
		mcs::transition_list <
			mcs::transition <no_op_run,		indel_run,		0.5>,
			mcs::transition <no_op_run,		non_indel_run,	0.5>,
			mcs::transition <indel_run,		boundary_run,	1.0>,
			mcs::transition <non_indel_run,	boundary_run,	1.0>,
			mcs::transition <boundary_run,	indel_run,		0.5>,
			mcs::transition <boundary_run,	non_indel_run,	0.5>
		>
	> test_input_markov_chain_type;
}


namespace rc {
	
	template <>
	struct Arbitrary <empty_cigar>
	{
		static Gen <empty_cigar> arbitrary()
		{
			return gen::just(empty_cigar{});
		}
	};
	
	
	template <typename t_op>
	struct Arbitrary <cigar_ <t_op>>
	{
		static Gen <cigar_ <t_op>> arbitrary()
		{
			return gen::construct <cigar_ <t_op>>(
				gen::elementOf(t_op::allowed_operations),
				gen::inRange <panvc3::cigar_count_type>(1, 15)
			);
		}
	};
	
	
	template <>
	struct Arbitrary <no_op_run>
	{
		static Gen <no_op_run> arbitrary()
		{
			return gen::just(no_op_run{});
		}
	};
	
	
	template <>
	struct Arbitrary <boundary_run>
	{
		static Gen <boundary_run> arbitrary()
		{
			// Make sure that the run is non-empty.
			return gen::apply(
				[](auto &&run, cigar_boundary const boundary_op){
					if (run.content.empty())
						run.content.emplace_back(boundary_op.value);
					return std::move(run);
				},
				arbitrary_simple_run <boundary_run>(),
				gen::arbitrary <cigar_boundary>()
			);
		}
	};
	
	
	template <>
	struct Arbitrary <indel_run>
	{
		enum {
			HAS_NONE		= 0x0,
			HAS_INSERTION	= 0x1,
			HAS_DELETION	= 0x2,
			HAS_BOTH		= 0x3
		};
		
		static Gen <indel_run> arbitrary()
		{
			// We want to make sure that there are in fact both insertions and deletions in the run.
			return gen::apply(
				[](auto &&run, cigar_i const ins_op, cigar_d const del_op){
					// Check the sequence of operations.
					std::uint8_t status{HAS_NONE};
					for (auto const cc : run.content)
					{
						using seqan3::get;
						auto const operation(get <1>(cc));
						switch (operation.to_char())
						{
							case 'I':
								status |= HAS_INSERTION;
								break;
								
							case 'D':
								status |= HAS_DELETION;
								break;
								
							default:
								libbio_fail("Unexpected CIGAR operation");
						}
					}
					
					// Update the run if needed.
					switch (status)
					{
						case HAS_NONE:
							run.content.emplace_back(ins_op.value);
							run.content.emplace_back(del_op.value);
							break;
							
						case HAS_INSERTION:
							run.content.emplace_back(del_op.value);
							break;
							
						case HAS_DELETION:
							run.content.emplace_back(ins_op.value);
							break;
						
						case HAS_BOTH:
							break;
						
						default:
							libbio_fail("Unexpected status");
					}
					
					return std::move(run);
				},
				arbitrary_simple_run <indel_run>(),
				gen::arbitrary <cigar_i>(),
				gen::arbitrary <cigar_d>()
			);
		}
	};
	
	
	template <>
	struct Arbitrary <non_indel_run>
	{
		static Gen <non_indel_run> arbitrary()
		{
			// fmap the result of the Markov chain to panvc3::cigar_vector.
			// We can’t fmap the instances of cigar directly to seqan3::cigar b.c.
			// the Markov chain needs access to the previous element(s). 
			return gen::construct <non_indel_run>(
				gen::map(
					gen::scale(
						0.35,
						gen::arbitrary <non_indel_markov_chain_type>()
					),
					[](auto const &chain){
						panvc3::cigar_vector retval;
						if (chain.values.empty())
							return retval;
						
						retval.reserve(chain.values.size() - 1);
						ranges::copy(
							chain.values | rsv::tail | rsv::transform([](auto const vv){ return vv.value; }),
							ranges::back_inserter(retval)
						);
						
						return retval;
					}
				)
			);
		}
	};
	
	
	template <>
	struct Arbitrary <test_input>
	{
		static Gen <test_input> arbitrary()
		{
			return gen::construct <test_input>(
				gen::map(
					gen::arbitrary <test_input_markov_chain_type>(),
					[](auto &&chain){
						libbio_assert(!chain.values.empty());
						libbio_assert_eq(run_type::no_op, chain.values.front().type);
						chain.values.erase(chain.values.begin()); // Remove the no-op.
						return std::move(chain.values);
					}
				)
			);
		}
	};
}


TEST_CASE("indel_run_checker can process an arbitrary runs of CIGAR operations", "[indel_run_checker]")
{
	// FIXME: check the unaligned positions, too?
	rc::prop(
		"indel_run_checker() works as expected",
		[](test_input const &input){
			
			// Tag by indel run count.
			RC_TAG(ranges::count_if(input.runs, [](auto const &run){ return run_type::indel == run.type; }));
			
			panvc3::cigar_vector cigar_seq;
			input.copy_operations(cigar_seq);
			INFO("Runs:  " << input.runs);
			INFO("CIGAR: " << tests::to_readable(cigar_seq));
			
			panvc3::indel_run_checker checker;
			checker.reset(cigar_seq, 0);
			while (checker.find_next_range_for_realigning())
			{
				auto const realn_range(checker.cigar_realigned_range());
				auto const begin(input.length_cs.begin());
				auto const end(input.length_cs.end());
				
				auto const pos(std::distance(cigar_seq.cbegin(), realn_range.first));
				auto const end_pos(std::distance(cigar_seq.cbegin(), realn_range.second));
				INFO("Re-aln. range: [" << pos << ", " << end_pos << ')');
				
				// Find the correct run.
				auto const it(std::lower_bound(begin, end, pos));
				
				REQUIRE(it != end);
				auto const next(it + 1);
				
				// Compare to the position of the run.
				INFO("*it:     " << *it);
				INFO("pos:     " << pos);
				CHECK(lb::is_equal(*it, pos));
				REQUIRE(next != end);
				INFO("*next:   " << *next);
				INFO("end_pos: " << end_pos);
				CHECK(lb::is_equal(*next, end_pos));
				
				// Finally check that the run was an indel run.
				auto const run_idx(std::distance(begin, it));
				REQUIRE(lb::is_lt(run_idx, input.runs.size()));
				CHECK(input.runs[run_idx].type == run_type::indel);
			}
		}
	);
}
