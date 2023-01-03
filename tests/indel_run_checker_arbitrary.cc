/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <libbio/fmap.hh>
#include <libbio/tuple/group_by.hh>
#include <libbio/tuple/filter.hh>
#include <libbio/tuple/for.hh>
#include <libbio/tuple/unique.hh>
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
namespace rsv		= ranges::views;
namespace tests		= panvc3::tests;
namespace tuples	= libbio::tuples;

using seqan3::operator""_cigar_operation;


namespace {
	
	template <typename, typename> struct callback_table_builder {};
	
	template <typename t_callback, typename... t_args>
	struct callback_table_builder <t_callback, std::tuple <t_args...>>
	{
		constexpr static std::array const fns{(&t_callback::template operator() <t_args>)...};
	};
	
	
	template <typename... t_args>
	struct aggregate : public t_args... {};
}


// FIXME: rewrite this s.t. defining the type still produces (ultimately) a constinit or constexpr Markov chain but the chain may also be built at run time.
namespace markov_chains {
	
	template <typename t_src, typename t_dst, double t_probability>
	struct transition
	{
		typedef t_src					source_type;
		typedef t_dst					destination_type;
		constexpr static double const	probability{t_probability};
	};
	
	
	template <typename... t_transitions>
	struct transition_list
	{
		typedef std::tuple <t_transitions...>	transitions_type;
	};
	
	
	template <typename t_initial_state, typename... t_others>
	struct transition_to_any_other
	{
		typedef std::tuple <t_others...>		non_initial_states_type;
		constexpr static inline auto const other_count{sizeof...(t_others)};
		constexpr static inline auto const other_count_1{other_count - 1};
		
		template <typename t_type>
		using from_initial_t = transition <t_initial_state, t_type, 1.0 / other_count>;
		
		template <typename t_src, typename t_dst>
		using from_other_t = transition <t_src, t_dst, 1.0 / other_count_1>;
		
		typedef tuples::map_t <non_initial_states_type, from_initial_t>	initial_transitions_type;
		
		// Determine the cross product, then remove pairs of the same type.
		// Finally make the transitions.
		typedef tuples::map <
			tuples::filter <
				tuples::cross_product_t <non_initial_states_type, non_initial_states_type, std::tuple>,
				tuples::negation <std::is_same>::with
			>,
			from_other_t
		>																non_initial_transitions_type;
		
		typedef tuples::cat_t <
			initial_transitions_type,
			non_initial_transitions_type
		>																transitions_type;
	};
	
	
	struct uses_runtime_polymorphism_trait {};
}


namespace markov_chains::detail {
	
	typedef std::size_t	node_type;
	
	template <typename t_transition> using transition_src_t = typename t_transition::source_type;
	template <typename t_transition> using transition_dst_t = typename t_transition::destination_type;
	
	template <typename, typename>
	struct transition_states_have_base {};
	
	template <typename t_base, typename... t_transition>
	struct transition_states_have_base <t_base, std::tuple <t_transition...>> :
		public std::bool_constant <
			(std::is_base_of_v <t_base, typename t_transition::source_type> && ...) &&
			(std::is_base_of_v <t_base, typename t_transition::destination_type> && ...)
		> {};
	
	template <typename t_base, typename t_tuple>
	constexpr inline bool const transition_states_have_base_v = transition_states_have_base <t_base, t_tuple>::value;
	
	
	// We would like to calculate the cumulative sum of probabilities for each equivalence class.
	// transitions_by_source_type has tuples that represent equivalence classes like
	// {src, {{src, dst, probability}, …}}, so we need to process them.
	template <typename t_transition, double t_cs>
	using update_transition_probability_t = transition <
		typename t_transition::source_type,
		typename t_transition::destination_type,
		t_transition::probability + t_cs
	>;
	
	template <typename t_transitions_ = tuples::empty, double t_cs = 0.0>
	struct transition_probability_cs_acc
	{
		typedef t_transitions_ transitions;
		constexpr static inline double const cs{t_cs};
	};
	
	template <
		typename t_acc,
		typename t_transition,
		typename t_updated_transition = update_transition_probability_t <t_transition, t_acc::cs>
	>
	using transition_probability_cs_fold_fn = transition_probability_cs_acc <
		tuples::append_t <typename t_acc::transitions, t_updated_transition>,	// Append to the transition list.
		t_updated_transition::probability										// Store the new probability.
	>;
	
	// Check that the sum of transition probabilities is 1.0.
	template <typename t_transitions>
	struct transition_probability_cs_check
	{
		// FIXME: compare with some epsilon since we calculate a sum of floating point values?
		static_assert(1.0 == tuples::last_t <t_transitions>::probability);
		typedef t_transitions transitions_type;
	};
	
	// Map each equivalence class of transitions s.t. the probabilities are changed to their cumulative sum.
	// The equivalence class representatives are discarded.
	template <
		typename t_eq_class,
		typename t_transitions = typename tuples::foldl_t <
			transition_probability_cs_fold_fn,
			transition_probability_cs_acc <>,
			tuples::second_t <t_eq_class>
		>::transitions,
		typename t_check = transition_probability_cs_check <t_transitions>
	>
	using transition_probability_cs_map_fn = typename t_check::transitions_type;
	
	struct transition_key
	{
		double		probability_threshold{};
		node_type	node{};
		
		constexpr auto as_tuple() const { return std::make_tuple(node, probability_threshold); }
		constexpr bool operator<(transition_key const &other) const { return as_tuple() < other.as_tuple(); }
		constexpr bool operator==(transition_key const &other) const { return as_tuple() == other.as_tuple(); }
	};
	
	template <typename t_nodes, typename t_transitions>
	struct transition_map_builder
	{
		typedef t_nodes										nodes_type;
		typedef std::pair <transition_key, node_type>		map_value_type;
		
		// Intermediate type for creating map_value_type.
		template <std::size_t t_src, std::size_t t_dst, double t_probability>
		struct transition_
		{
			constexpr map_value_type to_map_value() const { return {{t_probability, t_src}, t_dst}; }
			constexpr /* implicit */ operator map_value_type() const { return to_map_value(); }
		};
	
		// FIXME: this is really inefficient. Come up with another solution, e.g. store the indices in sorted vectors.
		template <typename t_transition>
		using to_transition__t = transition_ <
			tuples::first_index_of_v <nodes_type, typename t_transition::source_type>,
			tuples::first_index_of_v <nodes_type, typename t_transition::destination_type>,
			t_transition::probability
		>;
	
		constexpr static auto make_transition_map()
		{
			// Group the transitions by the source node.
			typedef typename tuples::group_by_type <
				t_transitions,
				transition_src_t
			>::keyed_type										transitions_by_source_type;
		
			// Build the transition table.
			typedef tuples::cat_with_t <
				tuples::map_t <
					transitions_by_source_type,
					detail::transition_probability_cs_map_fn
				>
			>													transitions_with_probability_cs_type;
		
			typedef tuples::map_t <
				transitions_with_probability_cs_type,
				to_transition__t
			>													intermediate_transitions_type;
		
			constexpr intermediate_transitions_type intermediates{};
			auto retval(lb::map_to_array(
				std::make_index_sequence <std::tuple_size_v <t_transitions>>{},
				[&intermediates](auto const idx) {
					return std::get <idx()>(intermediates).to_map_value();
				}
			));
		
			std::sort(retval.begin(), retval.end());
		
			// Make sure that the keys are distinct.
			libbio_assert_eq(retval.cend(), std::adjacent_find(retval.cbegin(), retval.cend(), [](auto const &lhs, auto const &rhs){
				return lhs.first == rhs.first;
			}));
		
			return retval;
		}
	};
}


namespace markov_chains {
	
	// Since there is potentially a very large number of paths through the corresponding directed (cyclic) graph,
	// we can make use of runtime polymorphism if needed. The user can specify the base class of the instantiated
	// objects, though. This should result in an aggregate type as of C++20. (In particular, all the non-static
	// data members are public.)
	template <typename t_base, typename t_initial_state, typename t_transitions, typename... t_traits>
	struct chain : public t_traits...
	{
	private:
		typedef aggregate <t_traits...>						traits_type;
		
	public:
		typedef t_initial_state								initial_state_type;
		typedef typename t_transitions::transitions_type	transitions_type;
		typedef detail::node_type							node_type;
		
		constexpr static auto const NODE_MAX{std::numeric_limits <node_type>::max()}; // Max. value of node_type.
		
		// Number the nodes.
		// This can be done easily by determining the unique types over the concatenation of the lists
		// of source and destination node types.
		typedef tuples::unique_t <
			tuples::cat_t <
				std::tuple <t_initial_state>,
				tuples::map_t <transitions_type, detail::transition_src_t>,
				tuples::map_t <transitions_type, detail::transition_dst_t>
			>
		>													nodes_type;
		
		constexpr static std::size_t const					node_limit{std::tuple_size_v <nodes_type>}; // The number of nodes in this chain.
		constexpr static bool const							uses_runtime_polymorphism{std::is_base_of_v <uses_runtime_polymorphism_trait, traits_type>};
		static_assert(!uses_runtime_polymorphism || std::is_base_of_v <t_base, t_initial_state>);
		static_assert(!uses_runtime_polymorphism || detail::transition_states_have_base_v <t_base, transitions_type>);
		
		typedef detail::transition_map_builder <
			nodes_type,
			transitions_type
		>													transition_map_builder_type;
		typedef transition_map_builder_type::map_value_type	transition_map_value_type;
		
	private:
		// The actual chain.
		constexpr static auto const			initial_state{tuples::first_index_of_v <nodes_type, t_initial_state>};
		constexpr static auto const			transitions{transition_map_builder_type::make_transition_map()};
		
	public:
		typedef std::vector <
			std::conditional_t <
				uses_runtime_polymorphism,
				std::unique_ptr <t_base>,
				t_base
			>
		>									values_type;
		
		values_type							values{};
		
		template <typename t_probabilities, typename t_visitor>
		static void visit_node_types(t_probabilities &&probabilities, t_visitor &&visitor)
		{
			typedef detail::transition_key transition_key;
			
			auto const &fns(callback_table_builder <t_visitor, nodes_type>::fns);
			node_type current_node{initial_state};
			
			libbio_assert_lt(current_node, fns.size());
			(visitor.*(fns[current_node]))();
			for (auto const &prob : probabilities)
			{
				transition_key const key{prob, current_node};
				auto const it(std::upper_bound(
					transitions.begin(),
					transitions.end(),
					key,
					[](transition_key const &key, transition_map_value_type const &val){
						return key < val.first;
					}
				));
				libbio_assert_neq(transitions.cend(), it);
				current_node = it->second;
				libbio_assert_lt(current_node, fns.size());
				(visitor.*(fns[current_node]))();
			}
		}
		
		// No user-declared or inherited constructors allowed, so we construct with this static function template.
		template <typename t_probabilities>
		static chain from_probabilities(t_probabilities &&probabilities)
		requires uses_runtime_polymorphism
		{
			chain retval;
			retval.values.reserve(1 + probabilities.size());
			visit_node_types(probabilities, [&retval]<typename t_type> {
				if constexpr (uses_runtime_polymorphism)
					retval.values.emplace_back(std::make_unique <t_type>());
				else
					retval.values.emplace_back(t_type{}); // Needs a converting constructor.
			});
			return retval;
		}
	};
}


namespace mcs = markov_chains;


namespace rc {
	
	template <typename t_base, typename t_initial_state, typename t_transitions>
	struct Arbitrary <mcs::chain <t_base, t_initial_state, t_transitions>>
	{
		typedef mcs::chain <t_base, t_initial_state, t_transitions>	chain_type;
		
		static Gen <chain_type> arbitrary()
		{
			return gen::withSize([](int const size){
				// I don’t think RapidCheck knows how to shrink the chain, so we’ll use gen::shrink.
				// Also it seems to be quite difficult to functionally map a collection (of generators to a generator
				// that produces a collection), so we’ll resort to gen::exec. (We could use gen::mapCat to produce the
				// vector of probabilities but then using gen::arbitrary with the nodes / states would not be possible.)
				return gen::shrink(
					gen::exec([size](){
						auto const probabilities(
							*gen::container <std::vector <double>>(
								size,
								gen::map(
									// We need to use the half-open range here b.c. std::upper_bound is applied to find
									// the correct transition.
									gen::inRange <std::uint32_t>(0.0, UINT32_MAX),
									[](auto const val){
										return 1.0 * val / UINT32_MAX;
									}
								)
							)
						);
						
						chain_type retval;
						retval.values.reserve(1 + size);
						chain_type::visit_node_types(probabilities, [&retval]<typename t_type> {
							if constexpr (chain_type::uses_runtime_polymorphism)
							{
								// It seems that there is no other way to combine gen::arbitrary and operator new.
								retval.values.emplace_back(
									std::make_unique <t_type>(
										*gen::arbitrary <t_type>()
									)
								);
							}
							else
							{
								retval.values.emplace_back(*gen::arbitrary <t_type>()); // Needs a converting constructor.
							}
						});
						return retval;
					}),
					[](chain_type &&mc){
						return seq::takeWhile(
							seq::iterate(std::move(mc), [](chain_type &&mc){ mc.values.pop_back(); return std::move(mc); }),
							[](chain_type const &mc){ return !mc.values.empty(); }
						);
					}
				);
			});
		}
	};
}


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
