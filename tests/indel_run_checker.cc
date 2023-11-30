/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <libbio/file_handling.hh>
#include <libbio/generic_parser.hh>
#include <libbio/generic_parser/cigar_field_seqan.hh>
#include <panvc3/indel_run_checker.hh>
#include <panvc3/range.hh>
#include <seqan3/alphabet/quality/phred42.hpp>
#include "test_additions.hh"

namespace lb		= libbio;
namespace lbp		= libbio::parsing;
namespace rsv		= ranges::views;
namespace tests		= panvc3::tests;
namespace tuples	= libbio::tuples;


namespace {
	
	template <typename... t_args>
	using parser = lbp::parser <
		lbp::traits::delimited <lbp::delimiter <'\t'>, lbp::delimiter <'\n'>>,
		t_args...
	>;
	
	struct indel_run_tag : public lbp::empty_tag {};
	struct non_indel_run_tag : public lbp::empty_tag {};
	struct cigar_input_tag : public lbp::empty_tag {};
	
	
	struct run
	{
		typedef panvc3::range	range;
		
		panvc3::cigar_vector	expected_cigar;
		range					expected_query_range;
		range					expected_ref_range;
		bool					is_indel_run{};
		
		
		template <std::size_t t_idx, typename t_tuple>
		static decltype(auto) get(t_tuple const &tup)
		{
			using std::get;
			using tuples::get;
			return get <t_idx>(tup);
		}
		
		
		template <std::size_t t_idx, typename t_tuple>
		static range make_range(t_tuple const &tup)
		{
			return {get <t_idx>(tup), get <1 + t_idx>(tup) - get <t_idx>(tup)};
		}
		
		
		template <typename t_tuple>
		run(indel_run_tag const, t_tuple const &tup):
			expected_cigar(get <1>(tup)),
			expected_query_range(make_range <2>(tup)),
			expected_ref_range(make_range <4>(tup)),
			is_indel_run(true)
		{
		}
		
		template <typename t_tuple>
		run(non_indel_run_tag const, t_tuple const &tup):
			expected_cigar(get <1>(tup))
		{
		}
	};
	
	
	struct conditional_field
	{
		template <typename t_caller>
		bool parse(t_caller &&caller) const
		{
			auto &range(caller.range());
			switch (*range.it)
			{
				case 'I':
					caller.read_delimiter();
					return caller.continue_parsing(indel_run_tag{});
					
				case 'N':
					caller.read_delimiter();
					return caller.continue_parsing(non_indel_run_tag{});
				
				case 'C':
					caller.read_delimiter();
					return caller.continue_parsing(cigar_input_tag{});
				
				default:
					throw lbp::parse_error_tpl(lbp::errors::unexpected_character(*range.it));
			}
		}
	};
}


SCENARIO("indel_run_checker can handle predefined inputs")
{
	lb::file_istream stream;
	lb::open_file_for_reading("indel_run_checker_inputs.tsv", stream);
	
	typedef lbp::traits::delimited <lbp::delimiter <'\t'>>							parser_traits;
	typedef lbp::traits::delimited <lbp::delimiter <'\t'>, lbp::delimiter <'\n'>>	conditional_parser_traits;
	
	typedef lbp::parser <
		parser_traits,
		lbp::fields::make_conditional <
			conditional_field,
			conditional_parser_traits,
			lbp::fields::option <
				non_indel_run_tag,
				lbp::fields::cigar
			>,
			lbp::fields::option <
				indel_run_tag,
				lbp::fields::cigar,
				lbp::fields::numeric <std::uint32_t>,
				lbp::fields::numeric <std::uint32_t>,
				lbp::fields::numeric <std::uint32_t>,
				lbp::fields::numeric <std::uint32_t>
			>,
			lbp::fields::option <
				cigar_input_tag,
				lbp::fields::cigar,
				lbp::fields::boolean,
				lbp::fields::boolean,
				lbp::fields::numeric <std::uint32_t>,
				lbp::fields::numeric <std::uint32_t>,
				lbp::fields::text <>
			>
		>
	> parser_type;
	
	typedef seqan3::phred42	quality_alphabet;
	
	std::istreambuf_iterator it(stream);
	std::istreambuf_iterator <char> const sentinel;
	auto range(lbp::make_range(it, sentinel));
	
	parser_type parser;
	parser_type::record_type rec;
	parser_type::buffer_type buffer;
	
	std::vector <run> expected_runs;
	parser.parse_all(range, rec, buffer, tests::overloaded{
		[&expected_runs]<typename t_tuple>(t_tuple const &tup) requires tests::type_matches_tuple_element_v <indel_run_tag, t_tuple> {
			expected_runs.emplace_back(indel_run_tag{}, tup);
		},
		[&expected_runs]<typename t_tuple>(t_tuple const &tup) requires tests::type_matches_tuple_element_v <non_indel_run_tag, t_tuple> {
			expected_runs.emplace_back(non_indel_run_tag{}, tup);
		},
		[&expected_runs]<typename t_tuple>(t_tuple const &tup) requires tests::type_matches_tuple_element_v <cigar_input_tag, t_tuple> {
			auto const &[tag, cigar_seq, has_preceding_, has_tail, expected_query_pos, expected_ref_pos, section_name] = tup;
			
			// cigar_input_tag marks the end of input for one test.
			SECTION(section_name)
			{
				panvc3::indel_run_checker checker;
				checker.reset(cigar_seq, 0);
				
				std::size_t idx{};
				bool has_preceding(has_preceding_); // Convert the wrapper to bool.
				auto cigar_begin(cigar_seq.begin());
				
				// Run the checker.
				while (checker.find_next_range_for_realigning())
				{
					auto const realn_range(checker.cigar_realigned_range());
					
					// Check the preceding non-indel run if needed.
					if (has_preceding)
					{
						REQUIRE(idx < expected_runs.size());
						auto const &run(expected_runs[idx]);
						auto const next_begin(realn_range.first);
						CHECK(!run.is_indel_run);
						INFO("expected_cigar: " << tests::to_readable(run.expected_cigar));
						INFO("actual_cigar:   " << tests::to_readable(std::span(cigar_begin, next_begin)));
						CHECK(std::equal(run.expected_cigar.begin(), run.expected_cigar.end(), cigar_begin, next_begin));
						cigar_begin = next_begin;
						++idx;
					}
					
					// Check the reported indel run.
					REQUIRE(idx < expected_runs.size());
					auto const &run(expected_runs[idx]);
					auto const next_begin(realn_range.second);
					CHECK(run.is_indel_run);
					INFO("expected_cigar: " << tests::to_readable(run.expected_cigar));
					INFO("actual_cigar:   " << tests::to_readable(std::span(cigar_begin, next_begin)));
					CHECK(std::equal(run.expected_cigar.begin(), run.expected_cigar.end(), cigar_begin, next_begin));
					
					auto const query_range(checker.query_range());
					auto const ref_range(checker.reference_range());
					CHECK(run.expected_query_range == query_range);
					CHECK(run.expected_ref_range == ref_range);
					
					cigar_begin = next_begin;
					has_preceding = true; // There must be a non-indel run before the next indel run.
					++idx;
				}
				
				// If there is a succeeding non-indel run, check it, too.
				if (has_tail)
				{
					REQUIRE(idx < expected_runs.size());
					auto const &run(expected_runs[idx]);
					CHECK(!run.is_indel_run);
					INFO("expected_cigar: " << tests::to_readable(run.expected_cigar));
					INFO("actual_cigar:   " << tests::to_readable(std::span(cigar_begin, cigar_seq.end())));
					CHECK(std::equal(run.expected_cigar.begin(), run.expected_cigar.end(), cigar_begin, cigar_seq.end()));
					++idx;
				}
				
				// Verify that everything was processed.
				CHECK(idx == expected_runs.size());
				CHECK(expected_query_pos == checker.query_position());
				CHECK(expected_ref_pos == checker.reference_position());
			}
			
			// Prepare for next test.
			expected_runs.clear();
		}
	});
}
