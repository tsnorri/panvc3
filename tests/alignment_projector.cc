/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <libbio/file_handling.hh>
#include <libbio/generic_parser.hh>
#include <libbio/generic_parser/cigar_field_seqan.hh>
#include <libbio/tuple/slice.hh>
#include <panvc3/alignment_projector.hh>
#include <panvc3/cigar_eq.hh>
#include "test_additions.hh"

namespace lb	= libbio;
namespace lbp	= libbio::parsing;
namespace rsv	= ranges::views;
namespace tests	= panvc3::tests;


namespace {
	
	struct input_fixture_
	{
		std::string					src;
		std::string 				dst;
		std::vector <char>			dst_without_gaps;
		std::string					given;
		panvc3::sequence_entry_pair	seq_entry_pair;
		
		input_fixture_() = default;
		
		template <typename t_src, typename t_dst, typename t_given>
		input_fixture_(
			t_src &&src_,
			t_dst &&dst_,
			t_given &&given_
		):
			src(std::forward <t_src>(src_)),
			dst(std::forward <t_dst>(dst_)),
			dst_without_gaps(tests::copy_without_gaps <std::vector <char>>(dst)),
			given(std::forward <t_given>(given_))
		{
			panvc3::make_sequence_entry_pair(src, dst, seq_entry_pair);
		}
	};
	
	struct input_fixture : public input_fixture_
	{
		using input_fixture_::input_fixture_;
		
		input_fixture &operator=(input_fixture &&other) &
		{
			input_fixture_::operator=(std::move(other));
			seq_entry_pair.first.fix_rank_select_pointers();
			seq_entry_pair.second.fix_rank_select_pointers();
			return *this;
		}
	};
	
	
	template <typename t_cigar, typename t_query_seq, typename t_dst_seq>
	std::size_t project_alignment(
		panvc3::alignment_projector_seqan3 &projector,
		std::size_t const src_pos,							// Position in the source reference.
		panvc3::sequence_entry_pair const &seq_entry_pair,
		t_query_seq const &query_seq,						// Typically std::vector <seqan3::dna5>
		std::vector <seqan3::phred42> const &quality_seq,
		t_cigar const &cigar_seq,
		t_dst_seq const &dst_seq
	)
	{
		return projector.project_alignment(
			src_pos,
			std::get <0>(seq_entry_pair),
			std::get <1>(seq_entry_pair),
			dst_seq,
			query_seq,
			cigar_seq,
			quality_seq,
			0,
			-1
		);
	}
	
	
	template <typename t_alphabet, typename t_character_filter = lbp::filters::no_op>
	struct seqan_alphabet_field
	{
		typedef t_alphabet					alphabet_type;
		typedef std::vector <t_alphabet>	vector_type;
		
		template <bool t_should_copy>
		using value_type = vector_type;
		

		constexpr void clear_value(vector_type &dst) const { dst.clear(); }

		
		template <typename t_delimiter, lbp::field_position t_field_position = lbp::field_position::middle_, typename t_range>
		constexpr inline bool parse(t_range &range, vector_type &dst) const
		{
			struct helper
			{
				vector_type &dst;
				
				void clear() { dst.clear(); }
				
				void handle_character(char const cc)
				{
					alphabet_type nt;
					nt.assign_char(cc);
					dst.push_back(nt);
				}
			};
			
			helper hh{dst};
			return lbp::fields::parse_sequential <t_delimiter, t_character_filter, t_field_position, std::string_view>(range, hh);
		}
	};
	
	
	template <typename... t_args>
	using parser = lbp::parser <
		lbp::traits::delimited <lbp::delimiter <'\t'>, lbp::delimiter <'\n'>>,
		t_args...
	>;
	
	
	struct sequence_pair_tag : public lbp::empty_tag {};
	struct query_tag : public lbp::empty_tag {};
	
	struct conditional_field
	{
		template <typename t_caller>
		bool parse(t_caller &&caller) const
		{
			auto &range(caller.range());
			switch (*range.it)
			{
				case 'S':
					caller.read_delimiter();
					return caller.continue_parsing(sequence_pair_tag{});
					
				case 'Q':
					caller.read_delimiter();
					return caller.continue_parsing(query_tag{});
					
				default:
					throw lbp::parse_error_tpl(lbp::errors::unexpected_character(*range.it));
			}
		}
	};
}


// FIXME: this is essentially the same as the test in rewrite_cigar.cc, so they could be combined.
SCENARIO("alignment_projector can re-align predefined inputs")
{
	lb::file_istream stream;
	lb::open_file_for_reading("alignment_projector_inputs.tsv", stream);
	
	typedef lbp::traits::delimited <lbp::delimiter <'\t'>>							parser_traits;
	typedef lbp::traits::delimited <lbp::delimiter <'\t'>, lbp::delimiter <'\n'>>	conditional_parser_traits;
	
	typedef lbp::parser <
		parser_traits,
		lbp::fields::make_conditional <
			conditional_field,
			conditional_parser_traits,
			lbp::fields::option <
				sequence_pair_tag,
				lbp::fields::text <>,
				lbp::fields::text <>,
				lbp::fields::text <>
			>,
			lbp::fields::option <
				query_tag,
				seqan_alphabet_field <seqan3::dna5>,
				lbp::fields::numeric <std::uint16_t>,
				lbp::fields::numeric <std::uint16_t>,
				lbp::fields::cigar,
				lbp::fields::cigar,
				lbp::fields::text <>, 
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
	std::vector <quality_alphabet> quality_seq;
	
	input_fixture fixture;
	parser.parse_all(range, rec, buffer, tests::overloaded{
		[&fixture]<typename t_tuple>(t_tuple const &tup) requires tests::type_matches_tuple_element_v <sequence_pair_tag, t_tuple> {
			auto const &[tag, inputs_given, src_seq, dst_seq] = tup;
			
			// Update the fixture.
			fixture = input_fixture(src_seq, dst_seq, inputs_given);
		},
		[&fixture, &quality_seq]<typename t_tuple>(t_tuple const &tup) requires tests::type_matches_tuple_element_v <query_tag, t_tuple> {
			auto const &[tag, query, src_pos, expected_dst_pos, cigar, expected_cigar, section_name, query_given] = tup;
			
			SECTION(section_name)
			{
				GIVEN(fixture.given)
				{
					AND_GIVEN(query_given)
					{
						auto const &seq_entry_pair(fixture.seq_entry_pair);
						auto const &dst(fixture.dst_without_gaps);
						
						WHEN("the segment is re-aligned")
						{
							quality_seq.clear();
							quality_seq.resize(query.size(), panvc3::max_letter <quality_alphabet>());
							
							panvc3::alignment_projector_seqan3 projector;
							auto const actual_dst_pos(project_alignment(projector, src_pos, seq_entry_pair, query, quality_seq, cigar, dst));
							
							THEN("the new position and the rewritten CIGAR match the expected")
							{
								auto const &actual_cigar(projector.alignment());
								
								INFO("expected_dst_pos: " << expected_dst_pos);
								INFO("actual_dst_pos:   " << actual_dst_pos);
								INFO("expected_cigar:   " << tests::to_readable(expected_cigar));
								INFO("actual_cigar:     " << tests::to_readable(actual_cigar));
								CHECK(expected_dst_pos == actual_dst_pos);
								CHECK(panvc3::cigar_eq_seqan3 <true>(expected_cigar, actual_cigar));
							}
						}
					}
				}
			}
		}
	});
	
	
}
