/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <array>
#include <catch2/catch.hpp>
#include <libbio/assert.hh>
#include <libbio/file_handling.hh>
#include <libbio/generic_parser.hh>
#include <libbio/utility/tuple_slice.hh>
#include <panvc3/rewrite_cigar.hh>
#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>


namespace lb	= libbio;
namespace lbp	= libbio::parsing;
namespace rsv	= ranges::views;

using seqan3::operator""_cigar_operation;


namespace {
	
	// Output helpers.
	struct cigar
	{
		seqan3::cigar	op;
		
		cigar(seqan3::cigar op_): op(op_) {}
	};
	
	std::ostream &operator<<(std::ostream &os, cigar const cc)
	{
		using seqan3::get;
		os << '(' << get <0>(cc.op) << ", " << get <1>(cc.op).to_char() << ')';
		return os;
	}
	
	auto to_readable(std::span <seqan3::cigar const> span)
	{
		return span
		| rsv::transform([](auto const op) -> cigar {
			return {op};
		});
	}


	std::string copy_without_gaps(std::string_view const src)
	{
		auto rng(src | rsv::filter([](auto const cc){ return '-' != cc; }));
		return {ranges::begin(rng), ranges::end(rng)};
	}
	
	
	struct input_fixture
	{
		std::string					src;
		std::string 				dst;
		std::string					dst_without_gaps;
		panvc3::sequence_entry_pair	seq_entry_pair;
		
		input_fixture() = default;
		
		template <typename t_src, typename t_dst>
		input_fixture(
			t_src &&src_,
			t_dst &&dst_
		):
			src(std::forward <t_src>(src_)),
			dst(std::forward <t_dst>(dst_)),
			dst_without_gaps(copy_without_gaps(dst))
		{
			panvc3::make_sequence_entry_pair(src, dst, seq_entry_pair);
		}
	};
	
	
	template <typename t_cigar, typename t_query_seq, typename t_dst_seq>
	std::size_t rewrite_cigar(
		std::size_t const src_pos,					// Position in the source reference.
		t_cigar const &cigar_seq,
		panvc3::sequence_entry_pair const &seq_entry_pair,
		t_query_seq const &query_seq,				// Typically std::vector <seqan3::dna5>
		t_dst_seq const &dst_seq,
		panvc3::cigar_buffer &destination
	)
	{
		return panvc3::rewrite_cigar(src_pos, cigar_seq, std::get <0>(seq_entry_pair), std::get <1>(seq_entry_pair), query_seq, dst_seq, destination);
	}
}


namespace {
	
	template <typename... t_args>
	using parser = lbp::parser <
		lbp::traits::delimited <lbp::delimiter <'\t'>, lbp::delimiter <'\n'>>,
		t_args...
	>;
	
	
	struct cigar_field
	{
		template <bool t_should_copy>
		using value_type = std::vector <seqan3::cigar>;
		
		
		template <typename t_range>
		constexpr inline seqan3::cigar parse_one(t_range &range) const
		{
			seqan3::cigar retval{};
			panvc3::cigar_count_type count{};
			char cigar_op{};
			lbp::fields::integer <panvc3::cigar_count_type> integer_field;
			lbp::fields::character character_field;
			
			integer_field.parse_value <lbp::field_position::middle_>(range, count);
			retval = count;
			
			character_field.parse_value <lbp::field_position::middle_>(range, cigar_op);
			switch (cigar_op)
			{
				case 'M':
					retval = 'M'_cigar_operation;
					break;
				case 'I':
					retval = 'I'_cigar_operation;
					break;
				case 'D':
					retval = 'D'_cigar_operation;
					break;
				case 'N':
					retval = 'N'_cigar_operation;
					break;
				case 'S':
					retval = 'S'_cigar_operation;
					break;
				case 'H':
					retval = 'H'_cigar_operation;
					break;
				case 'P':
					retval = 'P'_cigar_operation;
					break;
				case '=':
					retval = '='_cigar_operation;
					break;
				case 'X':
					retval = 'X'_cigar_operation;
					break;
				default:
					throw lbp::parse_error_tpl(lbp::errors::unexpected_character(cigar_op));
			}
			
			return retval;
		}
		
		
		template <typename t_delimiter, lbp::field_position t_field_position = lbp::field_position::middle_, typename t_range>
		constexpr inline bool parse(t_range &range, std::vector <seqan3::cigar> &dst) const
		{
			dst.clear();
			
			if constexpr (any(lbp::field_position::initial_ & t_field_position))
			{
				if (range.is_at_end())
					return false;
				goto continue_parsing;
			}
			
			while (!range.is_at_end())
			{
			continue_parsing:
				dst.emplace_back(parse_one(range));
				
				auto const cc(*range.it);
				if (t_delimiter::matches(cc))
				{
					++range.it;
					return true;
				}
			}
			
			if constexpr (t_field_position == lbp::field_position::final_)
			{
				if (range.is_at_end())
					return true;
			}
			else
			{
				if (range.is_at_end())
					throw lbp::parse_error_tpl(lbp::errors::unexpected_eof());
			}
			
			libbio_fail("Should not be reached");
			return false;
		}
	};
}


// Not particularly good scenario name but in Catch2, there does not seem to be a container for scenarios.
// A better option would be to have the scenario names in the input TSV.
SCENARIO("rewrite_cigar() can handle predefined inputs", "[rewrite_cigar]")
{
	lb::file_istream stream;
	lb::open_file_for_reading("rewrite_cigar_inputs.tsv", stream);
	
	typedef lbp::parser <
		lbp::traits::delimited <lbp::delimiter <'\t'>>,
		lbp::fields::character
	> record_type_column_parser_type;
	
	typedef parser <lbp::fields::text <>, lbp::fields::text <>, lbp::fields::text <>>	sequence_pair_parser_type;
	typedef parser <
		lbp::fields::text <>,
		lbp::fields::numeric <std::uint16_t>,
		lbp::fields::numeric <std::uint16_t>,
		cigar_field,
		cigar_field,
		lbp::fields::text <>, 
		lbp::fields::text <>
	> query_parser_type;
	
	std::istreambuf_iterator it(stream);
	std::istreambuf_iterator <char> const sentinel;
	auto range(lbp::make_range(it, sentinel));
	
	record_type_column_parser_type record_type_parser;
	sequence_pair_parser_type sequence_pair_parser;
	query_parser_type query_parser;
	
	record_type_column_parser_type::record_type record_type_tuple;
	sequence_pair_parser_type::record_type sequence_pair_tuple;
	query_parser_type::record_type query_tuple;
	
	// Initial record.
	{
		auto const status(record_type_parser.parse(range, record_type_tuple));
		CHECK(status);
		CHECK(std::get <0>(record_type_tuple) == 'S');
	}
	
	while (true)
	{
		// Parse the sequence pair record.
		{
			auto const status(sequence_pair_parser.parse(range, sequence_pair_tuple));
			REQUIRE(status);
		}
		
		// Set up the sequence fixture.
		input_fixture fixture(std::get <1>(sequence_pair_tuple), std::get <2>(sequence_pair_tuple));
		
		while (true)
		{
			// Parse the inputs for the fixture.
			{
				auto const status(record_type_parser.parse(range, record_type_tuple));
			
				// Stop if we got an EOF.
				if (!status)
					return;
			}
			
			auto const op(std::get <0>(record_type_tuple));
			switch (op)
			{
				case 'S':
					goto continue_outer_loop;
				
				case 'Q':
				{
					{
						auto const status(query_parser.parse(range, query_tuple));
						REQUIRE(status);
					}
				
					SECTION(std::get <5>(query_tuple))
					{
						GIVEN(std::get <0>(sequence_pair_tuple))
						{
							AND_GIVEN(std::get <6>(query_tuple))
							{
								std::apply(
									[&fixture]
									(auto const &query, auto const src_pos, auto const expected_dst_pos, auto const &cigar, auto const &expected_cigar){
										auto const &seq_entry_pair(fixture.seq_entry_pair);
										auto const &dst(fixture.dst_without_gaps);
										
										WHEN("the CIGAR string is rewritten")
										{
											panvc3::cigar_buffer buffer;
											auto const dst_pos(rewrite_cigar(src_pos, cigar, seq_entry_pair, query, dst, buffer));
											
											THEN("the new position and the rewritten CIGAR match the expected")
											{
												INFO("expected_dst_pos: " << expected_dst_pos);
												INFO("actual:           " << dst_pos);
												INFO("expected_cigar:   " << to_readable(expected_cigar));
												INFO("actual:           " << to_readable(buffer.operations()));
												CHECK(expected_dst_pos == dst_pos);
												CHECK(panvc3::cigar_eq <true>(expected_cigar, buffer.operations()));
											}
										}
									},
									lb::tuple_slice <0, 5, true>(query_tuple)
								);
							}
						}
					}
					
					break;
				}
			
				default:
					FAIL("Unexpected input");
			}
		}
		
	continue_outer_loop:
		;
	}
}
