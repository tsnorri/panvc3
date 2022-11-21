/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <array>
#include <catch2/catch.hpp>
#include <libbio/assert.hh>
#include <panvc3/rewrite_cigar.hh>
#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>


namespace lb	= libbio;
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
	
	
	void rewrite_cigar_and_check(
		std::string const &query,
		std::size_t const src_pos,
		std::string const &dst,
		std::size_t const expected_dst_pos,
		panvc3::sequence_entry_pair const &seq_entry_pair,
		std::span <seqan3::cigar const> const cigar,
		std::span <seqan3::cigar const> const expected_cigar
	)
	{
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
	}
	
	
	inline void rewrite_cigar_and_check(
		std::string const &query,
		std::size_t const src_pos,
		std::size_t const expected_dst_pos,
		input_fixture const &input,
		std::span <seqan3::cigar const> const cigar,
		std::span <seqan3::cigar const> const expected_cigar
	)
	{
		rewrite_cigar_and_check(query, src_pos, input.dst_without_gaps, expected_dst_pos, input.seq_entry_pair, cigar, expected_cigar);
	}
}


namespace {
	input_fixture const &simple_input_fixture_1()
	{
		static input_fixture retval{
			"GAT-ACA",
			"GATTACA"
		};
		return retval;
	}
}


// FIXME: All the test cases do basically the same thing. Instead of having separate scenarios, the input could be read from a configuration file and passed to a parsing function.
SCENARIO("rewrite_cigar() can handle deletions in source sequence", "[rewrite_cigar]")
{
	GIVEN("an alignment with a gap in the source and an aligned segment")
	{
		std::string const query("TA");
		auto const cigar(std::to_array <seqan3::cigar>({
			{2, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, '='_cigar_operation},
			{1, 'D'_cigar_operation},
			{1, '='_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 2, 2, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can handle a no-op in the CIGAR", "[rewrite_cigar]")
{
	GIVEN("an alignment with no-ops (“H” or “P”)")
	{
		std::string const query("TA");
		auto const cigar(std::to_array <seqan3::cigar>({
			{1, 'H'_cigar_operation},
			{2, 'M'_cigar_operation},
			{1, 'P'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, 'H'_cigar_operation},
			{1, '='_cigar_operation},
			{1, 'D'_cigar_operation},
			{1, '='_cigar_operation},
			{1, 'P'_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 2, 2, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can handle an insertion", "[rewrite_cigar]")
{
	GIVEN("an alignment with an insertion")
	{
		std::string const query("CC");
		auto const cigar(std::to_array <seqan3::cigar>({
			{2, 'I'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{2, 'I'_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 1, 1, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can handle an insertion after a gap in the source", "[rewrite_cigar]")
{	
	GIVEN("an alignment with an insertion after a gap in the source")
	{
		std::string const query("GG");
		auto const cigar(std::to_array <seqan3::cigar>({
			{2, 'I'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{2, 'I'_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 4, 5, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can handle soft clipping", "[rewrite_cigar]")
{
	GIVEN("an alignment with soft clipping")
	{
		std::string const query("CTAC");
		auto const cigar(std::to_array <seqan3::cigar>({
			{1, 'S'_cigar_operation},
			{2, 'M'_cigar_operation},
			{1, 'S'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, 'S'_cigar_operation},
			{1, '='_cigar_operation},
			{1, 'D'_cigar_operation},
			{1, '='_cigar_operation},
			{1, 'S'_cigar_operation}
		}));
	
		rewrite_cigar_and_check(query, 2, 2, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can handle soft clipping just before a deletion in the source", "[rewrite_cigar]")
{
	GIVEN("an alignment with soft clipping")
	{
		std::string const query("TC");
		auto const cigar(std::to_array <seqan3::cigar>({
			{1, 'M'_cigar_operation},
			{1, 'S'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, '='_cigar_operation},
			{1, 'S'_cigar_operation}
		}));
	
		rewrite_cigar_and_check(query, 2, 2, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can handle soft clipping just after a deletion in the source", "[rewrite_cigar]")
{
	GIVEN("an alignment with soft clipping")
	{
		std::string const query("CA");
		auto const cigar(std::to_array <seqan3::cigar>({
			{1, 'S'_cigar_operation},
			{1, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, 'S'_cigar_operation},
			{1, '='_cigar_operation}
		}));
	
		rewrite_cigar_and_check(query, 3, 4, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can handle deletion in the query", "[rewrite_cigar]")
{
	GIVEN("an alignment with a deletion")
	{
		std::string const query("T");
		auto const cigar(std::to_array <seqan3::cigar>({
			{1, 'D'_cigar_operation},
			{1, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, 'D'_cigar_operation},
			{1, '='_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 1, 1, simple_input_fixture_1(), cigar, expected_cigar);
	}
}


namespace {
	input_fixture const &simple_input_fixture_2()
	{
		static input_fixture retval{
			"GATTACA",
			"GAT-ACA"
		};
		return retval;
	}
}


SCENARIO("rewrite_cigar() can remap a match to an insertion if the destination has a deletion", "[rewrite_cigar]")
{
	GIVEN("an alignment with a match")
	{
		std::string const query("TTA");
		auto const cigar(std::to_array <seqan3::cigar>({
			{3, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, '='_cigar_operation},
			{1, 'I'_cigar_operation},
			{1, '='_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 2, 2, simple_input_fixture_2(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can remove a deletion if the destination has a deletion", "[rewrite_cigar]")
{
	GIVEN("an alignment with a match")
	{
		std::string const query("TA");
		auto const cigar(std::to_array <seqan3::cigar>({
			{1, 'M'_cigar_operation},
			{1, 'D'_cigar_operation},
			{1, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{2, '='_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 2, 2, simple_input_fixture_2(), cigar, expected_cigar);
	}
}


namespace {
	
	// Everything const.
	input_fixture const &simple_input_fixture_3()
	{
		static input_fixture retval{
		//	 01234567890123
			"GATTACAGATTACA",
			"GAT-ACAG--T-CA"
		};
		return retval;
	}
}


SCENARIO("rewrite_cigar() can handle an alignment that overlaps with a long, non-contiguous deletion in the destination", "[rewrite_cigar]")
{
	GIVEN("an alignment with a match")
	{
		std::string const query("TTACAGATTAC");
		auto const cigar(std::to_array <seqan3::cigar>({
			{11, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, '='_cigar_operation},
			{1, 'I'_cigar_operation},
			{4, '='_cigar_operation},
			{2, 'I'_cigar_operation},
			{1, '='_cigar_operation},
			{1, 'I'_cigar_operation},
			{1, '='_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 2, 2, simple_input_fixture_3(), cigar, expected_cigar);
	}
}


SCENARIO("rewrite_cigar() can convert the co-ordinate after multiple deletions", "[rewrite_cigar]")
{
	GIVEN("an alignment with a match")
	{
		std::string const query("CA");
		auto const cigar(std::to_array <seqan3::cigar>({
			{2, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{2, '='_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 12, 8, simple_input_fixture_3(), cigar, expected_cigar);
	}
}


namespace {
	
	// Everything const.
	input_fixture const &badly_aligned_input_fixture_1()
	{
		static input_fixture retval{
		//   0123456789012345
			"A-C-G-T-A-C-G-T-",
			"-G-T-A-C-G-T-A-C"
		};
		return retval;
	}
}


SCENARIO("rewrite_cigar() can handle badly aligned sections", "[rewrite_cigar]")
{
	GIVEN("an alignment with a match")
	{
		std::string const query("ACG");
		auto const cigar(std::to_array <seqan3::cigar>({
			{3, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, 'I'_cigar_operation},
			{1, 'D'_cigar_operation},
			{1, 'I'_cigar_operation},
			{1, 'D'_cigar_operation},
			{1, 'I'_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 4, 4, badly_aligned_input_fixture_1(), cigar, expected_cigar);
	}
}


namespace {
	input_fixture const &very_simple_input_fixture_1()
	{
		static input_fixture retval(
			"A-C",
			"ATC"
		);
		return retval;
	}
}


SCENARIO("rewrite_cigar() can handle padding next to deletion", "[rewrite_cigar]")
{
	GIVEN("a simple alignment")
	{
		std::string const query("AC");
		auto const cigar(std::to_array <seqan3::cigar>({
			{1, 'M'_cigar_operation},
			{1, 'P'_cigar_operation},
			{1, 'M'_cigar_operation}
		}));
		auto const expected_cigar(std::to_array <seqan3::cigar>({
			{1, '='_cigar_operation},
			{1, 'D'_cigar_operation},
			{1, 'P'_cigar_operation},
			{1, '='_cigar_operation}
		}));
		
		rewrite_cigar_and_check(query, 0, 0, very_simple_input_fixture_1(), cigar, expected_cigar);
	}
}
