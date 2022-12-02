/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_TEST_ADDITIONS_HH
#define PANVC3_TEST_ADDITIONS_HH

#include <ostream>
#include <panvc3/cigar.hh>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <span>
#include <string_view>


namespace panvc3::tests {
	
	// Output helpers.
	struct cigar
	{
		seqan3::cigar	op;
		
		cigar(seqan3::cigar op_): op(op_) {}
	};
	
	inline std::ostream &operator<<(std::ostream &os, cigar const cc)
	{
		using seqan3::get;
		os << '(' << get <0>(cc.op) << ", " << get <1>(cc.op).to_char() << ')';
		return os;
	}
	
	inline auto to_readable(std::span <seqan3::cigar const> span)
	{
		return span
		| ranges::views::transform([](auto const op) -> cigar {
			return {op};
		});
	}


	inline std::string copy_without_gaps(std::string_view const src)
	{
		auto rng(src | ranges::views::filter([](auto const cc){ return '-' != cc; }));
		return {ranges::begin(rng), ranges::end(rng)};
	}
}

#endif
