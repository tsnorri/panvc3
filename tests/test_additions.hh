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
		os << get <0>(cc.op) << get <1>(cc.op).to_char();
		return os;
	}
	
	inline auto to_readable(std::span <seqan3::cigar const> span)
	{
		return span
		| ranges::views::transform([](auto const op) -> cigar {
			return {op};
		});
	}
	
	inline auto to_readable(std::vector <seqan3::cigar> const &vec)
	{
		std::span <seqan3::cigar const> span(vec.begin(), vec.end());
		return to_readable(span);
	}


	template <typename t_type>
	inline t_type copy_without_gaps(std::string_view const src)
	{
		auto rng(src | ranges::views::filter([](auto const cc){ return '-' != cc; }));
		return {ranges::begin(rng), ranges::end(rng)};
	}
	
	
	template <typename... t_base>
	struct overloaded : public t_base...
	{
		using t_base::operator()...;
	};
	
	
	template <typename t_type, typename t_tuple, std::size_t t_idx = 0>
	constexpr static auto const type_matches_tuple_element_v{
		std::is_same_v <t_type, std::tuple_element_t <t_idx, t_tuple>>
	};
}

#endif
