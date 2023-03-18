/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_UTILITY_HH
#define PANVC3_UTILITY_HH

#include <seqan3/utility/type_list/type_list.hpp>
#include <tuple>

namespace panvc3 {

	// I’m not sure how to get the types from seqan3::type_list, so here is a helper for that.
	// (Also I don’t know if there is syntax for generalising this for any template that takes a parameter pack.)
	
	template <typename t_type>
	struct type_list_to_tuple {};
	
	template <typename... t_types>
	struct type_list_to_tuple <seqan3::type_list <t_types...>>
	{
		typedef std::tuple <t_types...> type;
	};
	
	template <typename t_type>
	using type_list_to_tuple_t = type_list_to_tuple <t_type>::type;

	
	// FIXME: Move to libbio?
	template <typename t_type, typename t_proj, typename t_cmp = std::less <>>
	struct cmp_proj
	{
	private:
		template <typename, typename = void> struct is_transparent_t : public std::false_type {};

		template <typename t_type_>
		struct is_transparent_t <t_type_, std::void_t <typename t_type_::is_transparent>> : public t_type_::is_transparent {};

	public:
		typedef is_transparent_t <t_cmp>	is_transparent;


		constexpr bool operator()(t_type const &lhs, t_type const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(proj(lhs), proj(rhs));
		}

		template <typename t_other>
		constexpr bool operator()(t_type const &lhs, t_other const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(proj(lhs), rhs);
		}

		template <typename t_other>
		constexpr bool operator()(t_other const &lhs, t_type const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(lhs, proj(rhs));
		}
	};
	

	// FIXME: check that t_alphabet is one of SeqAn3’s alphabets.
	template <typename t_alphabet>
	constexpr inline t_alphabet max_letter()
	{
		t_alphabet retval;
		retval.assign_rank(t_alphabet::alphabet_size - 1);
		return retval;
	}
	
	
	template <typename t_value>
	constexpr bool is_power_of_2(t_value const val)
	{
		return 0 < val && !(val & (val - 1));
	}
}

#endif
