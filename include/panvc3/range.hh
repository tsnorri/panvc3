/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_RANGE_HH
#define PANVC3_RANGE_HH

#include <ostream>
#include <range/v3/view/slice.hpp>
#include <tuple>
#include <vector>


namespace panvc3 {
	
	struct range
	{
		typedef std::size_t	size_type;
		
		size_type location{};
		size_type length{};
		
		range() = default;
		
		range(size_type const location_, size_type const length_):
			location(location_),
			length(length_)
		{
		}
		
		auto slice() const { return ranges::views::slice(location, location + length); }
		void update_length(size_type const end_location) { length = end_location - location; }
		auto to_tuple() const { return std::make_tuple(location, length); }
		
		bool operator==(range const &other) const { return location == other.location && length == other.length; }
	};

	typedef std::vector <range>	range_vector;


	inline std::ostream &operator<<(std::ostream &os, range const &rr)
	{
		os << '[' << rr.location << ", " << (rr.location + rr.length) << ')';
		return os;
	}
}

#endif
