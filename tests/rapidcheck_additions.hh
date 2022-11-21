/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_RAPIDCHECK_ADDITIONS_HH
#define PANVC3_RAPIDCHECK_ADDITIONS_HH

#include <rapidcheck.h>

// gen::inRange takes a half-open range but we would like to be able to generate e.g. UINT64_MAX.
namespace rc::gen {
	
	template <typename t_type>
	Gen <t_type> inClosedRange(t_type const min, t_type const max)
	{
		return [=](Random const &random, int size){
			if (max < min)
				throw GenerationFailure("Invalid range");
			
			auto const rangeSize(detail::scaleInteger(static_cast <Random::Number>(max) - static_cast <Random::Number>(min), size) + 1);
			auto const value(static_cast <t_type>((Random(random).next() % rangeSize) + min));
			assert(min <= value && value <= max);
			return shrinkable::shrinkRecur(value, [=](t_type const x){ return shrink::towards <t_type>(x, min); });
		};
	}
}

#endif