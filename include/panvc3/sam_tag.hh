/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_SAM_TAG_HH
#define PANVC3_SAM_TAG_HH

#include <array>
#include <libbio/assert.hh>


namespace panvc3 {

	typedef std::uint16_t						seqan3_sam_tag_type;
	typedef std::vector <seqan3_sam_tag_type>	seqan3_sam_tag_vector;
	

	// Convert a SeqAn 3 SAM tag to a std::array <char, N> where 2 ≤ N.
	// I don't think SeqAn 3 has this utility function.
	// Compare to operator""_tag() in <seqan3/io/sam_file/sam_tag_dictionary.hpp>.
	template <std::size_t t_size>
	constexpr void from_tag(seqan3_sam_tag_type const val, std::array <char, t_size> &buffer)
	requires (2 <= t_size)
	{
		char const char0(val / 256); // Narrowed automatically when () (not {}) are used.
		char const char1(val % 256);
		std::get <0>(buffer) = char0;
		std::get <1>(buffer) = char1;
	}
	
	
	// Convert a std::array <char, 2> to a SeqAn 3 SAM tag.
	constexpr seqan3_sam_tag_type to_tag(std::array <char, 2> const &buffer)
	{
		// The tag needs to match /[A-Za-z][A-Za-z0-9]/ (SAMv1, Section 1.5 The alignment section: optional fields),
		// so the values will not be negative.
		libbio_always_assert_lte(std::get <0>(buffer), 127);
		libbio_always_assert_lte(std::get <1>(buffer), 127);
		seqan3_sam_tag_type retval(std::get <0>(buffer)); // Narrows when () (not {}) are used.
		retval *= 256;
		retval += std::get <1>(buffer);
		return retval;
	}
}

#endif
