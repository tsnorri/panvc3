/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DNA11_ALPHABET_HH
#define PANVC3_DNA11_ALPHABET_HH

#include <seqan3/alphabet/nucleotide/dna5.hpp>

namespace panvc3 {
	
	constexpr inline seqan3::dna5 dna5_from_char(char const cc)
	{
		seqan3::dna5 retval;
		retval.assign_char(cc);
		return retval;
	}
	
	
	// Alphabet for representing soft-clipped bases with lowercase letters.
	// FIXME: this could be rewritten as a decorator.
	class dna11 : public seqan3::nucleotide_base <dna11, 11>
	{
	private:
		typedef seqan3::nucleotide_base <dna11, 11>	base_t;
		friend base_t;
		friend base_t::base_t;
		
	public:
		constexpr dna11()							noexcept = default;
		constexpr dna11(dna11 const &)				noexcept = default;
		constexpr dna11(dna11 &&)					noexcept = default;
		constexpr dna11 &operator=(dna11 const &)	noexcept = default;
		constexpr dna11 &operator=(dna11 &&)		noexcept = default;
		~dna11()									noexcept = default;
		
		/* implicit */ constexpr dna11(seqan3::dna5 const cc, bool const is_clipped = false) noexcept
		{
			assign_rank(cc.to_rank() + is_clipped * alphabet_size / 2);
		}
		
		/* implicit */ constexpr dna11(char const cc, bool const is_clipped = false) noexcept:
			dna11(dna5_from_char(cc), is_clipped)
		{
		}
		
		using base_t::base_t;
		
	protected:
		constexpr static char_type rank_to_char_table[alphabet_size]
		{
			'A', // 0
			'C', // 1
			'G', // 2
			'N', // 3
			'T', // 4
			'a', // 5
			'c', // 6
			'g', // 7
			'n', // 8
			't', // 9
			'~', // 10
		};
		
		constexpr static rank_type rank_complement_table[alphabet_size]{
			4,	// A
			2,	// C
			1,	// G
			3,	// N
			0,	// T
			9,	// a
			7,	// c
			6,	// g
			8,	// n
			5,	// t
			10	// ~
		};
		
		constexpr static rank_type rank_complement(rank_type const rank)
		{
			return rank_complement_table[rank];
		}
		
		constexpr static char_type rank_to_char(rank_type const rank)
		{
			return rank_to_char_table[rank];
		}
		
		constexpr static rank_type char_to_rank(char_type const chr)
		{
			typedef std::make_unsigned_t <char_type> index_t;
			return char_to_rank_table[static_cast <index_t>(chr)];
		}
		
		constexpr static std::array <rank_type, 256> char_to_rank_table
		{
			[]() constexpr {
				std::array <rank_type, 256> retval;
				retval.fill(3); // N’s rank.
				
				for (std::size_t rank(0); rank < alphabet_size; ++rank)
					retval[rank_to_char_table[rank]] = rank;
				
				// For the sake of completeness; we don’t currently handle RNA.
				retval['U'] = retval['T'];
				retval['u'] = retval['t'];
				
				return retval;
			}()
		};
	};
	
	
	typedef std::vector <dna11>	dna11_vector;
	
	
	constexpr dna11 operator""_dna11(char const cc) noexcept
	{
		return dna11{}.assign_char(cc);
	}
	
	
	std::ostream &operator<<(std::ostream &os, dna11 const cc)
	{
		os << cc.to_char();
		return os;
	}
}

#endif
