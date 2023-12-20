/*
 * Copyright (c) 2022–2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_CIGAR_EQ_HH
#define PANVC3_CIGAR_EQ_HH

#include <panvc3/cigar_adapter.hh>
#include <range/v3/range/access.hpp>


namespace panvc3 {
	
	template <bool t_should_count_padding = false, typename t_lhs_rng, typename t_rhs_rng, typename t_adapter>
	bool cigar_eq_tpl(t_lhs_rng &&lhs, t_rhs_rng &&rhs, t_adapter &&adapter)
	{
	
		// We would like to treat runs of indels as equivalent as long as the operation counts match.
		// Otherwise we expected the CIGAR sequences to have been normalised, i.e. consecutive operations are of distinct types.
	
		struct counts
		{
			std::size_t insertions{};
			std::size_t deletions{};
			std::size_t padding{};
		
			void reset() { insertions = 0; deletions = 0; padding = 0; }
			bool operator==(counts const &other) const { return insertions == other.insertions && deletions == other.deletions; }
		};
	
		// Count the indels in the current run, return wheter the iterator can be dereferenced.
		auto count_indels([&adapter](auto &it, auto const end, auto &counts){
			while (it != end)
			{
				auto const op(adapter.operation(*it));
			
				if (t_adapter::insertion_op == op)
				{
					counts.insertions += adapter.count(*it);
					++it;
					continue;
				}
			
				if (t_adapter::deletion_op == op)
				{
					counts.deletions += adapter.count(*it);
					++it;
					continue;
				}
				
				if constexpr (t_should_count_padding)
				{
					if (t_adapter::padding_op == op)
					{
						counts.padding += adapter.count(*it);
						++it;
						continue;
					}
				}
			
				return it != end;
			}
		
			return false;
		});
	
		auto lhs_it(ranges::begin(lhs));
		auto rhs_it(ranges::begin(rhs));
		auto const lhs_end(ranges::end(lhs));
		auto const rhs_end(ranges::end(rhs));
		counts lhs_counts;
		counts rhs_counts;
		while (true)
		{
			lhs_counts.reset();
			rhs_counts.reset();
		
			auto const r1(count_indels(lhs_it, lhs_end, lhs_counts));
			auto const r2(count_indels(rhs_it, rhs_end, rhs_counts));
		
			// If the current run has different counts, the sequences differ.
			if (lhs_counts != rhs_counts) return false;
		
			// Check the current character.
			auto const res(r1 | (r2 << 0x1));
			switch (res)
			{
				case 0x0:
					// At the end of both sequences.
					return true;
			
				case 0x1:
				case 0x2:
					// Only one of the sequences has an extra item, so it does not match.
					return false;
			
				case 0x3:
					// Both sequences have remaining characters
					if (*lhs_it != *rhs_it)
						return false;
					break;
			
				default:
					libbio_fail("Unexpected state");
			}
		
			++lhs_it;
			++rhs_it;
		}
	}
	
	
	template <bool t_should_count_padding = false, typename t_lhs_rng, typename t_rhs_rng>
	bool cigar_eq_libbio(t_lhs_rng &&lhs, t_rhs_rng &&rhs)
	{
		return cigar_eq_tpl <t_should_count_padding>(std::forward <t_lhs_rng>(lhs), std::forward <t_rhs_rng>(rhs), cigar_adapter_libbio{});
	}
	
	
	template <bool t_should_count_padding = false, typename t_lhs_rng, typename t_rhs_rng>
	bool cigar_eq_seqan3(t_lhs_rng &&lhs, t_rhs_rng &&rhs)
	{
		return cigar_eq_tpl <t_should_count_padding>(std::forward <t_lhs_rng>(lhs), std::forward <t_rhs_rng>(rhs), cigar_adapter_seqan3{});
	}
}

#endif
