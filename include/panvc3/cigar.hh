/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_CIGAR_HH
#define PANVC3_CIGAR_HH

#include <libbio/assert.hh>
#include <range/v3/range/access.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <vector>


namespace panvc3 {
	
	typedef std::vector <seqan3::cigar>	cigar_vector;
	
	
	class cigar_buffer
	{
	protected:
		cigar_vector	m_operations;
		seqan3::cigar	m_current_op;
		bool			m_is_first{true};
		
	public:
		void push_back(seqan3::cigar::operation const op, std::uint32_t const count = 1);
		void finish();
		void clear();
		
		cigar_vector &operations() { return m_operations; }
		cigar_vector const &operations() const { return m_operations; }
	};
	
	
	template <bool t_should_count_padding = false, typename t_lhs_rng, typename t_rhs_rng>
	bool cigar_eq(t_lhs_rng &&lhs, t_rhs_rng &&rhs)
	{
		using seqan3::get;
		using seqan3::operator""_cigar_operation;
	
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
	
		auto const cigar_count([](auto const cigar_item){
			return get <0>(cigar_item);
		});
	
		auto const cigar_op([](auto const cigar_item){
			return get <1>(cigar_item);
		});
	
		// Count the indels in the current run, return wheter the iterator can be dereferenced.
		auto count_indels([&cigar_op, &cigar_count](auto &it, auto const end, auto &counts){
			while (it != end)
			{
				auto const op(cigar_op(*it));
			
				if ('I'_cigar_operation == op)
				{
					counts.insertions += cigar_count(*it);
					++it;
					continue;
				}
			
				if ('D'_cigar_operation == op)
				{
					counts.deletions += cigar_count(*it);
					++it;
					continue;
				}
				
				if constexpr (t_should_count_padding)
				{
					if ('P'_cigar_operation == op)
					{
						counts.padding += cigar_count(*it);
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
}

#endif
