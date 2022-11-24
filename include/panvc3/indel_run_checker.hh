/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_INDEL_RUN_CHECKER_HH
#define PANVC3_INDEL_RUN_CHECKER_HH

#include <panvc3/cigar.hh>
#include <range/v3/view/slice.hpp>
#include <vector>


namespace panvc3 {
	
	class indel_run_checker
	{
	public:
		typedef std::vector <char>	sequence_vector;
		
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
		};
		
	protected:
		enum
		{
			RUN_HAS_INSERTIONS	= 0x1,
			RUN_HAS_DELETIONS	= 0x2,
			RUN_HAS_BOTH		= RUN_HAS_INSERTIONS | RUN_HAS_DELETIONS,
			RUN_HAS_NEITHER		= 0x0
		};
		
	protected:
		panvc3::cigar_vector::const_iterator	m_cigar_it{};
		panvc3::cigar_vector::const_iterator	m_cigar_end{};
		panvc3::cigar_vector::const_iterator	m_cigar_realigned_range_begin{};
		range									m_ref_range{};
		range									m_query_range{};
		std::size_t								m_ref_pos{};
		std::size_t								m_query_pos{};
		std::uint8_t							m_run_type{};
		
	public:
		void reset(panvc3::cigar_vector const &cigar_vector, std::size_t const ref_pos);
		bool find_next_range_for_realigning();
		range reference_range() const { return m_ref_range; }
		range query_range() const { return m_query_range; }
		panvc3::cigar_vector::const_iterator cigar_current_pos() const { return m_cigar_it; }
		panvc3::cigar_vector::const_iterator cigar_realigned_range_begin() const { return m_cigar_realigned_range_begin; }
		
	protected:
		void update_ranges(std::size_t const ref_pos, std::size_t const query_pos);
	};
}

#endif
