/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/cigar.hh>

namespace lb	= libbio;


namespace panvc3 {
	
	void cigar_buffer::clear()
	{
		m_is_first = true;
		m_operations.clear();
	}
	
	
	void cigar_buffer::push_back(seqan3::cigar::operation const op, cigar_count_type count)
	{
		using seqan3::get;
		
		if (0 == count)
			return;
		
		if (m_is_first)
		{
			m_is_first = false;
			m_current_op = op;		// Assigns to one component.
			m_current_op = count;
		}
		else if (get <1>(m_current_op) == op)
		{
			// The alphabet instance cannot be modified directly.
			count += get <0>(m_current_op);		// Converts magically to integral?
			m_current_op = count;				// Assigns to one component.
		}
		else
		{
			if (get <0>(m_current_op))
				m_operations.emplace_back(m_current_op);
			m_current_op = op;
			m_current_op = count;
		}
	}
	
	
	void cigar_buffer::finish()
	{
		using seqan3::get;
		if (!m_is_first && get <0>(m_current_op))
			m_operations.emplace_back(m_current_op);
	}
}
