/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/cigar.hh>

namespace lb	= libbio;


namespace panvc3 {

	void collapse_cigar_operations(/* inout */ cigar_vector &ops)
	{
		// Essentially a conditional copy.
		auto const begin(ops.begin());
		auto it(begin);
		auto it_(it + 1);
		auto const end(ops.end());

		// If the operation of the item at the reading position (*it_) matches that of
		// the current item at the writing position (*it), add to the count. Otherwise,
		// copy to the next place (i.e. *(it + 1) = *it_).
		while (it_ != end)
		{
			using seqan3::get;

			auto &curr_item(*it);
			auto const curr_op(get <1>(curr_item));
			auto const next_item(*it_);
			auto const next_op(get <1>(next_item));

			if (curr_op == next_op)
			{
				auto const curr_count(get <0>(curr_item));
				auto const next_count(get <0>(next_item));
				curr_item = curr_count + next_count;
				++it_;
			}
			else
			{
				auto &dst_item(*(it + 1));
				dst_item = next_item;
				++it;
				++it_;
			}
		}

		// Remove the excess items.
		ops.resize(std::distance(begin, it + 1));
	}

	
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
