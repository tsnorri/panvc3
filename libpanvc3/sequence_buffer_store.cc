/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/assert.hh>
#include <panvc3/sequence_buffer_store.hh>

namespace lb	= libbio;


namespace panvc3 {

	// Call from thread 1.
	auto sequence_buffer_store::acquire_buffer(std::size_t const idx) -> sequence_buffer &
	{
		libbio_assert_lt(idx, m_items.size());
		auto &retval(m_items[idx]);
		// Synchronise with release_buffer() since we would like everything done in thread 2
		// before that to become visible here.
		auto const prev_count(retval.m_use_count.fetch_add(1, std::memory_order_acq_rel));

		// Find a spare buffer if possible. If 0 < prev_count, the buffer is still in use
		// in thread 2.
		if (0 == prev_count)
		{
			retval.m_buffer.clear();
			if (0 == retval.m_buffer.capacity() && !m_spare_buffers.empty())
			{
				retval.m_buffer = std::move(m_spare_buffers.back());
				m_spare_buffers.pop_back();
			}
		}

		return retval;
	}


	// Call from thread 1 after calling acquire_buffer() on some items.
	void sequence_buffer_store::recycle_buffers()
	{
		for (auto &item : m_items)
		{
			// Synchronise with release_buffer().
			auto const use_count(item.m_use_count.load(std::memory_order_acquire));
			if (0 == use_count)
			{
				item.m_buffer.clear();
				if (m_spare_buffers.size() < SPARE_BUFFER_COUNT)
					m_spare_buffers.emplace_back(std::move(item.m_buffer));
				else
					item.m_buffer = buffer_type(); // Frees memory.
			}
		}
	}


	// Call from any thread.
	auto sequence_buffer_store::buffer(std::size_t const idx) const -> sequence_buffer const &
	{
		auto const &retval(m_items[idx]);

		// Synchronise with acquire_buffer() and release_buffer().
		libbio_assert_lt(0, retval.m_use_count.load(std::memory_order_acquire));

		return retval;
	}


	// Call from thread 2.
	void sequence_buffer_store::release_buffer(std::size_t const idx)
	{
		libbio_assert_lt(idx, m_items.size());
		auto const prev_count(m_items[idx].m_use_count.fetch_sub(1, std::memory_order_release));
		libbio_assert_lt(0, prev_count);
	}
}
