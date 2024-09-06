/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/assert.hh>
#include <panvc3/dispatch.hh>


namespace panvc3::dispatch {
	
	void group::wait()
	{
		exit();
		
		std::unique_lock lock(m_mutex);
		m_cv.wait(lock, [this]{ return m_should_stop_waiting; });
		m_should_stop_waiting = false;
		m_count.fetch_add(1, std::memory_order_relaxed); // Restore the group’s initial state.
	}
	
	
	void group::notify(struct queue &queue, task tt)
	{
		m_queue = &queue;
		m_task = std::move(tt);
		// Relaxed should be enough b.c. exit() uses std::memory_order_acq_rel.
		m_count.fetch_or(NOTIFY_MASK, std::memory_order_relaxed);
		
		exit();
	}
	
	
	void group::exit()
	{
		auto const res(m_count.fetch_sub(1, std::memory_order_acq_rel));
		libbio_assert_neq(0, ~NOTIFY_MASK & res);
		if ((NOTIFY_MASK | 1) == res)
		{
			m_queue->async(std::move(m_task));
			m_queue = nullptr;
			m_task = task{};
			m_count.fetch_add(1, std::memory_order_relaxed); // Restore the group’s initial state.
		}
		else if (1 == res)
		{
			{
				std::lock_guard lock(m_mutex);
				m_should_stop_waiting = true;
			}
			
			m_cv.notify_all();
		}
	}
}
