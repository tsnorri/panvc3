/*
 * Copyright (c) 2023–2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/dispatch.hh>


namespace panvc3::dispatch {
	
	void thread_local_queue::async_(task &&tt)
	{
		{
			std::lock_guard lock(m_mutex);
			m_task_queue.emplace_back(std::move(tt), nullptr);
		}

		m_cv.notify_one();
	}
	
	
	void thread_local_queue::group_async(group &gg, task tt)
	{
		gg.enter();

		{
			std::lock_guard lock(m_mutex);
			m_task_queue.emplace_back(std::move(tt), &gg);
		}

		m_cv.notify_one();
	}
	
	
	void thread_local_queue::clear()
	{
		std::lock_guard lock(m_mutex);
		m_task_queue.clear();
	}
	
	
	void thread_local_queue::stop()
	{
		{
			std::lock_guard lock(m_mutex);
			m_should_continue = false;
		}
		
		m_cv.notify_one();
	}
	
	
	bool thread_local_queue::run()
	{
		std::unique_lock lock(m_mutex);
		while (true)
		{
			// Critical section; we now have the lock.
			if (!m_should_continue)
				return m_task_queue.empty(); // std::unique_lock unlocks automatically.
			
			if (m_task_queue.empty())
			{
				// Task queue is empty and we still hold the lock.
				m_cv.wait(lock); // Unlocks, waits, locks again.
				continue;
			}

			{
				// Get the next item from the queue.
				queue_item item;
				
				{
					using std::swap;
					swap(item, m_task_queue.front());
					m_task_queue.pop_front();
				}
				
				lock.unlock();
				
				// Non-critical section; execute the task.
				item.task_();
				if (item.group_)
					item.group_->exit();
				
				lock.lock();
			}
		}
	}
	
	
	thread_local_queue &main_queue()
	{
		static thread_local_queue main_queue;
		return main_queue;
	}
}
