/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cmath>				// std::floor
#include <libbio/assert.hh>
#include <panvc3/dispatch.hh>
#include <thread>

namespace chrono	= std::chrono;


namespace panvc3::dispatch {
	
	thread_pool::thread_count_type const thread_pool::default_max_worker_threads = thread_pool::thread_count_type(
		std::floor(1.5 * std::thread::hardware_concurrency())
	);
	
	
	class worker_thread_runner
	{
	private:
		typedef	chrono::steady_clock	clock_type;
		typedef clock_type::time_point	time_point_type;
		typedef clock_type::duration	duration_type;
		
	private:
		thread_pool		*m_thread_pool{};
		duration_type	m_max_idle_time{};
		
	public:
		worker_thread_runner(thread_pool &pool, duration_type const max_idle_time):
			m_thread_pool(&pool),
			m_max_idle_time(max_idle_time)
		{
		}
		
		void run();
		void operator()() { run(); }
	};
	
	
	void worker_thread_runner::run()
	{
		libbio_assert(m_thread_pool);
		auto &pool(*m_thread_pool);
		
		std::unique_lock lock(pool.m_mutex, std::defer_lock_t{});
		auto last_wake_up_time{clock_type::now()};
		while (true)
		{
			parallel_queue::queue_item item{};
			std::uint64_t executed_queue_items{};
			while (true)
			{
				auto const prev_executed_queue_items(executed_queue_items);
				for (auto queue : pool.m_queues) // Does not modify m_queues, hence thread-safe (as long as the queues have been pre-added).
				{
					if (queue->m_task_queue.try_dequeue(item))
					{
						++executed_queue_items;
						
#if PANVC3_ENABLE_DISPATCH_BARRIER
						{
							libbio_assert(item.barrier_);
							auto &bb(*item.barrier_);
							barrier::status_underlying_type state{barrier::NOT_EXECUTED};
							if (bb.m_state.compare_exchange_strong(state, barrier::EXECUTING, std::memory_order_acq_rel, std::memory_order_acquire))
							{
								// Wait for the previous tasks and the previous barrier to complete.
								bb.m_previous_has_finished.wait(false, std::memory_order_acquire);
								
								bb.m_task();
								bb.m_task = task{}; // Deallocate memory.
								
								// Not in critical section but no atomicity required b.c. m_task (that potentially calls m_pool.stop())
								// is executed in the same thread. (Otherwise it cannot be expected that the pending tasks would not continue.)
								if (pool.m_should_continue)
								{
									bb.m_state.store(barrier::DONE, std::memory_order_release);
									bb.m_state.notify_all();
								}
								else
								{
									bb.m_state.store(barrier::DO_STOP, std::memory_order_release);
									bb.m_state.notify_all();
									return;
								}
							}
							else
							{
								// The barrier task is either currently being executed or has already been finished.
								switch (state)
								{
									case barrier::EXECUTING:
									{
										bb.m_state.wait(barrier::EXECUTING, std::memory_order::acquire);
										// The acquire operation above should make the modification visible here.
										if (barrier::DO_STOP == bb.m_state.load(std::memory_order_relaxed))
											return;
										
										break;
									}
									
									case barrier::DONE:
										break;
									
									// Stop if the barrier’s task called m_pool.stop().
									case barrier::DO_STOP:
										return;
										
									case barrier::NOT_EXECUTED:
										// Unexpected.
										std::abort();
								}
							}
						}
#endif
						
						item.task_();
						if (item.group_)
							item.group_->exit(); // Important to do only after executing the task, since it can add new tasks to the group.
						break;
					}
				}
				
				if (executed_queue_items == prev_executed_queue_items)
					break;
			}
			
			{
				// Check the last wake-up time.
				auto const now{clock_type::now()};
				auto const diff(now - last_wake_up_time);
				if (0 == executed_queue_items && m_max_idle_time <= diff)
					break;
				
				last_wake_up_time = now;
			}
			
			// Critical section.
			lock.lock();
			// See the note about m_mutex after the switch statement.
			pool.m_sleeping_workers.fetch_add(1, std::memory_order_relaxed);
			
			switch (pool.m_cv.wait_until(lock, last_wake_up_time + m_max_idle_time))
			{
				case std::cv_status::no_timeout:
					break;
					
				case std::cv_status::timeout:
					goto end_worker_loop;
			}
			
			// Critical section.
			// wait_until() locks m_mutex, which is an acquire operation and hence the modification to pool.m_should_continue
			// done in the critical section in stop() is visible here. Similarly the unlock is a release operation and
			// makes the fetch_sub operation below visible elsewhere.
			if (!pool.m_should_continue)
				goto end_worker_loop;
			
			lock.unlock();
		}
		
	end_worker_loop:
		pool.m_workers.fetch_sub(1, std::memory_order_release);
		pool.m_workers.notify_one();
	}
	
	
	void thread_pool::start_worker()
	{
		m_workers.fetch_add(1, std::memory_order_relaxed);
		std::thread thread(worker_thread_runner{*this, m_max_idle_time});
		thread.detach();
	}
	
	
	void thread_pool::start_workers_if_needed()
	{
		auto const prev_sleeping_threads(m_sleeping_workers.fetch_sub(1, std::memory_order_acquire));
		if (0 < prev_sleeping_threads)
		{
			m_cv.notify_one();
			return;
		}
		
		m_sleeping_workers.fetch_add(1, std::memory_order_release);
		auto const workers(m_workers.load(std::memory_order_relaxed));
		if (workers < m_max_workers)
			start_worker();
	}
	
	
	void thread_pool::stop(bool should_wait)
	{
		{
			std::lock_guard lock(m_mutex);
			m_should_continue = false;
		}
		
		m_cv.notify_all();
		
		if (should_wait)
			wait();
	}
	
	
	void thread_pool::wait()
	{
		// FIXME: Not particularly efficient.
		while (true)
		{
			auto const count(m_workers.load(std::memory_order_acquire));
			if (0 == count)
				break;
			m_workers.wait(count, std::memory_order_acquire); // Wait until the value is no longer count.
		}
	}
}
