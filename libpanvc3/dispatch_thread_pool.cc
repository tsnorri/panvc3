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


	void block_signals()
	{
		// Block all signals.
		sigset_t mask{};
		if (-1 == sigfillset(&mask))
			throw std::runtime_error(::strerror(errno));
		if (-1 == ::pthread_sigmask(SIG_SETMASK, &mask, nullptr))
			throw std::runtime_error(::strerror(errno));
	}
	
	
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
		
		void remove_from_pool(std::int64_t const executed_tasks);
		void begin_idle(std::int64_t const executed_tasks);
	};
	
	
	void worker_thread_runner::remove_from_pool(std::int64_t const executed_tasks)
	{
		auto &pool(*m_thread_pool);
		pool.m_waiting_tasks -= executed_tasks;
		pool.remove_worker();
	}
	
	
	void worker_thread_runner::begin_idle(std::int64_t const executed_tasks)
	{
		auto &pool(*m_thread_pool);
		pool.m_waiting_tasks -= executed_tasks;
		++pool.m_idle_workers;
	}
	
	
	void worker_thread_runner::run()
	{
		block_signals();

		libbio_assert(m_thread_pool);
		auto &pool(*m_thread_pool);
		
		auto last_wake_up_time{clock_type::now()};
		parallel_queue::queue_item queue_item{}; // Initialise just once to save a bit of time.
		
		{
			std::unique_lock lock(pool.m_mutex, std::defer_lock_t{});
			while (true)
			{
				std::int64_t executed_tasks{}; // Total over the iterations of the loop below (but not the enclosing loop).
			
				{
					// Critical section 1.
					// We need the queues to persist while tasks are being executed.
					// Acquiring m_queue_mutex is sufficient b.c. the queue will not get deallocated before
					// it has been removed from the thread pool.
					std::shared_lock const queue_lock{pool.m_queue_mutex};
					while (true)
					{
						auto const prev_executed_tasks(executed_tasks);
						for (auto queue : pool.m_queues)
						{
							if (queue->m_task_queue.try_dequeue(queue_item))
							{
								++executed_tasks;
						
#if PANVC3_ENABLE_DISPATCH_BARRIER
								{
									libbio_assert(queue_item.barrier_);
									auto &bb(*queue_item.barrier_);
									barrier::status_underlying_type state{barrier::NOT_EXECUTED};
									if (bb.m_state.compare_exchange_strong(state, barrier::EXECUTING, std::memory_order_acq_rel, std::memory_order_acquire))
									{
										// Wait for the previous tasks and the previous barrier to complete.
										bb.m_previous_has_finished.wait(false, std::memory_order_acquire);
								
										bb.m_task();
										bb.m_task = task{}; // Deallocate memory.
								
										bool should_continue{};
										
										{
											std::lock_guard const lock_{lock};
											should_continue = pool.m_should_continue;
											if (!should_continue)
												remove_from_pool(executed_tasks);
										}
								
										if (should_continue)
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
												{
													std::lock_guard const lock_{lock};
													remove_from_pool(executed_tasks);
													return;
												}
										
												break;
											}
									
											case barrier::DONE:
												break;
									
											// Stop if the barrier’s task called m_pool.stop().
											case barrier::DO_STOP:
											{
												std::lock_guard const lock_{lock};
												remove_from_pool(executed_tasks);
												return;
											}
										
											case barrier::NOT_EXECUTED:
												// Unexpected.
												std::abort();
										}
									}
								}
#endif
						
								queue_item.task_();
								if (queue_item.group_)
									queue_item.group_->exit(); // Important to do only after executing the task, since it can add new tasks to the group.
							}
						} // Queue loop
					
						if (executed_tasks == prev_executed_tasks)
							break;
					} // Inner while (true)
				} // Critical section 1
			
				{
					// Check the last wake-up time.
					auto const now{clock_type::now()};
					auto const diff(now - last_wake_up_time);
					if (0 == executed_tasks && m_max_idle_time <= diff)
					{
						// Critical section 2.
						lock.lock();
						remove_from_pool(executed_tasks); // zero but does not matter.
						return;
					}
				
					last_wake_up_time = now;
				}
			
				// Critical section 2.
				{
					lock.lock();
					begin_idle(executed_tasks);
				
					// Handle spurious wake-ups by repeatedly calling wait_for().
					while (true)
					{
						// m_mutex is locked when wait_until() returns.
						switch (pool.m_cv.wait_for(lock, m_max_idle_time))
						{
							case std::cv_status::no_timeout:
								break;
					
							case std::cv_status::timeout:
								// Still marked idle.
								pool.remove_idle_worker();
								return;
						}
						
						if (!pool.m_should_continue)
						{
							pool.remove_idle_worker();
							return;
						}
						
						if (pool.m_notified_workers)
						{
							--pool.m_notified_workers;
							break;
						}
					}
					
					lock.unlock();
				}
			} // Outer while (true)
		}
	}
	
	
	void thread_pool::add_queue(parallel_queue &queue)
	{
		std::lock_guard const lock{m_queue_mutex}; // Exclusive
		m_queues.push_back(&queue);
	}
	
	
	void thread_pool::remove_queue(parallel_queue const &queue)
	{
		std::lock_guard const lock{m_queue_mutex}; // Exclusive
		auto const it(std::find(m_queues.begin(), m_queues.end(), &queue));
		if (it == m_queues.end())
			return;
		
		m_queues.erase(it);
	}
	
	
	void thread_pool::start_worker_()
	{
		++m_current_workers;
		std::thread thread(worker_thread_runner{*this, m_max_idle_time});
		thread.detach();
	}
	
	
	void thread_pool::remove_worker()
	{
		libbio_assert_lt(0, m_current_workers);
		if (0 == --m_current_workers)
			m_stop_cv.notify_one();
	}
	
	
	void thread_pool::remove_idle_worker()
	{
		// Worker still marked idle.
		libbio_assert_lt(0, m_current_workers);
		if (0 == --m_current_workers)
			m_stop_cv.notify_one();
	}
	
	
	void thread_pool::start_worker()
	{
		std::lock_guard const lock{m_mutex};
		start_worker_();
	}
	
	
	void thread_pool::notify()
	{
		{
			std::lock_guard const lock{m_mutex};
			++m_waiting_tasks;
			if (m_idle_workers)
			{
				--m_idle_workers;
				++m_notified_workers;
				goto do_notify;
			}
			
			if (m_max_workers <= m_current_workers)
				return;
			
			// Can start a new thread.
			start_worker_();
		}
		
	do_notify:
		m_cv.notify_one();
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
		std::unique_lock lock{m_mutex};
		if (0 < m_current_workers)
			m_stop_cv.wait(lock, [this]{ return 0 == m_current_workers; });
	}
}
