/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cmath>				// std::floor
#include <libbio/assert.hh>
#include <panvc3/dispatch.hh>

namespace chrono	= std::chrono;


namespace panvc3::dispatch::detail {
	
	struct serial_queue_executor_callable : public parametrised <>::callable
	{
		friend class dispatch::serial_queue;
	
	private:
		serial_queue	*queue{};
		
	private:
		serial_queue_executor_callable(serial_queue &queue_):
			queue(&queue_)
		{
		}
	
	public:
		serial_queue_executor_callable() = delete;
		void execute() override;
		void move_to(std::byte *buffer) override { new(buffer) serial_queue_executor_callable{*queue}; }
		inline void enqueue_transient_async(struct queue &qq) override;
	};
	
	
	void serial_queue_executor_callable::enqueue_transient_async(struct queue &qq)
	{
		// For the sake of completeness; just copy *this.
		qq.async(detail::serial_queue_executor_callable{*queue});
	}
}


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
				for (auto queue : pool.m_queues) // Does not modify m_queues, hence thread-safe.
				{
					if (queue->m_task_queue.try_dequeue(item))
					{
						++executed_queue_items;
						libbio_assert(item.barrier_);
					
						{
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
			// See the note about m_mutex in the beginning of this loop.
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
	
	
	void thread_pool::stop()
	{
		{
			std::lock_guard lock(m_mutex);
			m_should_continue = false;
		}
		
		m_cv.notify_all();
	}
	
	
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
		assert(0 != (~NOTIFY_MASK & res));
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
	
	
	void detail::serial_queue_executor_callable::execute()
	{
		assert(queue);
		
		serial_queue::queue_item item;
		while (queue->fetch_next_task(item))
		{
			item.task_();
			if (item.group_)
				item.group_->exit();
		}
	}
	
	
	void parallel_queue::enqueue(queue_item &&item)
	{
		m_task_queue.enqueue(std::move(item));
		m_thread_pool->start_workers_if_needed();
	}
	
	
	void parallel_queue::async(task tt)
	{
		enqueue(queue_item{
			std::move(tt),
			nullptr,
			current_barrier()
		});
	}
	
	
	void parallel_queue::group_async(group &gg, task tt)
	{
		gg.enter();
		enqueue(queue_item{std::move(tt), &gg, current_barrier()});
	}
	
	
	void parallel_queue::barrier(task tt)
	{
		// Prepare the barrier.
		auto bb(std::make_shared <class barrier>(std::move(tt)));
	
		// Store the barrier and set up the linked list.
		{
			auto old_barrier(m_current_barrier.exchange(bb, std::memory_order_acq_rel)); // Returns std::shared_ptr
			old_barrier->m_next.store(bb, std::memory_order_release); // Safe b.c. m_next is atomic and this is the only place where it is modified.
			// old_barrier (the std::shared_ptr) is now deallocated, and the object it points to may be too.
		}
		
		// Make sure the barrier’s task gets executed at some point by adding an empty task.
		enqueue(queue_item{
			task{},
			nullptr,
			std::move(bb)
		});
	}
	
	
	void serial_queue::enqueue(queue_item &&item)
	{
		bool has_thread{};
		
		{
			std::lock_guard lock(m_mutex);
			m_task_queue.emplace(std::move(item));
			has_thread = m_has_thread;
			m_has_thread = true;
		}
		
		if (!has_thread)
			m_parent_queue->async(detail::serial_queue_executor_callable{*this});
	}
	
	
	void serial_queue::async_(task &&tt)
	{
		enqueue(queue_item{std::move(tt), nullptr});
	}
	
	
	void serial_queue::group_async(group &gg, task tt)
	{
		gg.enter();
		enqueue(queue_item{std::move(tt), &gg});
	}
	
	
	bool serial_queue::fetch_next_task(queue_item &item)
	{
		std::lock_guard lock(m_mutex);
		
		if (m_task_queue.empty())
		{
			m_has_thread = false;
			return false;
		}
		
		item = std::move(m_task_queue.front());
		m_task_queue.pop();
		return true;
	}
}
