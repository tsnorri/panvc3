/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_THREAD_POOL_HH
#define PANVC3_DISPATCH_THREAD_POOL_HH

#include <atomic>
#include <chrono>					// std::chrono::steady_clock etc.
#include <cstdint>
#include <mutex>
#include <panvc3/dispatch/fwd.hh>
#include <shared_mutex>
#include <vector>


namespace panvc3::dispatch {

	void block_signals();

	
	class thread_pool
	{
		friend class worker_thread_runner;
		
	public:
		typedef std::uint32_t thread_count_type;
		constexpr static inline auto const default_max_idle_time{std::chrono::seconds(15)};
		static thread_count_type const default_max_worker_threads;
		
	private:
		typedef std::chrono::steady_clock	clock_type;
		typedef clock_type::duration		duration_type;
		
	private:
		std::vector <parallel_queue *>	m_queues;									// Non-owning.
		std::int64_t					m_waiting_tasks{};
		duration_type					m_max_idle_time{default_max_idle_time};
		thread_count_type				m_max_workers{default_max_worker_threads};
		thread_count_type				m_current_workers{};
		thread_count_type				m_idle_workers{};
		thread_count_type				m_notified_workers{};						// For detecting spurious wake-ups.
		std::condition_variable			m_cv{};										// For pausing the workers.
		std::condition_variable			m_stop_cv{};								// For stopping the thread pool.
		std::shared_mutex				m_queue_mutex{};							// Protects m_queues
		std::mutex						m_mutex{};									// Protects m_waiting_tasks, m_current_workers, m_idle_workers, m_should_continue.
		bool							m_should_continue{true};
		
	private:
		void start_worker_();
		void remove_worker();
		void remove_idle_worker();
		
	public:
		static inline thread_pool &shared_pool();
		void add_queue(parallel_queue &queue);			// Thread-safe.
		void remove_queue(parallel_queue const &queue);	// Thread-safe.
		void stop(bool should_wait = true);				// Thread-safe.

		void set_max_workers(thread_count_type const count) { m_max_workers = count; }
		
		void notify();									// Task was added to an observed queue. Thread-safe.
		void wait();
		void start_worker();							// Thread-safe.
		~thread_pool() { stop(); } // parallel_queue expects its thread_pool to persist until the queue has been deallocated.
	};
	
	
	thread_pool &thread_pool::shared_pool()
	{
		static thread_pool pool;
		return pool;
	}
}

#endif
