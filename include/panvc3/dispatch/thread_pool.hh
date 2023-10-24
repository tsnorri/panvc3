/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_THREAD_POOL_HH
#define PANVC3_DISPATCH_THREAD_POOL_HH

#include <atomic>
#include <chrono>					// std::chrono::steady_clock etc.
#include <cstdint>
#include <mutex>
#include <panvc3/dispatch/fwd.hh>
#include <vector>


namespace panvc3::dispatch {
	
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
		duration_type					m_max_idle_time{default_max_idle_time};
		thread_count_type				m_max_workers{default_max_worker_threads};
		std::condition_variable			m_cv{};										// For stopping the workers.
		std::mutex						m_mutex{};
		std::atomic <thread_count_type>	m_workers{};
		std::atomic <thread_count_type>	m_sleeping_workers{};
		bool							m_should_continue{true};
		
	private:
		void start_worker();
		
	public:
		static inline thread_pool &shared_pool();
		void add_queue(parallel_queue &queue) { m_queues.push_back(&queue); }	// Not thread-safe.
		void start_workers_if_needed();
		void stop();
		~thread_pool() { stop(); }
	};
	
	
	thread_pool &thread_pool::shared_pool()
	{
		static thread_pool pool;
		return pool;
	}
}

#endif
