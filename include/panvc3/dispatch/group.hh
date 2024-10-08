/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_GROUP_HH
#define PANVC3_DISPATCH_GROUP_HH

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <panvc3/dispatch/fwd.hh>
#include <panvc3/dispatch/task_decl.hh>


namespace panvc3::dispatch {
	
	class group
	{
		friend class parallel_queue;
		friend class serial_queue;
		friend class thread_local_queue;
		friend class worker_thread_runner;
		friend struct detail::serial_queue_executor_callable;
		
	private:
		typedef std::uint32_t						counter_type;
		constexpr static inline counter_type const	NOTIFY_MASK{0x80000000};
		
	private:
		task						m_task;		// For notify.
		queue						*m_queue{};	// For notify.
		std::atomic <counter_type>	m_count{1};
		std::condition_variable		m_cv{};
		std::mutex					m_mutex{};
		bool						m_should_stop_waiting{};
		
	public:
		void notify(queue &qq, task tt);
		void wait();
		
		void enter() { m_count.fetch_add(1, std::memory_order_relaxed); } // Relaxed ordering in increment should be enough; see https://en.cppreference.com/w/cpp/atomic/memory_order#Relaxed_ordering
		void exit();
	};
}

#endif
