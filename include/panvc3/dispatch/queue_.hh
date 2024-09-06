/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_QUEUE__HH
#define PANVC3_DISPATCH_QUEUE__HH

#ifndef DISPATCH_USE_CONCURRENT_TASK_QUEUE
#	define DISPATCH_USE_CONCURRENT_TASK_QUEUE 0
#endif

#include <list>
#include <mutex>
#include <queue>

#if DISPATCH_USE_CONCURRENT_TASK_QUEUE
#	include <concurrentqueue/concurrentqueue.h>
#endif


namespace panvc3::dispatch::detail {
	
	template <typename t_item>
	class blocking_queue
	{
	private:
		typedef std::list <t_item>	queue_type;
		
	private:
		queue_type	m_queue;
		std::mutex	m_mutex;
		
	public:
		inline void clear();
		inline void enqueue(t_item &&item);
		inline bool try_dequeue(t_item &item);
	};
	
	
	template <typename t_item>
	void blocking_queue <t_item>::clear()
	{
		std::lock_guard lock(m_mutex);
		m_queue.clear();
	}
	
	
	template <typename t_item>
	void blocking_queue <t_item>::enqueue(t_item &&item)
	{
		std::lock_guard lock(m_mutex);
		m_queue.emplace_back(std::move(item));
	}
	
	
	template <typename t_item>
	bool blocking_queue <t_item>::try_dequeue(t_item &item)
	{
		std::lock_guard lock(m_mutex);
		
		if (m_queue.empty())
			return false;
		
		using std::swap;
		swap(m_queue.front(), item);
		m_queue.pop_front();
		return true;
	}
	
	
#if DISPATCH_USE_CONCURRENT_TASK_QUEUE
#	error clear() not implemented.
	template <typename t_item>
	using queue_t = moodycamel::ConcurrentQueue <t_item>;
#else
	template <typename t_item>
	using queue_t = blocking_queue <t_item>;
#endif
}

#endif
