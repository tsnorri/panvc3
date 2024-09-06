/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_QUEUE_HH
#define PANVC3_DISPATCH_QUEUE_HH

#include <atomic>
#include <panvc3/dispatch/fwd.hh>
#include <panvc3/dispatch/barrier.hh>
#include <panvc3/dispatch/group.hh>
#include <panvc3/dispatch/queue_.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <panvc3/dispatch/thread_pool.hh>


namespace panvc3::dispatch {
	
	struct queue
	{
		virtual ~queue() {}
		virtual void clear() = 0;
		virtual void async(task tt) = 0;
		virtual void group_async(group &gg, task tt) = 0;
		
#if PANVC3_ENABLE_DISPATCH_BARRIER
		virtual void barrier(task tt) = 0;
#endif
	};
	
	
	class parallel_queue final : public queue
	{
		friend class worker_thread_runner;
		
	private:
		struct queue_item
		{
			task				task_{};
			group				*group_{};
#if PANVC3_ENABLE_DISPATCH_BARRIER
			shared_barrier_ptr	barrier_{};
#endif
		};
		
	private:
		typedef detail::queue_t <queue_item>	queue_type;
		
	private:
		thread_pool							*m_thread_pool{};
		queue_type							m_task_queue;
#if PANVC3_ENABLE_DISPATCH_BARRIER
		std::atomic <shared_barrier_ptr>	m_current_barrier{std::make_shared <class barrier>()};	// To get atomic store() etc.
#endif
		
	public:
		parallel_queue(): // Not thread-safe.
			m_thread_pool(&thread_pool::shared_pool())
		{
			m_thread_pool->add_queue(*this);
		}
		
		parallel_queue(thread_pool &pool): // Not thread-safe.
			m_thread_pool(&pool)
		{
			m_thread_pool->add_queue(*this);
		}
		
		~parallel_queue();
		
		static parallel_queue &shared_queue();
		
		parallel_queue(parallel_queue const &) = delete;
		parallel_queue(parallel_queue &&) = delete;
		parallel_queue &operator=(parallel_queue const &) = delete;
		parallel_queue &operator=(parallel_queue &&) = delete;
		
		void clear() override;
			
		void async(task tt) override;
		void group_async(group &gg, task tt) override;
		
#if PANVC3_ENABLE_DISPATCH_BARRIER
		void barrier(task tt) override;
#endif
		
	private:
		inline void enqueue(queue_item &&item);
		
#if PANVC3_ENABLE_DISPATCH_BARRIER
		shared_barrier_ptr current_barrier() const { return m_current_barrier.load(std::memory_order_acquire); }
#endif
	};
	
	
	class serial_queue final : public queue
	{
		friend struct detail::serial_queue_executor_callable;
		
	private:
		struct queue_item
		{
			task	task_{};
			group	*group_{};
		};
		
	private:
		typedef std::deque <queue_item> queue_type;
		
	private:
		parallel_queue			*m_parent_queue{};
		queue_type				m_task_queue;
		std::condition_variable	m_cv;				// For waiting for the serial_queue_executor_callable to finish.
		std::mutex				m_mutex;
		bool					m_has_thread{};
		
	public:
		serial_queue(parallel_queue &parent_queue):
			m_parent_queue(&parent_queue)
		{
		}
		
		~serial_queue();
		
		void clear() override;

		void async(task tt) override { async_(std::move(tt)); }
		void group_async(group &gg, task tt) override;
		
#if PANVC3_ENABLE_DISPATCH_BARRIER
		void barrier(task tt) override { async_(std::move(tt)); }
#endif
	private:
		void async_(task &&tt);
		inline void enqueue(queue_item &&item);
		bool fetch_next_task(queue_item &item);
	};
	
	
	// Implementation of queue’s interface for running tasks on a given thread.
	class thread_local_queue final : public queue
	{
	private:
		struct queue_item
		{
			task	task_{};
			group	*group_{};
		};
		
	private:
		typedef std::deque <queue_item> queue_type;
		
	private:
		queue_type				m_task_queue;
		std::condition_variable	m_cv;
		std::mutex				m_mutex;
		bool					m_should_continue{true};
		
	private:
		void async_(task &&tt);
		
	public:
		void clear() override;
		
		void async(task tt) override { async_(std::move(tt)); }
		void group_async(group &gg, task tt) override;
		
#if PANVC3_ENABLE_DISPATCH_BARRIER
		void barrier(task tt) override { async_(std::move(tt)); }
#endif
		
		bool run();
		void stop();
	};
	
	
	thread_local_queue &main_queue();
}

#endif
