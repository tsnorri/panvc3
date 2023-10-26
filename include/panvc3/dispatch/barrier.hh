/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_BARRIER_HH
#define PANVC3_DISPATCH_BARRIER_HH

#if PANVC3_ENABLE_DISPATCH_BARRIER

#include <atomic>
#include <cstdint>
#include <limits>	// std::numeric_limits
#include <memory>	// std::shared_ptr
#include <panvc3/dispatch/fwd.hh>
#include <panvc3/dispatch/task_decl.hh>


namespace panvc3::dispatch {
	
	// When a barrier B is added to the queue, it is not added to the lock-free queue but instead:
	//	– B’ = m_current_barrier.exchange(B, std::memory_order_release) is used.
	//	– B’.m_next is set to B
	// When a thread completes a barrier B, it sets m_next’s m_previous_has_finished to true.
	// When a regular task is added to a queue, m_current_barrier.load(std::memory_order_acquire) is called and the result is moved to the task.
	// When a task is executed:
	//	– m_state.compare_exchange_strong(NOT_EXECUTED, EXECUTING, std::memory_order_acq_rel, std::memory_order_acquire) is used.
	//		– If successful:
	//			– The thread waits until the barrier’s m_previous_has_finished is true,
	//			– The barrier is executed,
	//			– m_state.store(DONE, std::memory_order_release) is called,
	//			– m_state.notify_all() is called,
	//			– The original task is executed.
	//		– If S == EXECUTING:
	//			– m_state.wait(EXECUTING, std::memory_order::acquire) is called (waits iff. m_state == EXECUTING),
	//	– The original task is executed.
	// When a barrier is deallocated (neither parallel_queue’s m_current_barrier nor any queue_item points to it), m_next’s m_previous_has_finished is set to true and the waiting thread is notified.
	class barrier
	{
		friend class parallel_queue;
		friend class worker_thread_runner;
		
	public:
		typedef std::shared_ptr <barrier>			shared_barrier_ptr;
		typedef std::atomic <shared_barrier_ptr>	atomic_shared_barrier_ptr;
		
	private:
		//typedef std::atomic_unsigned_lock_free		status_type;
		typedef std::atomic <std::uint32_t>			status_type;
		typedef status_type::value_type				status_underlying_type;
		
		enum state
		{
			NOT_EXECUTED	= 0,
			EXECUTING		= 1,
			DONE			= 2,
			DO_STOP			= 3,
			
			STATE_MIN		= NOT_EXECUTED,
			STATE_MAX		= DO_STOP
		};
		
		static_assert(std::numeric_limits <status_underlying_type>::min() <= STATE_MIN);
		static_assert(STATE_MAX <= std::numeric_limits <status_underlying_type>::max());
		
		task						m_task;
		atomic_shared_barrier_ptr	m_next{};
		status_type					m_state{DONE}; // FIXME: Consider cache line size?
		status_type					m_previous_has_finished{true};
		
	public:
		barrier() = default;
		
		barrier(task tt):
			m_task(std::move(tt)),
			m_state(NOT_EXECUTED),
			m_previous_has_finished(false)
		{
		}
		
		inline ~barrier();
	};
	
	typedef barrier::shared_barrier_ptr		shared_barrier_ptr;
	
	
	barrier::~barrier()
	{
		auto next(m_next.load(std::memory_order_acquire)); // std::shared_ptr
		if (next)
		{
			// Safe b.c. next is a shared_ptr.
			auto &previous_has_finished(next->m_previous_has_finished);
			previous_has_finished.store(true, std::memory_order_release);
			previous_has_finished.notify_one();
		}
	}
}

#endif

#endif
