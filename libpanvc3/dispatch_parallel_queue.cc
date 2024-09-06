/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/dispatch.hh>


namespace panvc3::dispatch {
	
	parallel_queue &parallel_queue::shared_queue()
	{
		static parallel_queue queue;
		return queue;
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
			nullptr
#if PANVC3_ENABLE_DISPATCH_BARRIER
			, current_barrier()
#endif
		});
	}
	
	
	void parallel_queue::group_async(group &gg, task tt)
	{
		gg.enter();
		enqueue(queue_item{
			std::move(tt),
			&gg
#if PANVC3_ENABLE_DISPATCH_BARRIER
			, current_barrier()
#endif
		});
	}
	
	
#if PANVC3_ENABLE_DISPATCH_BARRIER
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
#endif
}
