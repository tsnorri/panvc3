/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/assert.hh>
#include <panvc3/dispatch.hh>


namespace panvc3::dispatch::detail {
	
	struct serial_queue_executor_callable : public callable <>
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
}


namespace panvc3::dispatch {
	
	void detail::serial_queue_executor_callable::enqueue_transient_async(struct queue &qq)
	{
		// For the sake of completeness; just copy *this.
		qq.async(detail::serial_queue_executor_callable{*queue});
	}
	
	
	void detail::serial_queue_executor_callable::execute()
	{
		libbio_assert(queue);
		
		serial_queue::queue_item item;
		while (queue->fetch_next_task(item))
		{
			item.task_();
			if (item.group_)
				item.group_->exit();
		}
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
