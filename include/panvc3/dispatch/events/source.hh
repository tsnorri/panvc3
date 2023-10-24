/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_SOURCE_HH
#define PANVC3_DISPATCH_EVENTS_SOURCE_HH

#include <atomic>
#include <cstdint>
#include <libbio/assert.hh>
#include <memory>
#include <panvc3/dispatch/fwd.hh>
#include <panvc3/dispatch/queue.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <utility>						// std::forward


namespace panvc3::dispatch::events {
	
	typedef std::int16_t	filter_type;				// Type from struct kevent.
	typedef int				file_descriptor_type;
	
	
	struct source
	{
		virtual ~source() {}
		virtual void disable() = 0;
		virtual void fire() = 0;
	};


	template <typename t_task>
	struct source_tpl_ :	public source
	{
	protected:
		t_task				m_task{};
		queue				*m_queue{};
		std::atomic_bool	m_is_enabled{true};
	
		template <typename t_task_>
		source_tpl_(queue &qq, t_task_ &&tt):
			m_task(std::forward <t_task_>(tt)),
			m_queue(&qq)
		{
		}
	
	public:
		bool is_enabled() const { return m_is_enabled.load(std::memory_order_acquire); } // Possibly relaxed would be enough.
		void disable() override { m_is_enabled.store(false, std::memory_order_release); }
	};
	
	
	template <typename t_self, typename t_task = parametrised <t_self &>::task>
	class source_tpl :	public source_tpl_ <t_task>,
						private std::enable_shared_from_this <t_self>
	{
	public:
		using source_tpl_ <t_task>::source_tpl_;
		using source_tpl_ <t_task>::is_enabled;
		void fire() final;
	
	public:
		void operator()() { if (is_enabled()) m_task(*this); }
	};
	
	
	// Partial specialisation for task without any parameters.
	template <typename t_self>
	struct source_tpl <t_self, task> : public source_tpl_ <task>
	{
		using source_tpl_ <task>::source_tpl_;
		void fire() final;
	};
	
	
	template <typename t_self, typename t_task>
	void source_tpl <t_self, t_task>::fire()
	{
		libbio_assert(this->m_queue);
		// operator() checks is_enabled().
		this->m_queue->async(std::enable_shared_from_this <t_self>::shared_from_this());
	}
	
	
	template <typename t_self>
	void source_tpl <t_self, task>::fire()
	{
		libbio_assert(m_queue);
		if (is_enabled())
			m_task.enqueue_transient_async(*m_queue);
	}
}

#endif
