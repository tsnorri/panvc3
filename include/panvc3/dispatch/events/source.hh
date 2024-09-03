/*
 * Copyright (c) 2023-2024 Tuukka Norri
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


namespace panvc3::dispatch::events::detail {
	
	class enabled_status
	{
		std::atomic_bool	m_status{true};
		
	public:
		operator bool() const { return m_status.load(std::memory_order_acquire); }
		void disable() { m_status.store(false, std::memory_order_release); }
	};
}


namespace panvc3::dispatch::events {
	
	typedef int				file_descriptor_type;
	typedef int				signal_type;
	
	
	struct source
	{
		virtual ~source() {}
		virtual event_listener_identifier_type identifier() const = 0;
		virtual bool is_enabled() const = 0;
		virtual void disable() = 0;
		virtual void fire() = 0;
		virtual void fire_if_enabled() = 0;
		
		// Needed by the Linux implementation.
		virtual bool is_read_event_source() const { return false; }
		virtual bool is_write_event_source() const { return false; }
	};
	
	
	template <typename t_task, bool t_needs_queue>
	struct source_tpl_ : public source
	{
	protected:
		t_task							m_task{};
		detail::enabled_status			m_is_enabled;
		
		template <typename t_task_>
		source_tpl_(t_task_ &&tt):
			m_task(std::forward <t_task_>(tt))
		{
		}
		
		bool is_enabled() const final { return bool(m_is_enabled); }
		void disable() final { m_is_enabled.disable(); }
		void fire_if_enabled() final { if (is_enabled()) fire(); }
	};
	
	
	template <typename t_task>
	struct source_tpl_ <t_task, true> : public source
	{
	protected:
		t_task							m_task{};
		queue							*m_queue{};
		detail::enabled_status			m_is_enabled;
		
		template <typename t_task_>
		source_tpl_(queue &qq, t_task_ &&tt):
			m_task(std::forward <t_task_>(tt)),
			m_queue(&qq)
		{
		}
		
		bool is_enabled() const final { return bool(m_is_enabled); }
		void disable() final { m_is_enabled.disable(); }
		void fire_if_enabled() final { if (is_enabled()) fire(); }
	};


	template <typename t_self, bool t_needs_queue = true, typename t_task = task_t <t_self &>>
	class source_tpl :	public source_tpl_ <t_task, t_needs_queue>,
						private std::enable_shared_from_this <t_self>
	{
	public:
		using source_tpl_ <t_task>::source_tpl_;
		void fire() final;
		void operator()() { if (this->is_enabled()) this->m_task(static_cast <t_self &>(*this)); }
	};
	
	
	// Partial specialisation for task without any parameters.
	template <typename t_self, bool t_needs_queue>
	struct source_tpl <t_self, t_needs_queue, task> : public source_tpl_ <task, t_needs_queue>
	{
		using source_tpl_ <task>::source_tpl_;
		void fire() final;
		void operator()() { if (this->is_enabled()) this->m_task(); }
	};
	
	
	template <typename t_self, bool t_needs_queue, typename t_task>
	void source_tpl <t_self, t_needs_queue, t_task>::fire()
	{
		if constexpr (t_needs_queue)
			(*this)();
		else
		{
			libbio_assert(this->m_queue);
			// operator() checks is_enabled().
			this->m_queue->async(std::enable_shared_from_this <t_self>::shared_from_this());
		}
	}
	
	
	template <typename t_self, bool t_needs_queue>
	void source_tpl <t_self, t_needs_queue, task>::fire()
	{
		if constexpr (t_needs_queue)
			(*this)();
		else
		{
			libbio_assert(m_queue);
			if (is_enabled())
				m_task.enqueue_transient_async(*m_queue);
		}
	}
}

#endif
