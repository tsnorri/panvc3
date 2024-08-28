/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_SIGNAL_SOURCE_HH
#define PANVC3_DISPATCH_EVENTS_SIGNAL_SOURCE_HH

#include <memory>							// std::make_shared, std::shared_ptr
#include <panvc3/dispatch/fwd.hh>
#include <panvc3/dispatch/events/source.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <sys/event.h>						// EVFILT_SIGNAL
#include <utility>							// std::forward


namespace panvc3::dispatch::events {
	
	class manager; // Fwd.
	
	
	class signal_source final : public source_tpl <signal_source>
	{
		friend class manager;
		
	public:
		typedef task_t <signal_source &>	task_type;
	
	private:
		signal_type	m_signal{-1};
	
	public:
		signal_source(
			queue &qq,
			task_type &&tt,
			event_listener_identifier_type const identifier,
			signal_type const sig
		):
			source_tpl(qq, std::move(tt), identifier),
			m_signal(sig)
		{
		}
		
		template <typename ... t_args>
		static std::shared_ptr <signal_source>
		make_shared(t_args && ... args) { return std::make_shared <signal_source>(std::forward <t_args>(args)...); }
		
		signal_type signal() const { return m_signal; }
		
		// Equivalence class in kqueue.
		signal_type ident() const { return m_signal; }
		filter_type filter() const { return EVFILT_SIGNAL; }
	};
	
	typedef signal_source::task_type signal_task;
	
	
	struct sigchld_handler
	{
		virtual ~sigchld_handler() {}
		virtual void child_did_exit_with_nonzero_status(pid_t const pid, int const exit_status, char const *reason) = 0;
		virtual void child_received_signal(pid_t const pid, int const signal_number) = 0;
		virtual void finish_handling(bool const did_report_error) = 0;
	};
	
	void install_sigchld_handler(events::manager &mgr, queue &qq, sigchld_handler &handler);
}


namespace panvc3::dispatch::detail {
	
	template <>
	struct member_callable_target <events::signal_source>
	{
		typedef events::source_tpl <events::signal_source>	type;
	};
}

#endif
