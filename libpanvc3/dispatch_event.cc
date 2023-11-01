/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>				// std::lower_bound
#include <chrono>
#include <cstring>					// ::strerror
#include <stdexcept>
#include <sys/types.h>
#include <sys/event.h>
#include <sys/time.h>
#include <panvc3/dispatch/event.hh>
#include <panvc3/dispatch/task_def.hh>

namespace chrono	= std::chrono;


namespace {
	
	void modify_kqueue(int const kqueue, struct kevent *event_list, int const count)
	{
		struct timespec ts{};
		auto const res(::kevent(kqueue, event_list, count, event_list, count, &ts));
		if (-1 == res)
			throw std::runtime_error{::strerror(errno)};
		
		if (count != res)
			throw std::logic_error{"Unexpected number of receipts"};
		
		for (int i(0); i < res; ++i)
		{
			auto const &rec(event_list[i]);
			
			if (! (EV_ERROR & rec.flags))
				throw std::logic_error{"Unexpected receipt"};
			
			if (0 != rec.data)
				throw std::runtime_error{::strerror(rec.data)};
		}
	}
}


namespace panvc3::dispatch::events::detail {
	
	void add_listener(file_descriptor_type const kqfd, std::uintptr_t const ident, std::uint16_t const filter, event_listener_identifier_type const identifier)
	{
		struct kevent kev{};
		EV_SET(&kev, ident, filter, EV_ADD | EV_ENABLE | EV_RECEIPT | EV_CLEAR, 0, 0, std::bit_cast <void *>(identifier));
		modify_kqueue(kqfd, &kev, 1);
	}
	
	
	void remove_listener(file_descriptor_type const kqfd, std::uintptr_t const ident, std::uint16_t const filter)
	{
		struct kevent kev{};
		EV_SET(&kev, ident, filter, EV_DELETE | EV_DISABLE | EV_RECEIPT, 0, 0, 0);
		modify_kqueue(kqfd, &kev, 1);
	}
}


namespace panvc3::dispatch::events {
	
	bool manager::check_timers(struct kevent &timer_event)
	{
		duration_type next_firing_time{timer::DURATION_MAX};
		
		{
			std::lock_guard lock{m_timer_mutex};
		
			if (m_timer_entries.empty())
				return false;
			
			do
			{
				auto &next_timer(m_timer_entries.front());
				if (!next_timer.timer->is_enabled())
				{
					std::pop_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
					m_timer_entries.pop_back(); // Safe b.c. the timer uses a shared_ptr_task.
					continue;
				}
				
				auto const duration(next_timer.firing_time - clock_type::now());
				auto const duration_ms(chrono::duration_cast <chrono::milliseconds>(duration).count());
				if (duration_ms <= 0)
				{
					// FIXME: Verify that we could fire while not holding m_timer_mutex?
					// FIXME: The concurrent queue could benefit from sorting the timers by queue and then bulk inserting the operations.
					next_timer.timer->fire();
				
					// If the timer repeats, re-schedule.
					if (next_timer.timer->repeats())
					{
						next_timer.firing_time += next_timer.timer->interval();
						std::make_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
					}
					else
					{
						// m_timers may be empty after popping.
						std::pop_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
						m_timer_entries.pop_back();
					}
				}
				else
				{
					// m_timer_entries not empty.
					next_firing_time = m_timer_entries.front().firing_time - clock_type::now();
					break;
				}
			}
			while (!m_timer_entries.empty());
		}
	
		// Schedule the kernel timer.
		if (timer::DURATION_MAX != next_firing_time)
		{
			auto const next_firing_time_ms(chrono::duration_cast <chrono::milliseconds>(next_firing_time).count());
			EV_SET(&timer_event, 0, EVFILT_TIMER, EV_ADD | EV_ENABLE | EV_ONESHOT | EV_CLEAR, 0, next_firing_time_ms, 0); // EV_CLEAR implicit for timers.
			return true;
		}
	
		return false;
	}


	void manager::setup(std::size_t const event_buffer_size)
	{
		m_event_buffer.resize(event_buffer_size);
		m_kqueue.fd = ::kqueue();
		if (-1 == m_kqueue.fd)
			throw std::runtime_error{::strerror(errno)};
	
		std::array <struct kevent, EVENT_COUNT> monitored_events{};
	
		// Handle the user events.
		// The kqueue can coalesce events with the same identifier (even if the user data or the lower 24 bits of the flags are different),
		// so we use different identifiers to identify the user event type.
		// Add signal handling here if needed.
		for (event_type_ i(EVENT_MIN); i < 1 + EVENT_MAX; ++i)
			EV_SET(&monitored_events[i], i, EVFILT_USER, EV_ADD | EV_ENABLE | EV_RECEIPT, NOTE_FFNOP, 0, 0);
	
		modify_kqueue(m_kqueue.fd, monitored_events.data(), monitored_events.size());
	}


	void manager::trigger_event(event_type const evt)
	{
		struct kevent kev{};
		struct timespec ts{};
		EV_SET(&kev, std::to_underlying(evt), EVFILT_USER, 0, NOTE_TRIGGER, 0, nullptr);
		auto const res(::kevent(m_kqueue.fd, &kev, 1, nullptr, 0, &ts));
		if (res < 0)
			throw std::runtime_error{::strerror(errno)};
	}


	void manager::schedule_timer(
		timer::duration_type const interval,
		bool repeats,
		queue &qq,
		timer::task_type tt
	)
	{
		{
			std::lock_guard lock(m_timer_mutex);
			auto const now(clock_type::now());
		
			// Insert to the correct position.
			m_timer_entries.emplace_back(now + interval, std::make_shared <timer>(qq, std::move(tt), interval, repeats));
			std::push_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
		}
	
		trigger_event(event_type::WAKE_UP);
	}


	void manager::run()
	{
		int modified_event_count{};
	
		while (true)
		{
			// Modify the events if needed and wait.
			auto const res(::kevent(m_kqueue.fd, m_event_buffer.data(), modified_event_count, m_event_buffer.data(), m_event_buffer.size(), nullptr));
		
			if (-1 == res)
				throw std::runtime_error{::strerror(errno)};
		
			for (int i(0); i < res; ++i)
			{
				auto const &rec(m_event_buffer[i]);
				switch (rec.filter)
				{
					case EVFILT_USER:
					{
						switch (static_cast <event_type>(rec.ident))
						{
							case event_type::EXIT:
								return;
						
							case event_type::WAKE_UP:
								continue;
						}
					}
					
					case EVFILT_READ:
					case EVFILT_WRITE:
						m_fd_source_handler.fire(std::bit_cast <event_listener_identifier_type>(rec.udata));
						continue;
					
					case EVFILT_SIGNAL:
						m_signal_source_handler.fire(std::bit_cast <event_listener_identifier_type>(rec.udata));
						continue;
					
					case EVFILT_TIMER:
					default:
						continue;
				}
			}
			
			modified_event_count = 0;
			modified_event_count += check_timers(m_event_buffer[modified_event_count]);
		}
	}
	
	
	void install_sigchld_handler(events::manager &mgr, queue &qq, sigchld_handler &handler)
	{
		mgr.add_signal_source(SIGCHLD, qq, signal_source::task_type::from_lambda([&handler](signal_source const &source){
			bool did_report_error(false);
			while (true)
			{
				int status(0);
				auto const pid(::waitpid(-1, &status, WNOHANG | WUNTRACED));
				if (pid <= 0)
					break;
		
				if (WIFEXITED(status))
				{
					auto const exit_status(WEXITSTATUS(status));
					if (0 != exit_status)
					{
						did_report_error = true;
						char const *reason{nullptr};
						switch (exit_status)
						{
							case 127:
								reason = "command not found";
								break;
					
							case 126:
								reason = "command invoked cannot execute";
								break;
					
							case 69: // EX_UNAVAILABLE
								reason = "service unavailable";
								break;
					
							case 71: // EX_OSERR
								reason = "unknown error from execvp()";
								break;
					
							case 74: // EX_IOERR
								reason = "an I/O error occurred";
								break;
					
							default:
								break;
						}
						
						handler.child_did_exit_with_nonzero_status(pid, exit_status, reason);
					}
				}
				else if (WIFSIGNALED(status))
				{
					did_report_error = true;
					auto const signal_number(WTERMSIG(status));
					handler.child_received_signal(pid, signal_number);
				}
			}
			
			handler.finish_handling(did_report_error);
		}));
	}
}
