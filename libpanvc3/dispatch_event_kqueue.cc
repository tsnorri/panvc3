/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/dispatch/events/platform/manager_kqueue.hh>
#include <sys/types.h>
#include <sys/event.h>
#include <sys/time.h>


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
	
	
	void add_listener(file_descriptor_type const kqfd, std::uintptr_t const ident, filter_type const filter)
	{
		struct kevent kev{};
		EV_SET(&kev, ident, filter, EV_ADD | EV_ENABLE | EV_RECEIPT | EV_CLEAR, 0, 0, nullptr);
		modify_kqueue(kqfd, &kev, 1);
	}
	
	
	void remove_listener(file_descriptor_type const kqfd, std::uintptr_t const ident, filter_type const filter)
	{
		struct kevent kev{};
		EV_SET(&kev, ident, filter, EV_DELETE | EV_DISABLE | EV_RECEIPT, 0, 0, 0);
		modify_kqueue(kqfd, &kev, 1);
	}
}


namespace panvc3::dispatch::events::platform::kqueue {
	
	void signal_monitor::listen(int const sig)
	{
		::sigaddset(&mask, sig);
		
		if (-1 == ::sigprocmask(SIG_BLOCK, &m_mask, m_is_empty ? &original_mask : nullptr))
			throw std::runtime_error(::strerror(errno));
	}
	
	
	int signal_monitor::unlisten(int sig)
	{
		if (!::sigismember(&m_original_mask, sig))
		{
			sigset_t mask{};
			::sigemptyset(&mask);
			::sigaddset(&mask, sig);
			::sigdelset(&m_mask, sig);
			
			if (-1 == ::sigprocmask(SIG_UNBLOCK, &mask, nullptr))
				throw std::runtime_error(::strerror(errno));
		}
	}
	
	
	void manager::setup() override
	{
		libbio_assert(!m_kqueue_handle);
		m_kqueue_handle.fd = ::kqueue();
		if (-1 == m_kqueue_handle.fd)
			throw std::runtime_error{::strerror(errno)};
		
		// Handle the user events.
		// The kqueue can coalesce events with the same identifier (even if the user data or the lower 24 bits of the flags are different),
		// so we use different identifiers to identify the user event type.
		std::array <struct kevent, EVENT_COUNT> events{};
		for (event_type_ i(EVENT_MIN); i < EVENT_LIMIT; ++i)
			EV_SET(&events[i], i, EVFILT_USER, EV_ADD | EV_ENABLE | EV_RECEIPT, NOTE_FFNOP, 0, 0);
	
		modify_kqueue(m_kqueue_handle.fd, monitored_events.data(), monitored_events.size());
	}
	
	
	void manager::run() override
	{
		int modified_event_count{};
		std::array <struct kevent, 16> event_buffer;
		while (true)
		{
			// Modify the events if needed and wait.
			auto const res(::kevent(
				m_kqueue_handle.fd,
				event_buffer.data(),
				modified_event_count,
				event_buffer.data(),
				event_buffer.size(),
				nullptr
			));
			
			if (-1 == res)
				throw std::runtime_error{::strerror(errno)};
		
			{
				std::lock_guard const lock{m_mutex};
				for (int i(0); i < res; ++i)
				{
					auto const &rec(m_event_buffer[i]);
					switch (rec.filter)
					{
						case EVFILT_USER:
						{
							switch (static_cast <event_type>(rec.ident))
							{
								case event_type::STOP:
									return;
						
								case event_type::WAKE_UP:
									continue;
							}
						}
						
						case EVFILT_READ:
						case EVFILT_WRITE:
						case EVFILT_SIGNAL:
						{
							// source_key having EVFILT_READ (EVFILT_WRITE, EVFILT_SIGNAL) implies that the
							// source monitors the file descriptor for reading (writing, signals).
							int const ident(rec.ident);
							auto range(m_sources.equal_range(source_key{ident, rec.filter}));
							while (range.first != range.second)
							{
								range.first->fire_if_enabled();
								++range.first;
							}
							
							break;
						}
						
						case EVFILT_TIMER:
						default:
							continue;
					}
				}
			}
			
			// Timers
			auto const next_firing_time(check_timers());
			if (events::timer::DURATION_MAX != next_firing_time)
			{
				auto const next_firing_time_ms(chrono::duration_cast <chrono::milliseconds>(next_firing_time).count());
				EV_SET(&timer_event, 0, EVFILT_TIMER, EV_ADD | EV_ENABLE | EV_ONESHOT | EV_CLEAR, 0, next_firing_time_ms, 0); // EV_CLEAR implicit for timers.
			}
		}
	}
	
	
	auto manager::add_file_descriptor_read_event_source(
		file_descriptor_type const fd,
		queue &qq,
		file_descriptor_source::task_type tt
	) -> file_descriptor_source &
	{
		source_key const key{fd, EVFILT_READ};
		bool should_add_listener{};
		
		{
			std::lock_guard const lock{m_mutex};
			should_add_listener = !m_sources.contains(key);
			auto it(m_sources.emplace(
				std::piecewise_construct{},
				std::forward_as_tuple(key),
				std::forward_as_tuple(std::make_shared <file_descriptor_source>(qq, std::move(tt), fd, file_descriptor_source::type::read_source))
			));
		}
		
		if (should_add_listener)
			add_listener(m_kqueue_handle.fd, fd, EVFILT_READ);
	}
	
	
	auto manager::add_file_descriptor_write_event_source(
		file_descriptor_type const fd,
		queue &qq,
		file_descriptor_source::task_type tt
	) -> file_descriptor_source &
	{
		source_key const key{fd, EVFILT_WRITE};
		bool should_add_listener{};
		
		{
			std::lock_guard const lock{m_mutex};
			should_add_listener = !m_sources.contains(key);
			auto it(m_sources.emplace(
				std::piecewise_construct{},
				std::forward_as_tuple(key),
				std::forward_as_tuple(std::make_shared <file_descriptor_source>(qq, std::move(tt), fd, file_descriptor_source::type::write_source))
			));
		}
		
		if (should_add_listener)
			add_listener(m_kqueue_handle.fd, fd, EVFILT_WRITE);
	}
	
	
	auto manager::add_signal_event_source(
		signal_type const sig,
		queue &qq,
		signal_source::task_type tt
	) -> signal_source &
	{
		source_key const key{sig, EVFILT_SIGNAL};
		bool should_add_listener{};
		
		{
			std::lock_guard const lock{m_mutex};
			should_add_listener = !m_sources.contains(key);
			auto const it(m_sources.emplace(
				std::piecewise_construct{},
				std::forward_as_tuple(key),
				std::forward_as_tuple(std::make_shared <signal_source>(qq, std::move(tt), sig))
			));
		}
		
		if (should_add_listener)
			add_listener(m_kqueue_handle.fd, sig, EVFILT_SIGNAL);
	}
	
	
	void manager::remove_event_source(source &es, source_key const key)
	{
		bool should_remove_listener{true};
		
		{
			std::lock_guard const lock{m_mutex};
			auto range(m_sources.equal_range(key));
			while (true)
			{
				if (&*range.first == &es)
				{
					es.disable();
					m_sources.erase(range.first);
					break;
				}
				else
				{
					should_remove_listener = false;
				}
				
				++range.first;
			}
			
			if (range.first != range.second)
				should_remove_listener = false;
		}
		
		if (should_remove_listener)
			remove_listener(m_kqueue_handle.fd, key.value, key.filter);
	}
	
	
	void manager::remove_file_descriptor_event_source(
		file_descriptor_source &es
	)
	{
		filter_type filter{};
		switch (es.file_descriptor_source_type())
		{
			case file_descriptor_source::type::read_source:
				filter = EVFILT_READ;
				break;
			
			case file_descriptor_source::type::write_source:
				filter = EVFILT_WRITE;
				break;
			
			default:
				return;
		}
		
		source_key const key{es.file_descriptor(), filter};
		remove_event_source(es, key);
	}
	
	
	void manager::remove_signal_event_source(
		signal_source &es
	)
	{
		source_key const key{es.signal(), EVFILT_SIGNAL};
		remove_event_source(es, key);
	}
}
