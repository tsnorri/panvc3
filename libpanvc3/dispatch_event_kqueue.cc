/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <chrono>
#include <panvc3/dispatch/events/platform/manager_kqueue.hh>
#include <panvc3/dispatch/task_def.hh>
#include <signal.h>
#include <sys/types.h>
#include <sys/event.h>
#include <sys/time.h>

namespace chrono	= std::chrono;
namespace events	= panvc3::dispatch::events;
namespace kq		= panvc3::dispatch::events::platform::kqueue;


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
	
	
	void add_listener(events::file_descriptor_type const kqfd, std::uintptr_t const ident, kq::filter_type const filter)
	{
		struct kevent kev{};
		EV_SET(&kev, ident, filter, EV_ADD | EV_ENABLE | EV_RECEIPT | EV_CLEAR, 0, 0, nullptr);
		modify_kqueue(kqfd, &kev, 1);
	}
	
	
	void remove_listener(events::file_descriptor_type const kqfd, std::uintptr_t const ident, kq::filter_type const filter)
	{
		struct kevent kev{};
		EV_SET(&kev, ident, filter, EV_DELETE | EV_DISABLE | EV_RECEIPT, 0, 0, 0);
		modify_kqueue(kqfd, &kev, 1);
	}
}


namespace panvc3::dispatch::events::platform::kqueue {
	
	void signal_mask::add(int const sig)
	{
		typedef struct sigaction sigaction_;
		auto res(m_old_actions.emplace(sig, sigaction_{}));
		libbio_assert(res.second);
		auto &old_action(res.first->second);
		
		sigaction_ new_action{};
		new_action.sa_handler = SIG_IGN;
		sigemptyset(&new_action.sa_mask);
		sigaddset(&new_action.sa_mask, sig);
		
		if (-1 == ::sigaction(sig, &new_action, &old_action))
			throw std::runtime_error(::strerror(errno));
	}
	
	
	void signal_mask::remove(int sig)
	{
		auto it(m_old_actions.find(sig));
		libbio_assert_neq(it, m_old_actions.end());
		if (-1 == ::sigaction(sig, &it->second, nullptr))
			throw std::runtime_error(::strerror(errno));
	}
	
	
	void signal_mask::remove_all_()
	{
		auto it(m_old_actions.begin());
		auto const end(m_old_actions.end());
		while (it != end)
		{
			auto const sig(it->first);
			auto &act(it->second);
			if (-1 != ::sigaction(sig, &act, nullptr))
			{
				auto const it_(it);
				++it;
				m_old_actions.erase(it_);
			}
			else
			{
				++it;
			}
		}
	}
	
	
	void signal_mask::remove_all()
	{
		auto it(m_old_actions.begin());
		auto const end(m_old_actions.end());
		while (it != end)
		{
			auto const sig(it->first);
			auto &act(it->second);
			if (-1 != ::sigaction(sig, &act, nullptr))
			{
				auto const it_(it);
				++it;
				m_old_actions.erase(it_);
			}
			else
			{
				throw std::runtime_error(::strerror(errno));
			}
		}
	}
	
	
	void manager::setup()
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
	
		modify_kqueue(m_kqueue_handle.fd, events.data(), events.size());
	}
	
	
	void manager::run_()
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
			
			modified_event_count = 0;
			
			if (-1 == res)
				throw std::runtime_error{::strerror(errno)};
		
			{
				std::lock_guard const lock{m_mutex};
				for (int i(0); i < res; ++i)
				{
					auto const &rec(event_buffer[i]);
					switch (rec.filter)
					{
						case EVFILT_USER:
						{
							switch (static_cast <event_type>(rec.ident))
							{
								case event_type::stop:
									return;
						
								case event_type::wake_up:
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
								range.first->second->fire_if_enabled();
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
				auto &timer_event(event_buffer[modified_event_count]); // Make sure that modified_event_count is less than size in case more modifications are needed.
				EV_SET(&timer_event, 0, EVFILT_TIMER, EV_ADD | EV_ENABLE | EV_ONESHOT | EV_CLEAR, 0, next_firing_time_ms, 0); // EV_CLEAR implicit for timers.
				++modified_event_count;
			}
		}
	}
	
	
	void manager::trigger_event(event_type const evt)
	{
		struct kevent kev{};
		struct timespec ts{};
		EV_SET(&kev, std::to_underlying(evt), EVFILT_USER, 0, NOTE_TRIGGER, 0, nullptr);
		auto const res(::kevent(m_kqueue_handle.fd, &kev, 1, nullptr, 0, &ts));
		if (res < 0)
			throw std::runtime_error{::strerror(errno)};
	}
	
	
	auto manager::add_file_descriptor_read_event_source(
		file_descriptor_type const fd,
		queue &qq,
		file_descriptor_source::task_type tt
	) -> file_descriptor_source &
	{
		source_key const key{fd, EVFILT_READ};
		bool should_add_listener{};
		
		auto &retval([&] -> file_descriptor_source & {
			std::lock_guard const lock{m_mutex};
			should_add_listener = !m_sources.contains(key);
			auto ptr(std::make_shared <file_descriptor_source>(qq, std::move(tt), fd, file_descriptor_source::type::read_source));
			auto &retval(*ptr);
			auto it(m_sources.emplace(
				std::piecewise_construct,
				std::forward_as_tuple(key),
				std::forward_as_tuple(std::move(ptr))
			));
			
			return retval;
		}());
		
		if (should_add_listener)
			add_listener(m_kqueue_handle.fd, fd, EVFILT_READ);
		
		return retval;
	}
	
	
	auto manager::add_file_descriptor_write_event_source(
		file_descriptor_type const fd,
		queue &qq,
		file_descriptor_source::task_type tt
	) -> file_descriptor_source &
	{
		source_key const key{fd, EVFILT_WRITE};
		bool should_add_listener{};
		
		auto &retval([&] -> file_descriptor_source & {
			std::lock_guard const lock{m_mutex};
			should_add_listener = !m_sources.contains(key);
			auto ptr(std::make_shared <file_descriptor_source>(qq, std::move(tt), fd, file_descriptor_source::type::write_source));
			auto &retval(*ptr);
			auto it(m_sources.emplace(
				std::piecewise_construct,
				std::forward_as_tuple(key),
				std::forward_as_tuple(std::move(ptr))
			));
			
			return retval;
		}());
		
		if (should_add_listener)
			add_listener(m_kqueue_handle.fd, fd, EVFILT_WRITE);
		
		return retval;
	}
	
	
	auto manager::add_signal_event_source(
		signal_type const sig,
		queue &qq,
		signal_source::task_type tt
	) -> signal_source &
	{
		source_key const key{sig, EVFILT_SIGNAL};
		
		auto &retval([&]() -> signal_source & {
			std::lock_guard const lock{m_mutex};
			libbio_assert(!m_sources.contains(key)); // Only one source per signal is currently supported.
			auto ptr(std::make_shared <signal_source>(qq, std::move(tt), sig));
			auto &retval(*ptr);
			auto const it(m_sources.emplace(
				std::piecewise_construct,
				std::forward_as_tuple(key),
				std::forward_as_tuple(std::move(ptr))
			));
			
			return retval;
		}());
		
		add_listener(m_kqueue_handle.fd, sig, EVFILT_SIGNAL);
		return retval;
	}
	
	
	auto manager::remove_event_source(source &es, source_key const key) -> remove_event_source_return_type
	{
		remove_event_source_return_type retval{m_mutex, true};
		
		auto range(m_sources.equal_range(key));
		while (true)
		{
			if (&*range.first->second == &es)
			{
				es.disable();
				m_sources.erase(range.first);
				
				++range.first;
				if (range.first != range.second)
					retval.second = false;
				
				return retval;
			}
			else
			{
				retval.second = false;
			}
			
			++range.first;
		}
		
		return retval;
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
		bool should_remove_listener{};
		
		{
			auto const res(remove_event_source(es, key));
			should_remove_listener = res.second;
		}
		
		// es no longer valid.
		if (should_remove_listener)
			remove_listener(m_kqueue_handle.fd, key.value, filter);
	}
	
	
	void manager::remove_signal_event_source(signal_source &es)
	{
		source_key const key{es.signal(), EVFILT_SIGNAL};
		bool should_remove_listener{};
		
		{
			auto const res(remove_event_source(es, key));
			libbio_assert(res.second);
		}
		
		// es no longer valid.
		remove_listener(m_kqueue_handle.fd, key.value, EVFILT_SIGNAL);
	}
}
