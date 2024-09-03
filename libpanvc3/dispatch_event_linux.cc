/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <chrono>
#include <panvc3/dispatch/events/platform/manager_linux.hh>
#include <range/v3/view/take_exactly.hpp>
#include <sys/epoll.h>
#include <sys/eventfd.h>
#include <sys/timerfd.h>

namespace chrono	= std::chrono;
namespace rsv		= ranges::views;


namespace {

	void add_read_event_listener(int const epoll_fd, int const fd, int const user_fd)
	{
		struct epoll_event ev{};
		ev.data.fd = user_fd;
		ev.events |= EPOLLIN;

		if (-1 == ::epoll_ctl(epoll_fd, EPOLL_CTL_ADD, fd, &ev))
			throw std::runtime_error(::strerror(errno));
	}
	
	
	void remove_fd_event_listener(int const epoll_fd, int const fd)
	{
		if (-1 == ::epoll_ctl(epoll_fd, EPOLL_CTL_DEL, fd, nullptr))
			throw std::runtime_error(::strerror(errno));
	}
}


namespace panvc3::dispatch::events::platform::linux {
	
	void manager::setup()
	{
		libbio_assert(!m_epoll_handle);
		
		m_epoll_handle.fd = ::epoll_create1(EPOLL_CLOEXEC);
		if (!m_epoll_handle)
			throw std::runtime_error(::strerror(errno));
		
		// Timers
		m_timer.prepare();
		add_read_event_listener(m_epoll_handle.fd, m_timer.handle.fd, -1);
		
		// User events
		m_event_monitor.prepare();
		add_read_event_listener(m_epoll_handle.fd, m_event_monitor.file_descriptor(), -1);
	}
	
	
	void manager::run()
	{
		std::array <struct epoll_event, 16> events;
		while (true)
		{
			auto const count(::epoll_wait(m_epoll_handle.fd, events.data(), events.size(), -1));
			if (-1 == count)
				throw std::runtime_error(::strerror(errno));
			
			// User events
			{
				auto const lock(m_event_monitor.lock());
				for (auto const ev : m_event_monitor.events())
				{
					switch (ev)
					{
						case event_type::stop:
							return;
						
						case event_type::wake_up:
							continue;
					}
				}
				
				m_event_monitor.clear();
			}
			
			// File descriptors
			{
				std::lock_guard const lock{m_mutex};
				for (auto const &event : events | rsv::take_exactly(count))
				{
					// Note that we always require ev.data’s fd member to be initialised.
					auto const it(m_sources.find(source_key{event.data.fd, source_key::fd_tag{}}));
					if (m_sources.end() == it)
						continue;
					
					auto &source(*it->second);
					if ((EPOLLIN & event.events) && source.is_read_event_source())
						source.fire_if_enabled();
					
					if ((EPOLLOUT & event.events) && source.is_write_event_source())
						source.fire_if_enabled();
				}
			}
			
			// Timers
			auto const next_firing_time(check_timers());
			if (events::timer::DURATION_MAX != next_firing_time)
				schedule_kernel_timer(next_firing_time);
		}
	}
	
	
	void manager::schedule_kernel_timer(duration_type dur)
	{
		auto const seconds(chrono::duration_cast <chrono::seconds>(dur));
		dur -= seconds;
		auto const nanoseconds(chrono::duration_cast <chrono::nanoseconds>(dur));
		
		struct timespec const ts{.tv_sec = seconds.count(), .tv_nsec = nanoseconds.count()};
		m_timer.start(ts);
	}
	
	
	auto manager::add_file_descriptor_read_event_source(
		file_descriptor_type const fd,
		queue &qq,
		file_descriptor_source::task_type tt
	) -> file_descriptor_source &				// Thread-safe.
	{
		auto ptr(std::make_shared <file_descriptor_source>(qq, std::move(tt), fd, file_descriptor_source::type::read_source));
		auto &ref(*ptr);

		std::lock_guard const lock{m_mutex};
		auto it(m_sources.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(fd, source_key::fd_tag{}),
			std::forward_as_tuple(std::move(ptr))
		));
		
		listen_for_fd_events(ref);
		return ref;
	}
	
	
	auto manager::add_file_descriptor_write_event_source(
		file_descriptor_type const fd,
		queue &qq,
		file_descriptor_source::task_type tt
	) -> file_descriptor_source &				// Thread-safe.
	{
		auto ptr(std::make_shared <file_descriptor_source>(qq, std::move(tt), fd, file_descriptor_source::type::write_source));
		auto &ref(*ptr);

		std::lock_guard const lock{m_mutex};
		auto const it(m_sources.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(fd, source_key::fd_tag{}),
			std::forward_as_tuple(std::move(ptr))
		));
		
		listen_for_fd_events(ref);
		return ref;
	}
	
	
	void manager::remove_file_descriptor_event_source(
		file_descriptor_source &es
	)											// Thread-safe.
	{
		std::lock_guard const lock{m_mutex};
		auto const fd(es.file_descriptor());
		source_key const key{fd, source_key::fd_tag{}};
		auto range(m_sources.equal_range(key));
		while (range.first != range.second)
		{
			if (&*range.first->second == &es)
			{
				es.disable();
				
				auto it(m_reader_writer_counts.find(fd));
				auto &count(it->second);
				libbio_assert_neq(m_reader_writer_counts.end(), it);
				
				if (es.is_read_event_source())
				{
					libbio_assert_lt(0, count.reader_count);
					--count.reader_count;
				}
				
				if (es.is_write_event_source())
				{
					libbio_assert_lt(0, count.writer_count);
					--count.writer_count;
				}
				
				if (!count)
				{
					// No listeners left for fd.
					m_reader_writer_counts.erase(fd);
					remove_fd_event_listener(m_epoll_handle.fd, fd);
				}
				
				m_sources.erase(range.first);
				break;
			}
			
			++range.first;
		}
	}
	
	
	auto manager::add_signal_event_source(
		signal_type const sig,
		queue &qq,
		signal_source::task_type tt
	) -> signal_source &						// Thread-safe.
	{
		source_key const key{sig, source_key::signal_tag{}};
		auto ptr(std::make_shared <signal_source>(qq, std::move(tt), sig));
		auto &ref(*ptr);

		std::lock_guard const lock{m_mutex};
		libbio_assert(!m_sources.contains(key)); // Our signalfd handle does not support multiple observers of the same signal.
		auto const it(m_sources.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(key),
			std::forward_as_tuple(std::move(ptr))
		));
		
		if (m_signal_monitor.listen(sig))
		{
			auto const sig_fd(m_signal_monitor.file_descriptor());
			auto const it(m_sources.emplace(
				std::piecewise_construct,
				std::forward_as_tuple(sig_fd, source_key::fd_tag{}),
				std::forward_as_tuple(std::make_shared <synchronous_source>([this](synchronous_source &){
					struct signalfd_siginfo buffer{};
					while (m_signal_monitor.read(buffer))
					{
						source_key const key{int(buffer.ssi_signo), source_key::signal_tag{}}; // Perhaps ssi_signo is std::uint32_t to get nicer struct layout?
						auto range(m_sources.equal_range(key)); // Extra safety.
						while (range.first != range.second)
						{
							auto &source(*range.first->second);
							source.fire_if_enabled();
							++range.first;
						}
					}
				}))
			));
			
			listen_for_fd_events(sig_fd, true, false);
		}

		return ref;
	}
	
	
	void manager::remove_signal_event_source(
		signal_source &es
	)											// Thread-safe.
	{
		std::lock_guard const lock{m_mutex};
		source_key const key{es.signal(), source_key::signal_tag{}};
		auto range(m_sources.equal_range(key)); // Extra safety.
		while (range.first != range.second)
		{
			if (&*range.first->second == &es)
			{
				auto const res(m_signal_monitor.unlisten(es.signal()));
				if (-1 != res)
				{
					// Remove the signalfd source.
					source_key const key_(res, source_key::fd_tag{});
					auto range_(m_sources.equal_range(key_));
					while (range_.first != range_.second)
					{
						range_.first->second->disable();
						m_sources.erase(range_.first);
					}
					
					// Since unlisten() called ::close(), epoll_ctl need not be called to delete the entry.
					// (We have not copied the file descriptor; it has close-on-exec and we don’t use ::dup() etc.)
				}
				
				es.disable();
				m_sources.erase(range.first);
				break;
			}
			
			++range.first;
		}
	}
	
	
	void manager::listen_for_fd_events(int const fd, bool const for_read, bool const for_write)
	{
		// Caller needs to acquire m_mutex.
		
		auto &counts(m_reader_writer_counts[fd]);
		auto counts_(counts);
		bool should_modify{};
		
		struct epoll_event ev{};
		ev.data.fd = fd;
		
		if (for_read)
		{
			++counts_.reader_count;
			ev.events |= EPOLLIN;
			should_modify = true;
		}
		
		if (for_write)
		{
			++counts_.writer_count;
			ev.events |= EPOLLOUT;
			should_modify = true;
		}
		
		if (counts) // Check the old counts.
		{
			// Modify an existing listener.
			if (should_modify && -1 == ::epoll_ctl(m_epoll_handle.fd, EPOLL_CTL_MOD, fd, &ev))
				throw std::runtime_error(::strerror(errno));
		}
		else
		{
			// Add a listener.
			if (-1 == ::epoll_ctl(m_epoll_handle.fd, EPOLL_CTL_ADD, fd, &ev))
				throw std::runtime_error(::strerror(errno));
		}
		
		counts = counts_;
	}
	
	
	void timer::prepare()
	{
		if (!handle)
			handle.fd = timerfd_create(CLOCK_MONOTONIC, TFD_CLOEXEC);
	}
	
	
	void timer::start(struct timespec const ts)
	{
		struct itimerspec const its{ .it_interval = { .tv_sec = 0, .tv_nsec = 0 }, .it_value = ts};
		if (0 != timerfd_settime(handle.fd, 0, &its, nullptr))
			throw std::runtime_error(::strerror(errno));
	}
	
	
	bool signal_monitor::listen(int const sig)
	{
		::sigaddset(&m_mask, sig);
		
		if (-1 == ::sigprocmask(SIG_BLOCK, &m_mask, m_handle ? nullptr : &m_original_mask))
			throw std::runtime_error(::strerror(errno));
		
		auto const res(::signalfd(m_handle.fd, &m_mask, SFD_NONBLOCK | SFD_CLOEXEC));
		if (-1 == res)
			throw std::runtime_error(::strerror(errno));
		
		if (-1 == m_handle.fd)
		{
			m_handle.fd = res;
			return true;
		}
		
		return false;
	}
	
	
	int signal_monitor::unlisten(int sig)
	{
		int retval{-1};
		sigset_t mask{};
		::sigemptyset(&mask);
		::sigaddset(&mask, sig);
		::sigdelset(&m_mask, sig);
		
		// Check if we still have signals to monitor.
		if (::sigisemptyset(&m_mask))
		{
			retval = m_handle.fd;
			m_handle.release();
		}
		else if (-1 == ::signalfd(m_handle.fd, &m_mask, 0))
		{
			throw std::runtime_error(::strerror(errno));
		}
		
		if (!::sigismember(&m_original_mask, sig) && -1 == ::sigprocmask(SIG_UNBLOCK, &mask, nullptr))
			throw std::runtime_error(::strerror(errno));
		
		return retval;
	}
	
	
	void signal_monitor::release()
	{
		m_handle.release();
		if (-1 == ::sigprocmask(SIG_SETMASK, &m_original_mask, nullptr))
			throw std::runtime_error(::strerror(errno));
	}


	bool signal_monitor::read(struct signalfd_siginfo &buffer)
	{
		if (-1 == ::read(m_handle.fd, &buffer, sizeof(struct signalfd_siginfo)))
		{
			if (EAGAIN == errno)
				return false;

			throw std::runtime_error(::strerror(errno));
		}

		return true;
	}


	void event_monitor::prepare()
	{
		m_handle.fd = ::eventfd(0, EFD_CLOEXEC | EFD_NONBLOCK);
		if (!m_handle)
			throw std::runtime_error(::strerror(errno));
	}


	void event_monitor::post(manager::event_type const event)
	{
		auto const ll(lock());
		m_events.push_back(event);

		std::uint64_t val{1};
		while (true)
		{
			switch (::write(m_handle.fd, &val, sizeof(std::uint64_t)))
			{
				case sizeof(std::uint64_t):
					return;

				case -1:
				{
					if (EAGAIN == errno)
					{
						// Should not happen?
						struct timespec const ts{.tv_sec = 0, .tv_nsec = 50};
						::nanosleep(&ts, nullptr);
						continue;
					}

					throw std::runtime_error(::strerror(errno));
				}

				default:
					throw std::runtime_error("Unexpected number of bytes written");
			};
		}
	}
	
	
	void event_monitor::clear()
	{
		m_events.clear();
		
		// Clear the fd.
		std::uint64_t buffer{};
		::read(m_handle.fd, &buffer, sizeof(std::uint64_t));
	}
}
