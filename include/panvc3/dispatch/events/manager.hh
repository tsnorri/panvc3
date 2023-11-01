/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_MANAGER_HH
#define PANVC3_DISPATCH_EVENTS_MANAGER_HH

#include <chrono>
#include <cstdint>
#include <memory>
#include <mutex>
#include <panvc3/dispatch/events/file_descriptor_source.hh>
#include <panvc3/dispatch/events/signal_source.hh>
#include <panvc3/dispatch/events/timer.hh>
#include <sys/event.h>										// struct kevent
#include <thread>											// std::jthread
#include <unistd.h>											// ::close
#include <utility>											// std::greater, std::to_underlying
#include <vector>


namespace panvc3::dispatch::events::detail {
	
	void add_listener(file_descriptor_type const kqfd, std::uintptr_t const ident, filter_type const filter, event_listener_identifier_type const identifier);
	void remove_listener(file_descriptor_type const kqfd, std::uintptr_t const ident, filter_type const filter);
	
	
	template <typename t_source>
	class source_handler_tpl
	{
	public:
		typedef t_source						source_type;
		typedef std::shared_ptr <source_type>	source_ptr;
		typedef typename source_type::task_type	task_type;
		
	private:
		struct cmp
		{
			bool operator()(source_ptr const &lhs, event_listener_identifier_type const rhs) const { return lhs->identifier() < rhs; }
			bool operator()(event_listener_identifier_type const lhs, source_ptr const &rhs) const { return lhs < rhs->identifier(); }
			bool operator()(source_ptr const &lhs, source_ptr const &rhs) const { return lhs->identifier() < rhs->identifier(); }
		};
		
	private:
		std::vector <source_ptr>				m_sources;
		event_listener_identifier_type			m_listener_idx{};
		std::mutex								m_mutex{};
		
	public:
		template <typename ... t_args>
		source_type &add_source(file_descriptor_type const kqfd, queue &qq, task_type &&tt, t_args && ... args);	// Thread-safe.
		
		void remove_source(file_descriptor_type const kqfd, source_type &source);									// Thread-safe.
		void fire(event_listener_identifier_type const identifier);													// Thread-safe.
	};
	
	
	typedef source_handler_tpl <file_descriptor_source>	file_descriptor_source_handler;
	typedef source_handler_tpl <signal_source>			signal_source_handler;
}


namespace panvc3::dispatch::events {
	
	class manager
	{
	public:
		typedef std::uint32_t event_type_;
		enum class event_type : event_type_
		{
			STOP	= 0x0,
			WAKE_UP	= 0x1
		};
	
	private:
		typedef timer::clock_type		clock_type;
		typedef clock_type::duration	duration_type;
		typedef clock_type::time_point	time_point_type;
		typedef std::greater <>			heap_cmp_type;
	
		constexpr static inline event_type_ const EVENT_MIN{std::to_underlying(event_type::STOP)};
		constexpr static inline event_type_ const EVENT_MAX{std::to_underlying(event_type::WAKE_UP)};
		constexpr static inline event_type_ const EVENT_COUNT{1 + EVENT_MAX - EVENT_MIN};
	
		struct kqueue_handle
		{
			int fd{-1};
			~kqueue_handle() { if (-1 != fd) ::close(fd); }
		};
	
		struct timer_entry
		{
			time_point_type					firing_time{};
			std::shared_ptr <class timer>	timer{};
		
			bool operator<(timer_entry const &other) const { return firing_time < other.firing_time; }
			bool operator>(timer_entry const &other) const { return firing_time > other.firing_time; }
		};
		
		typedef std::vector <timer_entry>					timer_entry_vector;
		
	private:
		std::vector <struct kevent>				m_event_buffer{};
		timer_entry_vector						m_timer_entries{};			// Protected by m_timer_mutex.
		detail::file_descriptor_source_handler	m_fd_source_handler;
		detail::signal_source_handler			m_signal_source_handler;
		kqueue_handle							m_kqueue{};
		std::mutex								m_timer_mutex{};
	
	private:
		bool check_timers(struct kevent &timer_event);
	
	public:
		void setup(std::size_t const event_buffer_size = 16);											// Not thread-safe.
		void run();																						// Not thread-safe.
		std::jthread start_thread_and_run() { return std::jthread([this]{ run(); }); }					// Not thread-safe.
		
		void trigger_event(event_type const evt);														// Thread-safe.
		void wake_up() { trigger_event(event_type::WAKE_UP); }
		void stop() { trigger_event(event_type::STOP); }
		
		// NOTE: Currently only one event source of the same type (read or write) is allowed per file descriptor.
		// FIXME: Fix the above issue e.g. the same way timers are handled, i.e. with listener <<-> kqueue item.
		file_descriptor_source &add_file_descriptor_read_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) { return m_fd_source_handler.add_source(m_kqueue.fd, qq, std::move(tt), fd, EVFILT_READ); }	// Thread-safe.
		
		file_descriptor_source &add_file_descriptor_write_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) { return m_fd_source_handler.add_source(m_kqueue.fd, qq, std::move(tt), fd, EVFILT_WRITE); };	// Thread-safe.
		
		signal_source &add_signal_source(
			signal_type const sig,
			queue &qq,
			signal_source::task_type &&tt
		) { return m_signal_source_handler.add_source(m_kqueue.fd, qq, std::move(tt), sig); };			// Thread-safe.
		
		void remove_file_descriptor_event_source(
			file_descriptor_source &es
		) { m_fd_source_handler.remove_source(m_kqueue.fd, es); }										// Thread-safe.
		
		void remove_signal_event_source(
			signal_source &es
		) { m_signal_source_handler.remove_source(m_kqueue.fd, es); };									// Thread-safe.
		
		void schedule_timer(
			timer::duration_type const interval,
			bool repeats,
			queue &qq,
			timer::task_type tt
		);																								// Thread-safe.
	};
}


namespace panvc3::dispatch::events::detail {
	
	template <typename t_source>
	template <typename ... t_args>
	auto source_handler_tpl <t_source>::add_source(file_descriptor_type const kqfd, queue &qq, task_type &&tt, t_args && ... args) -> source_type & // Thread-safe.
	{
		source_type *retval{};
		event_listener_identifier_type identifier{};
		
		{
			std::lock_guard lock{this->m_mutex};
			identifier = m_listener_idx++;
			retval = &*m_sources.emplace_back(
				source_type::make_shared(qq, std::move(tt), identifier, std::forward <t_args>(args)...)
			);
			
			// Move to the correct place. (Likely no need b.c. identifier is strictly increasing.)
			auto const it(std::upper_bound(m_sources.begin(), m_sources.end() - 1, identifier, cmp{}));
			std::rotate(it, m_sources.end() - 1, m_sources.end());
		}
		
		auto &retval_(*retval);
		add_listener(kqfd, retval_.ident(), retval_.filter(), retval_.identifier());
		return retval_;
	}
	
	
	template <typename t_source>
	void source_handler_tpl <t_source>::remove_source(file_descriptor_type const kqfd, source_type &source) // Thread-safe.
	{
		auto const ident{source.ident()};
		auto const filter{source.filter()};

		{
			std::lock_guard lock{m_mutex};
			auto const end(m_sources.end());
			auto const it(std::lower_bound(m_sources.begin(), end, ident, cmp{}));
			if (it != end && &source == &(**it))
				m_sources.erase(it);
		}
		
		remove_listener(kqfd, ident, filter);
	}
	
	
	template <typename t_source>
	void source_handler_tpl <t_source>::fire(event_listener_identifier_type const identifier) // Thread-safe.
	{
		std::lock_guard lock{m_mutex};
		auto const end(m_sources.end());
		auto const it(std::lower_bound(m_sources.begin(), end, identifier, cmp{}));
		if (it != end)
		{
			auto &event_source(**it);
			if (event_source.identifier() == identifier)
				event_source.fire();
		}
	}
}

#endif
