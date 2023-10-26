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
#include <unistd.h>											// ::close
#include <utility>											// std::greater, std::to_underlying
#include <vector>


namespace panvc3::dispatch::events {
	
	class manager
	{
	public:
		typedef std::uint32_t event_type_;
		enum class event_type : event_type_
		{
			EXIT	= 0x0,
			WAKE_UP	= 0x1
		};
	
	private:
		typedef timer::clock_type		clock_type;
		typedef clock_type::duration	duration_type;
		typedef clock_type::time_point	time_point_type;
		typedef std::greater <>			heap_cmp_type;
	
		constexpr static inline event_type_ const EVENT_MIN{std::to_underlying(event_type::EXIT)};
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
		
		typedef std::shared_ptr <file_descriptor_source>	fd_event_source_ptr;
		typedef std::shared_ptr <signal_source>				signal_event_source_ptr;
		typedef std::vector <timer_entry>					timer_entry_vector;
		typedef std::vector <fd_event_source_ptr>			fd_event_source_vector;
		typedef std::vector <signal_event_source_ptr>		signal_event_source_vector;
		
		struct fd_event_listener_cmp
		{
			bool operator()(fd_event_source_ptr const &lhs, event_listener_identifier_type const rhs) const { return lhs->m_identifier < rhs; }
			bool operator()(event_listener_identifier_type const lhs, fd_event_source_ptr const &rhs) const { return lhs < rhs->m_identifier; }
			bool operator()(fd_event_source_ptr const &lhs, fd_event_source_ptr const &rhs) const { return lhs->m_identifier < rhs->m_identifier; }
		};
		
		struct signal_event_listener_cmp
		{
			bool operator()(signal_event_source_ptr const &lhs, signal_type const rhs) const { return lhs->m_signal < rhs; }
			bool operator()(signal_type const lhs, signal_event_source_ptr const &rhs) const { return lhs < rhs->m_signal; }
			bool operator()(signal_event_source_ptr const &lhs, signal_event_source_ptr const &rhs) const { return lhs->m_signal < rhs->m_signal; }
		};
		
	private:
		std::vector <struct kevent>		m_event_buffer{};
		timer_entry_vector				m_timer_entries{};			// Protected by m_timer_mutex.
		fd_event_source_vector			m_fd_event_sources{};
		signal_event_source_vector		m_signal_event_sources{};
		kqueue_handle					m_kqueue{};
		event_listener_identifier_type	m_fd_event_listener_idx{};
	
		std::mutex						m_timer_mutex{};
		std::mutex						m_fd_es_mutex{};
		std::mutex						m_signal_es_mutex{}; // FIXME: combine some of the mutexes?
	
	private:
		bool check_timers(struct kevent &timer_event);
		inline file_descriptor_source &add_fd_event_source(
			file_descriptor_type const fd,
			filter_type const filter,
			queue &qq,
			file_descriptor_source::task_type &&tt
		);														// Thread-safe.
	
	public:
		void setup(std::size_t const event_buffer_size = 16);	// Not thread-safe.
		void run();												// Not thread-safe.
		
		void trigger_event(event_type const evt);				// Thread-safe.
		
		// NOTE: Currently only one event source of the same type (read or write) is allowed per file descriptor.
		file_descriptor_source &add_file_descriptor_read_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		);														// Thread-safe.
		
		file_descriptor_source &add_file_descriptor_write_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		);														// Thread-safe.
		
		signal_source &add_signal_source(
			signal_type const sig,
			queue &qq,
			signal_source::task_type &&tt
		);														// Thread-safe.
		
		void remove_file_descriptor_event_source(
			file_descriptor_source &fdes
		);														// Thread-safe.
		
		void remove_signal_event_source(
			signal_source &fdes
		);														// Thread-safe.
		
		void schedule_timer(
			timer::duration_type const interval,
			bool repeats,
			queue &qq,
			timer::task_type tt
		);														// Thread-safe.
	};
}

#endif
