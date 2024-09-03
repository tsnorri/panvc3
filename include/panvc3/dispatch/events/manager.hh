/*
 * Copyright (c) 2023-2024 Tuukka Norri
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
#include <panvc3/dispatch/events/synchronous_source.hh>
#include <panvc3/dispatch/events/timer.hh>
#include <thread>											// std::jthread
#include <unistd.h>											// ::close
#include <utility>											// std::greater, std::to_underlying
#include <vector>
#include <unordered_map>


namespace panvc3::dispatch::events::detail {
	
	struct file_handle
	{
		file_descriptor_type fd{-1};

		~file_handle() { if (-1 != fd) ::close(fd); }
		void release() { if (-1 != fd) { ::close(fd); fd = -1; }}
		constexpr operator bool() const { return -1 != fd; }
	};
}


namespace panvc3::dispatch::events {
	
	class manager_base
	{
	public:
		typedef std::uint32_t event_type_;
		enum class event_type : event_type_
		{
			stop	= 0x0,
			wake_up	= 0x1
		};
		
		constexpr static inline event_type_		EVENT_MIN{0x0};
		constexpr static inline event_type_		EVENT_LIMIT{0x2};
		constexpr static inline event_type_		EVENT_COUNT{0x2};
		
	protected:
		typedef timer::clock_type				clock_type;
		typedef clock_type::duration			duration_type;
		typedef clock_type::time_point			time_point_type;
		typedef std::greater <>					heap_cmp_type;
		
		struct timer_entry
		{
			time_point_type						firing_time{};
			std::shared_ptr <class timer>		timer{};
			
			constexpr bool operator<(timer_entry const &other) const { return firing_time < other.firing_time; }
			constexpr bool operator>(timer_entry const &other) const { return firing_time > other.firing_time; }
		};
		
		typedef std::vector <timer_entry>		timer_entry_vector;
		
	protected:
		timer_entry_vector						m_timer_entries{};	// Protected by m_timer_mutex.
		std::mutex								m_timer_mutex{};
		
	public:
		virtual ~manager_base() {}
		
		void setup() = 0;
		void run() = 0;
		std::jthread start_thread_and_run() { return std::jthread([this]{ run(); }); }
		
		void schedule_timer(
			timer::duration_type const interval,
			bool repeats,
			queue &qq,
			timer::task_type tt
		);											// Thread-safe.
		
		file_descriptor_source &add_file_descriptor_read_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) = 0;										// Thread-safe.
		
		file_descriptor_source &add_file_descriptor_write_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) = 0;										// Thread-safe.
		
		signal_source &add_signal_event_source(
			signal_type const sig,
			queue &qq,
			signal_source::task_type tt
		) = 0;										// Thread-safe.
		
		void remove_file_descriptor_event_source(
			file_descriptor_source &es
		) = 0;										// Thread-safe.
		
		void remove_signal_event_source(
			signal_source &es
		) = 0;										// Thread-safe.
		
	protected:
		duration_type check_timers();
	};
}

#endif
