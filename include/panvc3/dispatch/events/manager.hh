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
		
		typedef std::shared_ptr <timer>			timer_ptr;
		
	protected:
		typedef timer::clock_type				clock_type;
		typedef clock_type::duration			duration_type;
		typedef clock_type::time_point			time_point_type;
		typedef std::greater <>					heap_cmp_type;
		
		struct timer_entry
		{
			time_point_type						firing_time{};
			timer_ptr							timer{};
			
			constexpr bool operator<(timer_entry const &other) const { return firing_time < other.firing_time; }
			constexpr bool operator>(timer_entry const &other) const { return firing_time > other.firing_time; }
		};
		
		typedef std::vector <timer_entry>		timer_entry_vector;
		
	protected:
		timer_entry_vector						m_timer_entries{};	// Protected by m_timer_mutex.
		std::mutex								m_timer_mutex{};
		std::atomic_bool						m_is_running_worker{};
		
	public:
		virtual ~manager_base() {}
		
		virtual void setup() = 0;
		void run();
		void start_thread_and_run(std::jthread &);
		
		virtual void trigger_event(event_type const evt) = 0;
		void stop() { trigger_event(event_type::stop); }
		
		// The fact that the non-timer sources are stored in shared pointers
		// is an implementation detail b.c. they are guaranteed to persist
		// until one of the remove_*_event_source() member functions is called.
		// In case of timers this is not an implementation detail b.c. they are
		// in fact removed either after having been disabled or, for non-repeating
		// timers, after having been fired.
		
		timer_ptr schedule_timer(
			timer::duration_type const interval,
			bool repeats,
			queue &qq,
			timer::task_type tt
		);											// Thread-safe.
		
		virtual file_descriptor_source &add_file_descriptor_read_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) = 0;										// Thread-safe.
		
		virtual file_descriptor_source &add_file_descriptor_write_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) = 0;										// Thread-safe.
		
		virtual signal_source &add_signal_event_source(
			signal_type const sig,
			queue &qq,
			signal_source::task_type tt
		) = 0;										// Thread-safe.
		
		virtual void remove_file_descriptor_event_source(
			file_descriptor_source &es
		) = 0;										// Thread-safe.
		
		virtual void remove_signal_event_source(
			signal_source &es
		) = 0;										// Thread-safe.
		
	protected:
		void stop_and_wait();
		virtual void run_() = 0;
		duration_type check_timers();
	};
}

#endif
