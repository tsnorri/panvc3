/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_PLATFORM_MANAGER_LINUX_HH
#define PANVC3_DISPATCH_EVENTS_PLATFORM_MANAGER_LINUX_HH

#include <panvc3/dispatch/events/manager.hh>
#include <signal.h>
#include <sys/signalfd.h>	// struct signalfd_siginfo


namespace panvc3::dispatch::events::platform::linux {
	
	struct source_key
	{
		struct fd_tag {};
		struct signal_tag {};
		
		int value{}; // Both signals and file descriptors have type int in POSIX.
		bool is_signal{};
		
		source_key() = default;
		
		constexpr source_key(file_descriptor_type const fd, fd_tag):
			value(fd)
		{
		}
		
		constexpr source_key(signal_type const sig, signal_tag):
			value(sig),
			is_signal(true)
		{
		}
		
		constexpr inline bool operator==(source_key const other) const;
	};
}


template <>
struct std::hash <panvc3::dispatch::events::platform::linux::source_key>
{
	constexpr std::size_t operator()(panvc3::dispatch::events::platform::linux::source_key const &sk) const noexcept
	{
		std::size_t retval(std::abs(sk.value));
		retval <<= 2;
		retval |= std::size_t(sk.value < 0 ? 0x1 : 0x0) << 1;
		retval |= sk.is_signal;
		return retval;
	}
};


namespace panvc3::dispatch::events::platform::linux {

	typedef events::detail::file_handle	file_handle;

	
	struct fd_counter
	{
		std::uint8_t	reader_count{};
		std::uint8_t	writer_count{};
		
		constexpr operator bool() const { return reader_count || writer_count; }
	};
	
	
	struct timer
	{
		file_handle	handle;

		void prepare();
		void start(struct timespec const ts);
		void release() { handle.release(); }
	};


	class signal_monitor
	{
		file_handle m_handle;
		sigset_t m_mask{};
		sigset_t m_original_mask{};
		
	public:
		signal_monitor()
		{
			::sigemptyset(&m_mask);
			::sigemptyset(&m_original_mask);
		}
		
		int file_descriptor() const { return m_handle.fd; }
		bool listen(int sig);
		int unlisten(int sig); // Returns the old file descriptor.
		bool read(struct signalfd_siginfo &buffer);
		void release();
	};


	class event_monitor
	{
		typedef events::manager_base::event_type	event_type;

		file_handle									m_handle;
		std::vector <event_type>					m_events;
		std::mutex									m_mutex;

	public:
		void prepare();

		int file_descriptor() const { return m_handle.fd; }
		std::lock_guard <std::mutex> lock() { return std::lock_guard{m_mutex}; }
		void post(event_type const event);
		std::vector <event_type> const &events() const { return m_events; }
		void clear();
	};
	
	
	class manager final : public events::manager_base
	{
	private:
		typedef std::unordered_multimap <
			source_key,
			std::shared_ptr <events::source>
		>										source_map;
		
		typedef std::unordered_map <
			file_descriptor_type,
			fd_counter
		>										counter_map;
		
	private:
		file_handle								m_epoll_handle;
		
		source_map								m_sources;
		counter_map								m_reader_writer_counts;
		
		signal_monitor							m_signal_monitor;
		event_monitor							m_event_monitor;
		timer									m_timer;
		
		std::mutex								m_mutex{};
		
	public:
		void setup() override;
		void run() override;
		void trigger_event(event_type const evt) override { m_event_monitor.post(evt); }
		
		file_descriptor_source &add_file_descriptor_read_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) override;									// Thread-safe.
		
		file_descriptor_source &add_file_descriptor_write_event_source(
			file_descriptor_type const fd,
			queue &qq,
			file_descriptor_source::task_type tt
		) override;									// Thread-safe.
		
		signal_source &add_signal_event_source(
			signal_type const sig,
			queue &qq,
			signal_source::task_type tt
		) override;									// Thread-safe.
		
		void remove_file_descriptor_event_source(
			file_descriptor_source &es
		) override;									// Thread-safe.
		
		void remove_signal_event_source(
			signal_source &es
		) override;									// Thread-safe.
		
	private:
		void listen_for_fd_events(int fd, bool for_read, bool for_write);
		inline void listen_for_fd_events(file_descriptor_source const &source);
		void schedule_kernel_timer(duration_type const dur);
	};
	
	
	void manager::listen_for_fd_events(file_descriptor_source const &source)
	{
		listen_for_fd_events(
			source.file_descriptor(),
			source.is_read_event_source(),
			source.is_write_event_source()
		);
	}
	
	
	constexpr bool source_key::operator==(source_key const other) const
	{
		if (is_signal != other.is_signal)
			return false;
		
		return value == other.value;
	}
}

#endif
