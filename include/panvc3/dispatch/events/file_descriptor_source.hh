/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_FILE_DESCRIPTOR_SOURCE_HH
#define PANVC3_DISPATCH_EVENTS_FILE_DESCRIPTOR_SOURCE_HH

#include <memory>							// std::make_shared, std::shared_ptr
#include <panvc3/dispatch/fwd.hh>
#include <panvc3/dispatch/events/source.hh>
#include <panvc3/dispatch/task_decl.hh>
#include <utility>							// std::forward


namespace panvc3::dispatch::events {
	
	class file_descriptor_source final : public source_tpl <file_descriptor_source>
	{
		friend class manager;
		
	public:
		typedef task_t <file_descriptor_source &>	task_type;
	
	private:
		file_descriptor_type			m_fd{-1};
		filter_type						m_filter{};
	
	public:
		file_descriptor_source(
			queue &qq,
			task_type &&tt,
			event_listener_identifier_type const identifier,
			file_descriptor_type const fd,
			filter_type const filter
		):
			source_tpl(qq, std::move(tt), identifier),
			m_fd(fd),
			m_filter(filter)
		{
		}
		
		template <typename ... t_args>
		static std::shared_ptr <file_descriptor_source>
		make_shared(t_args && ... args) { return std::make_shared <file_descriptor_source>(std::forward <t_args>(args)...); }
		
		file_descriptor_type file_descriptor() const { return m_fd; }
		
		// Equivalence class in kqueue.
		file_descriptor_type ident() const { return m_fd; }
		filter_type filter() const { return m_filter; }
	};
	
	typedef file_descriptor_source::task_type file_descriptor_task;
}

#endif
