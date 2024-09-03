/*
 * Copyright (c) 2023-2024 Tuukka Norri
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
	public:
		typedef task_t <file_descriptor_source &>	task_type;
		
		enum class type
		{
			read_source,
			write_source
		};
	
	private:
		file_descriptor_type			m_fd{-1};
		type							m_source_type;
	
	public:
		file_descriptor_source(
			queue &qq,
			task_type &&tt,
			file_descriptor_type const fd,
			type const st
		):
			source_tpl(qq, std::move(tt)),
			m_fd(fd),
			m_source_type(st)
		{
		}
		
		file_descriptor_type file_descriptor() const { return m_fd; }
		type file_descriptor_source_type() const { return m_source_type; }
		bool is_read_event_source() const override { return type::read_source == m_source_type; }
		bool is_write_event_source() const override { return type::write_source == m_source_type; }
	};
	
	typedef file_descriptor_source::task_type file_descriptor_task;
}


namespace panvc3::dispatch::detail {
	
	template <>
	struct member_callable_target <events::file_descriptor_source>
	{
		typedef events::source_tpl <events::file_descriptor_source>	type;
	};
}

#endif
