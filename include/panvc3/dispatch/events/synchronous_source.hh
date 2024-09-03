/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_SYNCHRONOUS_SOURCE_HH
#define PANVC3_DISPATCH_EVENTS_SYNCHRONOUS_SOURCE_HH

#include <panvc3/dispatch/events/source.hh>


namespace panvc3::dispatch::events {
	
	// Used with Linux’s special file descriptors in such a way that
	// the event can be handled in the manager’s thread.
	class synchronous_source final : public source_tpl <synchronous_source, false>
	{
	public:
		typedef task_t <synchronous_source &>	task_type;
		
		using source_tpl <synchronous_source, false>::source_tpl;
		
		bool is_read_event_source() const override { return true; }
	};
}

#endif
