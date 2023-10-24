/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_FWD_HH
#define PANVC3_DISPATCH_FWD_HH

namespace panvc3::dispatch {
	
	// Fwd.
	class thread_pool;
	struct queue;
	class serial_queue;
	class parallel_queue;
	
	
	typedef std::uintptr_t	event_listener_identifier_type;	// Type from struct kevent.
	constexpr static inline event_listener_identifier_type const EVENT_LISTENER_IDENTIFIER_MAX{UINTPTR_MAX};
}


namespace panvc3::dispatch::detail {
	struct serial_queue_executor_callable;
}

#endif
