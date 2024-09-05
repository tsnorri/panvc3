/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_MANAGER_IMPL_HH
#define PANVC3_DISPATCH_EVENTS_MANAGER_IMPL_HH

#	ifdef __linux__
#		include <panvc3/dispatch/events/platform/manager_linux.hh>
namespace panvc3::dispatch::events {
	typedef platform::linux::manager		manager;
}
#	else
#		include <panvc3/dispatch/events/platform/manager_kqueue.hh>
namespace panvc3::dispatch::events {
	typedef platform::kqueue::manager		manager;
	typedef platform::kqueue::signal_mask	signal_mask;
}
#	endif


namespace panvc3::dispatch::events {
	
	void install_sigchld_handler(manager &mgr, queue &qq, sigchld_handler &handler);
}

#endif
