/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_DISPATCH_EVENTS_SIGNAL_MASK_HH
#define PANVC3_DISPATCH_EVENTS_SIGNAL_MASK_HH

#include <signal.h>


namespace panvc3::dispatch::events {
	
	class signal_mask
	{
		sigset_t	m_mask{};

	private:
		bool remove_all_();

	public:
		signal_mask()
		{
			sigemptyset(&m_mask);
		}
	
		~signal_mask() { remove_all_(); }
	
		void add(int sig);
		void remove(int sig);
		void remove_all();
	};
}

#endif
