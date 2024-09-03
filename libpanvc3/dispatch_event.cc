/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <chrono>
#include <panvc3/dispatch/event.hh>
#include <panvc3/dispatch/task_def.hh>

namespace chrono	= std::chrono;


namespace panvc3::dispatch::events
{
	auto manager_base::check_timers() -> duration_type
	{
		duration_type next_firing_time{timer::DURATION_MAX};
		
		{
			std::lock_guard lock{m_timer_mutex};
		
			if (m_timer_entries.empty())
				return next_firing_time;
			
			do
			{
				auto &next_timer(m_timer_entries.front());
				if (!next_timer.timer->is_enabled())
				{
					std::pop_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
					m_timer_entries.pop_back(); // Safe b.c. the timer uses a shared_ptr_task.
					continue;
				}
				
				auto const duration(next_timer.firing_time - clock_type::now());
				auto const duration_ms(chrono::duration_cast <chrono::milliseconds>(duration).count());
				if (duration_ms <= 0)
				{
					// FIXME: Verify that we could fire while not holding m_timer_mutex?
					// FIXME: The concurrent queue could benefit from sorting the timers by queue and then bulk inserting the operations.
					next_timer.timer->fire();
				
					// If the timer repeats, re-schedule.
					if (next_timer.timer->repeats())
					{
						next_timer.firing_time += next_timer.timer->interval();
						std::make_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
					}
					else
					{
						// m_timers may be empty after popping.
						std::pop_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
						m_timer_entries.pop_back();
					}
				}
				else
				{
					// m_timer_entries not empty.
					next_firing_time = m_timer_entries.front().firing_time - clock_type::now();
					break;
				}
			}
			while (!m_timer_entries.empty());
		}
		
		return next_firing_time;
	}
	
	
	auto manager_base::schedule_timer(
		timer::duration_type const interval,
		bool repeats,
		queue &qq,
		timer::task_type tt
	) -> timer_ptr										// Thread-safe.
	{
		auto retval([&]{
			std::lock_guard lock(m_timer_mutex);
			auto const now(clock_type::now());
		
			// Insert to the correct position.
			auto &entry(m_timer_entries.emplace_back(now + interval, std::make_shared <timer>(qq, std::move(tt), interval, repeats)));
			auto retval(entry.timer);
			std::push_heap(m_timer_entries.begin(), m_timer_entries.end(), heap_cmp_type{});
			return retval;
		}());
	
		trigger_event(event_type::wake_up);
		return retval;
	}
}
