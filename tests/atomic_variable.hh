/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_TEST_ATOMIC_VARIABLE_HH
#define PANVC3_TEST_ATOMIC_VARIABLE_HH

#include <chrono>
#include <condition_variable>
#include <mutex>


namespace panvc3::tests {
	
	template <typename t_value, t_value t_initial = t_value{}>
	class atomic_variable
	{
	public:
		typedef std::chrono::steady_clock	clock_type;
		
	private:
		t_value								m_value{t_initial};
		std::condition_variable				m_cv{};
		std::mutex							m_mutex{};
		
	public:
		template <typename t_value_>
		inline void assign(t_value_ &&value);
		
		std::lock_guard <std::mutex> lock() { return std::lock_guard{m_mutex}; }
		inline std::unique_lock <std::mutex> wait_and_lock(clock_type::duration const dur = std::chrono::seconds(1));
		t_value &value() { return m_value; }
		t_value const &value() const { return m_value; }
	};
	
	
	typedef atomic_variable <bool>			atomic_bool;
	typedef atomic_variable <std::uint32_t>	atomic_uint32_t;
	
	
	template <typename t_value, t_value t_initial>
	template <typename t_value_>
	void atomic_variable <t_value, t_initial>::assign(t_value_ &&value)
	{
		{
			std::lock_guard const lock{m_mutex};
			m_value = std::forward <t_value_>(value);
		}
		
		m_cv.notify_one();
	}
	
	
	template <typename t_value, t_value t_initial>
	std::unique_lock <std::mutex> atomic_variable <t_value, t_initial>::wait_and_lock(clock_type::duration const dur)
	{
		std::unique_lock lock{m_mutex};
		m_cv.wait_for(lock, dur, [this]{ return m_value != t_initial; });
		return lock;
	}
}

#endif
