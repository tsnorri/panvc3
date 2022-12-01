/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_SPSC_QUEUE_HH
#define PANVC3_SPSC_QUEUE_HH

#include <array>
#include <dispatch/dispatch.h>
#include <libbio/dispatch/dispatch_ptr.hh>
#include <numeric>							// std::iota
#include <panvc3/utility.hh>				// is_power_of_2()


namespace panvc3 {
	
	template <typename t_value, std::uint16_t t_size>
	class spsc_queue
	{
		static_assert(is_power_of_2(t_size));
		
	public:
		typedef t_value			value_type;
		typedef std::uint16_t	size_type;
		
	protected:
		struct index
		{
			alignas(std::hardware_destructive_interference_size) size_type value{}; // Not sure if the alignment is needed.
			
			// Assignment operator.
			template <typename t_value_>
			index &operator=(t_value_ const value_) { value = value_; return *this; }
		};
		
		typedef std::array <value_type, t_size>				value_array;
		typedef std::array <index, t_size>					index_array;
		typedef libbio::dispatch_ptr <dispatch_semaphore_t>	dispatch_semaphore_ptr;
		
	protected:
		value_array				m_values;
		index_array				m_indices;
		dispatch_semaphore_ptr	m_semaphore_ptr{};
		size_type				m_read_idx{};		// Used by thread 1.
		size_type				m_write_idx{};		// Used by thread 2.
		
	public:
		spsc_queue():
			m_semaphore_ptr(dispatch_semaphore_create(t_size))
		{
			std::iota(m_indices.begin(), m_indices.end(), 0);
		}
		
		value_array &values() { return m_values; }
		value_type &operator[](size_type const idx) { return m_values[idx]; }
		
		size_type pop_index();
		void push(value_type &value);
	};
	
	
	template <typename t_value, std::uint16_t t_size>
	auto spsc_queue <t_value, t_size>::pop_index() -> size_type
	{
		// Called from thread 1.
		dispatch_semaphore_wait(*m_semaphore_ptr, DISPATCH_TIME_FOREVER);
		
		auto const val_idx(m_indices[m_read_idx].value);
		++m_read_idx;
		m_read_idx %= t_size;
		
		return val_idx;
	}
	
	
	template <typename t_value, std::uint16_t t_size>
	void spsc_queue <t_value, t_size>::push(value_type &val)
	{
		// Called from thread 2.
		auto const val_idx(std::distance(&m_values.front(), &val));
		m_indices[m_write_idx].value = val_idx;
		++m_write_idx;
		m_write_idx %= t_size;
		
		dispatch_semaphore_signal(*m_semaphore_ptr);
	}
}

#endif
