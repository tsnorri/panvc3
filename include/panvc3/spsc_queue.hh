/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_SPSC_QUEUE_HH
#define PANVC3_SPSC_QUEUE_HH

#include <atomic>
#include <libbio/bits.hh>					// libbio::bits::gte_power_of_2
#include <limits>							// std::numeric_limits
#include <numeric>							// std::iota
#include <vector>

#if defined(PANVC3_USE_GCD) && PANVC3_USE_GCD
#	include <dispatch/dispatch.h>
#	include <libbio/dispatch/dispatch_ptr.hh>
#else
#	include <semaphore>
#endif


namespace panvc3 {
	
#if defined(PANVC3_USE_GCD) && PANVC3_USE_GCD
	class dispatch_semaphore_gcd
	{
	public:
		typedef libbio::dispatch_ptr <dispatch_semaphore_t>	dispatch_semaphore_ptr;
		
	private:
		dispatch_semaphore_ptr	m_semaphore_ptr{};
		
	public:
		explicit dispatch_semaphore_gcd(std::uint16_t const size):
			m_semaphore_ptr(dispatch_semaphore_create(size))
		{
		}
		
		void acquire() { dispatch_semaphore_wait(*m_semaphore_ptr, DISPATCH_TIME_FOREVER); }
		void release() { dispatch_semaphore_signal(*m_semaphore_ptr); }
	};
#endif
	
	
	template <typename t_value>
	class spsc_queue
	{
	public:
		typedef t_value			value_type;
		typedef std::uint16_t	size_type;

		constexpr static size_type const MAX_SIZE{std::numeric_limits <size_type>::max()};
		
	protected:
		class alignas(std::hardware_destructive_interference_size) index
		{
		private:
			std::atomic <size_type> m_value{};

		public:
			// Assignment operator.
			template <typename t_value_>
			index &operator=(t_value_ const value_) { m_value.store(value_, std::memory_order_relaxed); return *this; }

			/* implicit */ operator size_type() const { return m_value.load(std::memory_order_relaxed); }
		};
		
		typedef std::vector <value_type>					value_array;
		typedef std::vector <index>							index_array;
#if defined(PANVC3_USE_GCD) && PANVC3_USE_GCD		
		typedef dispatch_semaphore_gcd						semaphore_type;
#else
		typedef std::counting_semaphore <MAX_SIZE>			semaphore_type;
#endif
		
	protected:
		value_array				m_values;
		index_array				m_indices;
		semaphore_type			m_semaphore;
		size_type				m_index_mask{};
		size_type				m_read_idx{};		// Used by thread 1.
		size_type				m_write_idx{};		// Used by thread 2.
		
	private:
		inline std::size_t queue_size(size_type const size_);

	public:
		spsc_queue(size_type const size_, size_type const queue_size_):
			m_values(queue_size_),
			m_indices(queue_size_),
			m_semaphore(queue_size_),
			m_index_mask(queue_size_ - 1)
		{
			libbio_assert_lt(0, size());
			std::iota(m_indices.begin(), m_indices.end(), 0);
		}
		
		explicit spsc_queue(size_type const size_):
			spsc_queue(size_, queue_size(size_))
		{
		}

		constexpr size_type size() const { return m_values.size(); }
		
		value_array &values() { return m_values; }
		value_type &operator[](size_type const idx) { libbio_assert_lt(idx, m_values.size()); return m_values[idx]; }
		
		size_type pop_index();
		value_type &pop() { return m_values[pop_index()]; }
		void push(value_type &value);
	};
	
	
	template <typename t_value>
	inline std::size_t spsc_queue <t_value>::queue_size(size_type const size_)
	{
		auto const power(libbio::bits::gte_power_of_2(std::size_t(size_)));
		if (!power)
			throw std::range_error("Unable to construct a queue of the given size");
		
		auto const retval(*power);
		if (std::numeric_limits <size_type>::max() < retval - 1)
			throw std::range_error("Unable to construct a queue of the given size");
		
		return retval;
	}
	
	
	template <typename t_value>
	auto spsc_queue <t_value>::pop_index() -> size_type
	{
		// Called from thread 1.
		m_semaphore.acquire();

		libbio_assert_lt(m_read_idx, m_indices.size());
		size_type const val_idx(m_indices[m_read_idx]);
		++m_read_idx;
		m_read_idx &= m_index_mask;
		
		return val_idx;
	}
	
	
	template <typename t_value>
	void spsc_queue <t_value>::push(value_type &val)
	{
		// Called from thread 2.
		auto const val_idx(std::distance(&m_values.front(), &val));
		libbio_assert_lt(m_write_idx, m_indices.size());
		m_indices[m_write_idx] = val_idx;
		++m_write_idx;
		m_write_idx &= m_index_mask;

		m_semaphore.release();
	}
}

#endif
