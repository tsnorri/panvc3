/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_SEQUENCE_BUFFER_STORE_HH
#define PANVC3_SEQUENCE_BUFFER_STORE_HH

#include <atomic>
#include <vector>


namespace panvc3 {

	class sequence_buffer
	{
		friend class sequence_buffer_store;

	public:
		typedef std::vector <char>	buffer_type;
	
	private:
		buffer_type			m_buffer;
		std::atomic_size_t	m_use_count{};
	
	public:
		buffer_type &get() { return m_buffer; }
		buffer_type const &get() const { return m_buffer; }
	};


	class sequence_buffer_store
	{
	public:
		typedef sequence_buffer::buffer_type	buffer_type;
		typedef std::vector <buffer_type>		buffer_vector;
		typedef std::vector <sequence_buffer>	item_vector;
		constexpr static std::size_t			SPARE_BUFFER_COUNT{4};

	private:
		item_vector		m_items;
		buffer_vector	m_spare_buffers;

	public:
		explicit sequence_buffer_store(std::size_t const size):
			m_items(size)
		{
			m_spare_buffers.reserve(SPARE_BUFFER_COUNT);
		}

		sequence_buffer &acquire_buffer(std::size_t const idx);			// Call from thread 1.
		sequence_buffer const &buffer(std::size_t const idx) const;		// Call from any thread (1, 2, …, n).
		void recycle_buffers();											// Call from thread 1.
		void release_buffer(std::size_t const idx);						// Call from thread 2.
	};
}

#endif
