/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/file_handle.hh>
#include <range/v3/view/enumerate.hpp>
#include <span>
#include <string_view>
#include <sys/stat.h>
#include <vector>
#include "index_vector_builder.hh"

namespace lb	= libbio;
namespace rsv	= ranges::views;


namespace panvc3::msa_indices {

	void index_vector_builder::begin_sequence(std::string_view const chr_id, std::string_view const seq_id, std::uint64_t const size)
	{
		// Output the FASTA header for this sequence.
		if (m_should_output_seq)
			std::cout << '>' << chr_id << '/' << seq_id << '\n';

		m_bv.clear();
		m_bv.resize(size, 0);
		m_pos = 0;
		m_non_gap_count = 0;
	}


	void index_vector_builder::update_sequence(std::span <char const> const buffer)
	{
		for (auto const &[i, cc] : rsv::enumerate(buffer))
		{
			if ('-' == cc)
			{
				// Handle the case where sb.st_size == 0, i.e. reading from a pipe.
				if (m_bv.size() <= m_pos + i)
					m_bv.resize(1 + m_pos + i, 0);

				m_bv[m_pos + i] = 1;
			}
			else if (m_should_output_seq)
			{
				std::cout << cc;
				++m_non_gap_count;

				if (m_wrap_amt && 0 == m_non_gap_count % m_wrap_amt)
					std::cout << '\n';
			}
		}

		m_pos += buffer.size();
	}


	void index_vector_builder::end_sequence()
	{
		m_bv.resize(m_pos, 0);
		if (m_should_output_seq)
		{
			if (0 == m_wrap_amt || 0 != m_non_gap_count % m_wrap_amt)
				std::cout << '\n';
		}
	}


	void build_index_vector_one_sequence(
		index_vector_builder &builder,
		std::string_view const chr_id,
		std::string_view const seq_id,
		std::vector <char> &buffer,
		lb::file_handle &read_handle,
		std::size_t const input_size // 0 if not known
	)
	{
		// Read the input and output non-gapped fasta.

		// Prepare for reading.
		buffer.clear();
		struct stat sb;
		read_handle.stat(sb);
		buffer.resize(sb.st_blksize ?: 4096, 0);

		builder.begin_sequence(chr_id, seq_id, input_size ?: sb.st_size);

		std::size_t bytes_read{};
		while ((bytes_read = read_handle.read(buffer.size(), buffer.data())))
		{
			std::span <char> const buffer_(buffer.data(), bytes_read);
			builder.update_sequence(buffer_);
		}

		// Make the vector long enough in case sb.st_size was zero.
		if (input_size)
			libbio_always_assert_eq(builder.destination_bit_vector().size(), input_size);
	}


	bool index_vector_builder_a2m_input::handle_identifier(lb::fasta_reader_base &reader, std::string_view const &sv, std::vector <std::string_view> const &additional_info)
	{
		// FIXME: Add check that front() below exists.
		m_current_chrom = sv;
		m_current_seq = additional_info.front();
		m_delegate->index_vector_builder_will_process_sequence(*this, *m_builder, m_current_chrom, m_current_seq);
		m_builder->begin_sequence(sv, additional_info.front());
		return true;
	}


	bool index_vector_builder_a2m_input::handle_sequence_chunk(lb::fasta_reader_base &reader, std::string_view const &sv, bool has_newline)
	{
		std::span <char const> const buffer(sv.data(), sv.size());
		m_builder->update_sequence(buffer);
		return true;
	}


	bool index_vector_builder_a2m_input::handle_sequence_end(lb::fasta_reader_base &reader)
	{
		m_builder->end_sequence();
		m_delegate->index_vector_builder_did_process_sequence(*this, *m_builder, m_current_chrom, m_current_seq);
		return true;
	}


	void index_vector_builder_a2m_input::build(index_vector_builder &builder, lb::file_handle &read_handle, index_vector_builder_a2m_input_delegate &delegate)
	{
		m_builder = &builder;
		m_delegate = &delegate;

		lb::fasta_reader reader;
		reader.parse(read_handle, *this);

		m_builder = nullptr;
		m_delegate = nullptr;
	}
}
