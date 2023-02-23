/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/endian.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iterator>
#include <libbio/assert.hh>
#include <libbio/file_handling.hh>
#include <libbio/generic_parser.hh>
#include <panvc3/compressed_fasta_reader.hh>

namespace ios		= boost::iostreams;
namespace endian	= boost::endian;
namespace lb		= libbio;
namespace lbp		= libbio::parsing;


namespace {

	template <typename t_type>
	inline t_type read_le_and_advance(std::istreambuf_iterator <t_type> &it)
	{
		auto const val(*it);
		++it;
		return endian::little_to_native(val);
	}


	template <typename t_device, typename t_integer>
	class converting_source : public ios::device <ios::input, t_integer>
	{
	public:
		struct incomplete_data : public std::exception {};
	
	public:
		t_device		device{};

	private:
		t_integer		m_buffer{};
		std::uint8_t	m_buffered_bytes{};

	public:
		converting_source(t_device &&device_):
			device(std::move(device_))
		{
		}

		template <typename ... t_args>
		converting_source(t_args && ... args):
			device(std::forward <t_args>(args)...)
		{
		}

		std::streamsize read(t_integer * const dst, std::streamsize const count)
		{
			// Since the function prototype is required to have a t_integer * parameter, there should be no alignment mismatches.
			*dst = m_buffered_bytes;
			auto * const dst_(reinterpret_cast <char *>(dst));
			auto const res(ios::read(device, dst_ + m_buffered_bytes, sizeof(t_integer) * count - m_buffered_bytes));

			// Check for EOF.
			if (-1 == res)
			{
				if (m_buffered_bytes)
					throw incomplete_data{};

				return -1;
			}

			auto const total_bytes(res + m_buffered_bytes);
			auto const retval(total_bytes / sizeof(t_integer));

			// Copy the remaining bytes to m_buffer.
			m_buffered_bytes = total_bytes % sizeof(t_integer);
			if (m_buffered_bytes)
			{
				auto const mask{~t_integer{} << (CHAR_BIT * (sizeof(t_integer) - m_buffered_bytes))};
				m_buffer = *(dst + retval * sizeof(t_integer)) & mask;
			}

			return retval;
		}
	};
}


namespace panvc3 {

	void compressed_fasta_reader::read_sequence(faidx_entry const &seq_entry, sequence_vector &dst)
	{
		auto const block_offsets([&](){
			auto const begin(m_compressed_offsets.begin());
			auto const end(m_compressed_offsets.end());
			auto const it(std::upper_bound(
				begin,
				end,
				seq_entry.offset,
				gzi_offset_pair_cmp
				{}
			));

			libbio_always_assert_neq(begin, it);
			return *(it - 1);
		}());
			
		// Set up a chain for reading compressed input.
		ios::file_descriptor_source source(m_fasta_handle.get(), ios::never_close_handle);
		ios::filtering_streambuf <ios::input> in;
		in.push(ios::gzip_decompressor());
		in.push(source);

		auto remaining_length(seq_entry.length);
		std::size_t read_offset{seq_entry.offset - block_offsets.first};

		dst.clear();
		dst.resize(remaining_length);
		m_fasta_handle.seek(block_offsets.second);
		auto *dst_ptr(dst.data()); // Contiguous since std::vector.

		// Read one line at a time so that the newlines can be skipped.
		while (remaining_length)
		{
			auto const read_amt(std::min(remaining_length, seq_entry.line_bases));

			// I couldn’t find any other way to ignore characters.
			// ios::seek() does not work with non-random-access devices
			// even if the seeking direction is forward and the position
			// is given relative to the current position.
			for (std::size_t i(0); i < read_offset; ++i)
				ios::get(in);

			auto const read_amt_(ios::read(in, dst_ptr, read_amt));

			dst_ptr += read_amt_;
			remaining_length -= read_amt_;
			read_offset = 1; // Skip the newline.
		}
	}


	bool compressed_fasta_reader::read_sequence(std::string_view const seq_id, sequence_vector &dst)
	{
		faidx_entry_cmp cmp;
		auto const end(m_sequence_entries.end());
		auto const it(std::lower_bound(m_sequence_entries.begin(), end, seq_id, cmp));

		if (end == it)
			return false;

		if (it->name != seq_id)
			return false;

		read_sequence(*it, dst);
		return true;
	}


	void compressed_fasta_reader::open_path(std::string const &fasta_path)
	{
		{
			auto const fd(lb::open_file_for_reading(fasta_path));
			m_fasta_handle = lb::file_handle(fd);
		}

		{
			// Open the GZI file.
			auto const gzi_path(fasta_path + ".gzi");
			auto const fd(lb::open_file_for_reading(gzi_path));
			lb::file_handle const gzi_handle(fd);
			ios::stream <converting_source <ios::file_descriptor_source, std::uint64_t>> stream;
			stream.open(gzi_handle.get(), ios::never_close_handle);
			std::istreambuf_iterator <std::uint64_t> it(stream);
			std::istreambuf_iterator <std::uint64_t> const end;

			// Read the count.
			libbio_always_assert_neq(it, end);
			auto const count(read_le_and_advance(it));

			m_compressed_offsets.reserve(1 + count);

			// Read the offsets.
			bool is_first{true};
			bool did_add{false};
			while (it != end)
			{
				auto const compressed_offset(read_le_and_advance(it));
				libbio_always_assert_neq(it, end);
				auto const uncompressed_offset(read_le_and_advance(it));

				// The GZI file may not have an entry for the initial block.
				if (is_first && 0 != compressed_offset)
				{
					did_add = true;
					m_compressed_offsets.emplace_back(0, 0);
				}
				is_first = false;

				// Reverse.
				m_compressed_offsets.emplace_back(uncompressed_offset, compressed_offset);
			}

			libbio_always_assert(std::is_sorted(m_compressed_offsets.begin(), m_compressed_offsets.end(), gzi_offset_pair_cmp{}));
			libbio_always_assert_eq(count + did_add, m_compressed_offsets.size());
		}

		{
			// FIXME: come up with a way to convert the tuple representation of faidx_entry to a parser specification.
			typedef lbp::parser <
				lbp::traits::delimited <lbp::delimiter <'\t'>, lbp::delimiter <'\n'>>,
				lbp::fields::text <>,
				lbp::fields::numeric <std::size_t>,
				lbp::fields::numeric <std::size_t>,
				lbp::fields::numeric <std::size_t>,
				lbp::fields::numeric <std::size_t>
			> faidx_parser_type;
			typedef faidx_parser_type::record_type faidx_record_type;

			// Open the FASTA index.
			auto const fai_path(fasta_path + ".fai");
			lb::file_istream stream;
			lb::open_file_for_reading(fai_path, stream);
			std::istreambuf_iterator <char> it(stream);
			std::istreambuf_iterator <char> const end;

			// Parse.
			faidx_parser_type parser;
			faidx_record_type buffer;
			while (parser.parse(it, end, buffer))
			{
				auto &rec(m_sequence_entries.emplace_back());
				rec.to_ref_tuple() = buffer; // Copy since we would like to re-use the buffer anyway.
			}

			// Sort.
			std::sort(m_sequence_entries.begin(), m_sequence_entries.end(), faidx_entry_cmp{});
		}
	}
}
