/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_COMPRESSED_FASTA_READER_HH
#define PANVC3_COMPRESSED_FASTA_READER_HH

#include <libbio/file_handle.hh>
#include <libbio/utility/compare_strings_transparent.hh>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>


namespace panvc3::detail {
	
	// FIXME: Move to libbio?
	template <typename t_type, typename t_proj, typename t_cmp = std::less <>>
	struct cmp_proj
	{
	private:
		template <typename, typename = void> struct is_transparent_t : public std::false_type {};

		template <typename t_type_>
		struct is_transparent_t <t_type_, std::void_t <typename t_type_::is_transparent>> : public t_type_::is_transparent {};

	public:
		typedef is_transparent_t <t_cmp>	is_transparent;


		bool operator()(t_type const &lhs, t_type const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(proj(lhs), proj(rhs));
		}

		template <typename t_other>
		bool operator()(t_type const &lhs, t_other const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(proj(lhs), rhs);
		}

		template <typename t_other>
		bool operator()(t_other const &lhs, t_type const &rhs) const
		{
			t_cmp cmp;
			t_proj proj;
			return cmp(lhs, proj(rhs));
		}
	};
}


namespace panvc3 {

	struct faidx_entry
	{
		std::string	name;
		std::size_t	length{};
		std::size_t	offset{};
		std::size_t	line_bases{};
		std::size_t	line_width{};

		auto to_ref_tuple() { return std::tie(name, length, offset, line_bases, line_width); }

		struct proj_name
		{
			std::string const &operator()(faidx_entry const &entry) const { return entry.name; }
		};
	};

	typedef detail::cmp_proj <faidx_entry, faidx_entry::proj_name, libbio::compare_strings_transparent>	faidx_entry_cmp;


	class compressed_fasta_reader
	{
	private:
		template <typename>
		struct project_first {};

		template <typename t_lhs, typename t_rhs>
		struct project_first <std::pair <t_lhs, t_rhs>>
		{
			t_lhs const &operator()(std::pair <t_lhs, t_rhs> const &pair) const
			{
				return pair.first;
			}
		};

	public:
		template <typename t_type>
		using pair_t = std::pair <t_type, t_type>;

		typedef std::uint64_t					gzi_offset_type;
		typedef pair_t <gzi_offset_type>		gzi_offset_pair;
		typedef std::vector <gzi_offset_pair>	gzi_table_type;
		typedef std::vector <faidx_entry>		faidx_table_type;
		typedef std::vector <char>				sequence_vector;

		typedef detail::cmp_proj <
			gzi_offset_pair,
			project_first <gzi_offset_pair>
		>										gzi_offset_pair_cmp;

	private:
		libbio::file_handle	m_fasta_handle;
		faidx_table_type	m_sequence_entries;
		gzi_table_type		m_compressed_offsets; // Reversed w.r.t. the GZI file.

	public:
		faidx_table_type const &sequence_entries() { return m_sequence_entries; }
		void read_sequence(faidx_entry const &seq_entry, sequence_vector &dst);
		bool read_sequence(std::string_view const seq_id, sequence_vector &dst);
		void open_path(std::string const &fasta_path);
	};
}

#endif
