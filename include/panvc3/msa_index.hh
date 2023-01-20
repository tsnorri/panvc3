/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_MSA_INDEX_HH
#define PANVC3_MSA_INDEX_HH

#include <libbio/utility/compare_strings_transparent.hh>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <vector>

namespace panvc3::detail {
	
	template <typename t_rng>
	void fill_gaps(t_rng &&str, sdsl::bit_vector &bv)
	{
		namespace rsv = ranges::views;
		libbio_assert_eq(ranges::size(str), bv.size());
	
		auto rng(
			rsv::enumerate(str)
			| rsv::filter([](auto const &tup){ return '-' == std::get <1>(tup); })
			| rsv::transform([](auto const &tup) -> std::size_t { return std::get <0>(tup); })
		);
	
		for (auto const idx : rng)
			bv[idx] = 1;
	}
}

namespace panvc3 {
	
	typedef std::uint32_t cereal_version_type;
	
	
	struct msa_index
	{
		typedef sdsl::rrr_vector <15>		bit_vector;
		typedef bit_vector::rank_0_type		rank0_support_type;
		typedef bit_vector::select_0_type	select0_support_type;
		
		struct sequence_entry
		{
			std::string				seq_id;
			bit_vector				gap_positions;	// 1 iff. gap
			rank0_support_type		gap_positions_rank0_support;
			select0_support_type	gap_positions_select0_support;
			
			sequence_entry() = default;
			
			sequence_entry(std::string &&seq_id_, sdsl::bit_vector const &gap_positions_):
				seq_id(std::move(seq_id_)),
				gap_positions(gap_positions_),
				gap_positions_rank0_support(&gap_positions),
				gap_positions_select0_support(&gap_positions)
			{
			}
			
			bool operator<(sequence_entry const &other) const { return seq_id < other.seq_id; }
			template <typename t_archive> void load(t_archive &ar, cereal_version_type const version);
			template <typename t_archive> void save(t_archive &ar, cereal_version_type const version) const;
			inline void fix_rank_select_pointers();
		};
		
		struct sequence_entry_cmp
		{
			template <typename t_string> inline bool operator()(sequence_entry const &lhs, t_string const &rhs) const;
			template <typename t_string> inline bool operator()(t_string const &lhs, sequence_entry const &rhs) const;
			inline bool operator()(sequence_entry const &lhs, sequence_entry const &rhs) const;
		};
		
		typedef std::vector <sequence_entry>	sequence_entry_vector;
		
		struct chr_entry
		{
			std::string				chr_id;
			sequence_entry_vector	sequence_entries;
			
			chr_entry() = default;
			
			explicit chr_entry(std::string const &chr_id_):
				chr_id(chr_id_)
			{
			}
			
			bool operator<(chr_entry const &other) const { return chr_id < other.chr_id; }
			template <typename t_archive> void serialize(t_archive &ar, cereal_version_type const version);
		};
		
		struct chr_entry_cmp
		{
			template <typename t_string> inline bool operator()(chr_entry const &lhs, t_string const &rhs) const;
			template <typename t_string> inline bool operator()(t_string const &lhs, chr_entry const &rhs) const;
			inline bool operator()(chr_entry const &lhs, chr_entry const &rhs) const;
		};
		
		typedef std::vector <chr_entry>			chr_entry_vector;
		
		chr_entry_vector	chr_entries;
		
		msa_index() = default;
		
		template <typename t_archive> void serialize(t_archive &ar, cereal_version_type const version);
	};
	
	
	template <typename t_string>
	bool msa_index::sequence_entry_cmp::operator()(sequence_entry const &lhs, t_string const &rhs) const
	{
		libbio::compare_strings_transparent cmp;
		return cmp(lhs.seq_id, rhs);
	}
	
	template <typename t_string>
	bool msa_index::sequence_entry_cmp::operator()(t_string const &lhs, sequence_entry const &rhs) const
	{
		libbio::compare_strings_transparent cmp;
		return cmp(lhs, rhs.seq_id);
	}

	bool msa_index::sequence_entry_cmp::operator()(sequence_entry const &lhs, sequence_entry const &rhs) const
	{
		return lhs.seq_id < rhs.seq_id;
	}
	
	
	template <typename t_string>
	bool msa_index::chr_entry_cmp::operator()(chr_entry const &lhs, t_string const &rhs) const
	{
		libbio::compare_strings_transparent cmp;
		return cmp(lhs.chr_id, rhs);
	}
	
	template <typename t_string>
	bool msa_index::chr_entry_cmp::operator()(t_string const &lhs, chr_entry const &rhs) const
	{
		libbio::compare_strings_transparent cmp;
		return cmp(lhs, rhs.chr_id);
	}

	bool msa_index::chr_entry_cmp::operator()(chr_entry const &lhs, chr_entry const &rhs) const
	{
		return lhs.chr_id < rhs.chr_id;
	}
	
	
	template <typename t_archive>
	void msa_index::sequence_entry::load(t_archive &ar, cereal_version_type const version)
	{
		ar(seq_id);
		ar(gap_positions);
		ar(gap_positions_rank0_support);
		ar(gap_positions_select0_support);
		gap_positions_rank0_support.set_vector(&gap_positions);
		gap_positions_select0_support.set_vector(&gap_positions);
	}
	
	
	template <typename t_archive>
	void msa_index::sequence_entry::save(t_archive &ar, cereal_version_type const version) const
	{
		ar(seq_id);
		ar(gap_positions);
		ar(gap_positions_rank0_support);
		ar(gap_positions_select0_support);
	}
	
	
	template <typename t_archive>
	void msa_index::chr_entry::serialize(t_archive &ar, cereal_version_type const version)
	{
		ar(chr_id);
		ar(sequence_entries);
	}
	
	
	template <typename t_archive>
	void msa_index::serialize(t_archive &ar, cereal_version_type const version)
	{
		ar(chr_entries);
	}
	
	
	void msa_index::sequence_entry::fix_rank_select_pointers()
	{
		gap_positions_rank0_support.set_vector(&gap_positions);
		gap_positions_select0_support.set_vector(&gap_positions);
	}

	
	inline std::ostream &operator<<(std::ostream &os, msa_index::sequence_entry const &entry)
	{
		// FIXME: this is likely quite inefficient.
		for (auto const val : entry.gap_positions)
			os << (+val);
		return os;
	}
	
	
	typedef std::pair <msa_index::sequence_entry, msa_index::sequence_entry> sequence_entry_pair;
	
	// Helper function for writing tests.
	template <typename t_lhs_rng, typename t_rhs_rng>
	void make_sequence_entry_pair(
		t_lhs_rng &&lhs,
		t_rhs_rng &&rhs,
		sequence_entry_pair &dst
	)
	{
		auto const lhs_size(ranges::size(lhs));
		auto const rhs_size(ranges::size(rhs));
		libbio_assert_eq(lhs_size, rhs_size);
		
		sdsl::bit_vector lhsg(lhs_size, 0);
		sdsl::bit_vector rhsg(rhs_size, 0);
		
		detail::fill_gaps(lhs, lhsg);
		detail::fill_gaps(rhs, rhsg);
		
		dst.first = {"", lhsg};
		dst.second = {"", rhsg};
		
		dst.first.fix_rank_select_pointers();
		dst.second.fix_rank_select_pointers();
	}
}

#endif
