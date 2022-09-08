/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_MSA_INDEX_HH
#define PANVC3_MSA_INDEX_HH

#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <vector>

namespace panvc3 {
	
	struct msa_index
	{
		typedef sdsl::rrr_vector <15>		bit_vector;
		typedef bit_vector::rank_0_type		rank0_support_type;
		typedef bit_vector::select_0_type	select0_support_type;
		
		struct sequence_entry
		{
			std::string				seq_id;
			bit_vector				gap_positions;
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
			template <typename t_archive> void load(t_archive &ar, std::uint32_t const version);
			template <typename t_archive> void save(t_archive &ar, std::uint32_t const version) const;
		};
		
		struct sequence_entry_cmp
		{
			bool operator()(sequence_entry const &lhs, std::string const &rhs) const { return lhs.seq_id < rhs; }
			bool operator()(std::string const &lhs, sequence_entry const &rhs) const { return lhs < rhs.seq_id; }
		};
		
		typedef std::vector <sequence_entry>	sequence_entry_vector;
		
		struct chr_entry
		{
			std::string				chr_id;
			sequence_entry_vector	sequence_entries;
			
			explicit chr_entry(std::string const &chr_id_):
				chr_id(chr_id_)
			{
			}
			
			bool operator<(chr_entry const &other) const { return chr_id < other.chr_id; }
			template <typename t_archive> void serialize(t_archive &ar, std::uint32_t const version);
		};
		
		struct chr_entry_cmp
		{
			bool operator()(chr_entry const &lhs, std::string const &rhs) const { return lhs.chr_id < rhs; }
			bool operator()(std::string const &lhs, chr_entry const &rhs) const { return lhs < rhs.chr_id; }
		};
		
		typedef std::vector <chr_entry>			chr_entry_vector;
		
		chr_entry_vector	chr_entries;
		
		msa_index() = default;
		template <typename t_archive> void serialize(t_archive &ar, std::uint32_t const version);
	};
	
	
	template <typename t_archive>
	void msa_index::sequence_entry::load(t_archive &ar, std::uint32_t const version)
	{
		ar(seq_id);
		ar(gap_positions);
		ar(gap_positions_rank0_support);
		ar(gap_positions_select0_support);
		gap_positions_rank0_support.set_vector(&gap_positions);
		gap_positions_select0_support.set_vector(&gap_positions);
	}
	
	
	template <typename t_archive>
	void msa_index::sequence_entry::save(t_archive &ar, std::uint32_t const version) const
	{
		ar(seq_id);
		ar(gap_positions);
		ar(gap_positions_rank0_support);
		ar(gap_positions_select0_support);
	}
	
	
	template <typename t_archive>
	void msa_index::chr_entry::serialize(t_archive &ar, std::uint32_t const version)
	{
		ar(chr_id);
		ar(sequence_entries);
	}
	
	
	template <typename t_archive>
	void msa_index::serialize(t_archive &ar, std::uint32_t const version)
	{
		ar(chr_entries);
	}
}

#endif
