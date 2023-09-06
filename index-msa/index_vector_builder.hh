/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_MSA_INDICES_INDEX_VECTOR_BUILDER_HH
#define PANVC3_MSA_INDICES_INDEX_VECTOR_BUILDER_HH

#include <libbio/fasta_reader.hh>
#include <sdsl/int_vector.hpp>
#include <span>
#include <string>
#include <vector>


namespace panvc3::msa_indices {

	class index_vector_builder
	{
	protected:
		sdsl::bit_vector	m_bv{};
		std::size_t			m_pos{};
		std::size_t			m_non_gap_count{};
		std::size_t			m_wrap_amt{};
		bool				m_should_output_seq{};
	
	public:
		index_vector_builder(bool should_output_seq, std::size_t wrap_amt):
			m_wrap_amt(wrap_amt),
			m_should_output_seq(should_output_seq)
		{
		}
		
		sdsl::bit_vector &destination_bit_vector() { return m_bv; }
	
		void begin_sequence(std::string_view const chr_id, std::string_view const seq_id, std::uint64_t const size = 0);
		void update_sequence(std::span <char const> const buffer);
		void end_sequence();
	};
	
	
	void build_index_vector_one_sequence(
		index_vector_builder &builder,
		std::string_view const chr_id,
		std::string_view const seq_id,
		std::vector <char> &buffer,
		libbio::file_handle &read_handle,
		std::size_t const input_size // 0 if not known
	);
	
	
	class index_vector_builder_a2m_input; // Fwd.
	
	
	struct index_vector_builder_a2m_input_delegate
	{
		virtual ~index_vector_builder_a2m_input_delegate() {}

		virtual void index_vector_builder_will_process_sequence(
			index_vector_builder_a2m_input &input,
			index_vector_builder &builder,
			std::string const &chrom,
			std::string const &seq
		) = 0;
		
		virtual void index_vector_builder_did_process_sequence(
			index_vector_builder_a2m_input &input,
			index_vector_builder &builder,
			std::string const &chrom,
			std::string const &seq
		) = 0;
	};
	
	
	class index_vector_builder_a2m_input : libbio::fasta_reader_delegate
	{
	private:
		index_vector_builder					*m_builder{};
		index_vector_builder_a2m_input_delegate	*m_delegate{};
		std::string								m_current_chrom;
		std::string								m_current_seq;
		
	public:
		bool handle_comment_line(libbio::fasta_reader_base &reader, std::string_view const &sv) override { return true; }
		bool handle_identifier(libbio::fasta_reader_base &reader, std::string_view const &sv, std::vector <std::string_view> const &additional_info) override;
		bool handle_sequence_chunk(libbio::fasta_reader_base &reader, std::string_view const &sv, bool has_newline) override;
		bool handle_sequence_end(libbio::fasta_reader_base &reader) override;
		
		void build(index_vector_builder &builder, libbio::file_handle &read_handle, index_vector_builder_a2m_input_delegate &delegate);
	};
}

#endif
