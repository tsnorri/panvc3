/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_MSA_INDICES_INPUT_PROCESSOR_HH
#define PANVC3_MSA_INDICES_INPUT_PROCESSOR_HH

#include <cstddef>
#include <libbio/dispatch.hh>
#include <libbio/subprocess.hh>
#include <mutex>
#include <panvc3/msa_index.hh>
#include <string>
#include <vector>
#include "index_vector_builder.hh"
#include "input_handler.hh"


namespace panvc3::msa_indices {

	class input_processor
	{
	protected:
		index_vector_builder		m_index_vector_builder;
		std::string					m_input_path;
		std::string					m_msa_index_input_path;
		std::string					m_msa_index_output_path;
		std::vector <std::string>	m_pipe_command;

	public:
		input_processor(
			char const *input_path,
			char const *msa_index_input_path,
			char const *msa_index_output_path,
			char const *pipe_input_command,
			bool const should_output_fasta,
			std::size_t const fasta_line_width
		):
			m_index_vector_builder(should_output_fasta, fasta_line_width),
			m_input_path(input_path),
			m_msa_index_input_path(msa_index_input_path ?: ""),
			m_msa_index_output_path(msa_index_output_path),
			m_pipe_command(libbio::parse_command_arguments(pipe_input_command))
		{
		}

		virtual void process(input_handler &handler) = 0;

		void operator()()
		{
			if (m_pipe_command.empty())
			{
				file_input_handler handler;
				process(handler);
			}
			else
			{
				subprocess_input_handler handler(m_pipe_command);
				process(handler);
			}
		}
	};


	class sequence_list_input_processor final : public input_processor
	{
	public:
		using input_processor::input_processor;
		void operator()() { input_processor::operator()(); } // Needed by Clang (16.0.6) b.c. “non-type template argument of pointer-to-member type that refers to member of a different class is not supported yet”.
		void process(input_handler &handler) override;
	};


	class a2m_input_processor final : public input_processor, public index_vector_builder_a2m_input_delegate
	{
	private:
		std::mutex				m_msa_index_mutex{};
		panvc3::msa_index		m_msa_index;
		libbio::dispatch::group	m_main_group;

	public:
		using input_processor::input_processor;
		void operator()() { input_processor::operator()(); } // See above.
		void process(input_handler &handler) override;

		void index_vector_builder_will_process_sequence(
			index_vector_builder_a2m_input &input,
			index_vector_builder &builder,
			std::string const &chrom,
			std::string const &seq
		) override;

		void index_vector_builder_did_process_sequence(
			index_vector_builder_a2m_input &input,
			index_vector_builder &builder,
			std::string const &chrom,
			std::string const &seq
		) override;
	};
}

#endif
