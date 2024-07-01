/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_MSA_INDICES_INPUT_HANDLER_HH
#define PANVC3_MSA_INDICES_INPUT_HANDLER_HH

#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/subprocess.hh>
#include <string>
#include <vector>


namespace panvc3::msa_indices {

	struct input_handler
	{
		virtual ~input_handler() {}

		virtual void process_input(
			std::string const &path,
			std::function <void(libbio::file_handle &)> const &cb
		) = 0;
	};


	struct file_input_handler final : public input_handler
	{
		virtual void process_input(
			std::string const &path,
			std::function <void(libbio::file_handle &)> const &cb
		) override
		{
			libbio::file_handle handle(libbio::open_file_for_reading(path));
			cb(handle);
		}
	};


	class subprocess_input_handler final : public input_handler
	{
	public:
		typedef libbio::subprocess <libbio::subprocess_handle_spec::STDOUT> subprocess_type;

	protected:
		std::vector <std::string>	&m_pipe_command;

	public:
		subprocess_input_handler(std::vector <std::string> &pipe_command):
			m_pipe_command(pipe_command)
		{
			m_pipe_command.emplace_back(); // Add an element for the target path. (See process_input() below.)
		}

		virtual void process_input(
			std::string const &path,
			std::function <void(libbio::file_handle &)> const &cb
		) override
		{
			m_pipe_command.back() = path;
			subprocess_type subprocess;
			auto const status(subprocess.open(m_pipe_command));
			if (status)
			{
				cb(subprocess.stdout_handle());
				return;
			}
			
			status.output_status(std::cerr, true);
			std::exit(EXIT_FAILURE);
		}
	};
}

#endif
