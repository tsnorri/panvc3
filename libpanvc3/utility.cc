/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <charconv>				// std::from_chars
#include <panvc3/utility.hh>

namespace dispatch	= libbio::dispatch;
namespace sam		= libbio::sam;


namespace panvc3 {
	
	std::uint32_t check_sam_program_id_index(std::string_view const id_prefix, std::string_view const id)
	{
		if (id.starts_with(id_prefix))
		{
			auto const pos(id_prefix.size());
			auto const tail(id.substr(pos));
			std::uint32_t idx{};
			auto const res(std::from_chars(tail.data(), tail.data() + tail.size(), idx, 10));
			if (std::errc{} == res.ec)
				return idx;
		}
		
		return 0;
	}
	
	
	std::string command_line_call(
		int const argc,
		char const * const * const argv
	)
	{
		// CLI call
		std::string cmd;
		for (int i(0); i < argc; ++i)
		{
			if (0 != i)
				cmd.append(" ");
			cmd.append(argv[i]);
		}
		
		return cmd;
	}
	

	void append_sam_program_info(
		std::string_view const id_prefix,
		std::string_view const name,
		int const argc,
		char const * const * const argv,
		char const * const version,
		sam::program_entry_vector &dst
	)
	{
		// Identifier
		std::uint32_t idx{};
		for (auto const &pginfo : dst)
			idx = std::max(idx, check_sam_program_id_index(id_prefix, pginfo.id));
		
		++idx;
		std::string id(id_prefix);
		id += std::to_string(idx);

		// Build the record.
		sam::program_entry rec{
			.id{std::move(id)},
			.name{std::string(name)},
			.command_line{command_line_call(argc, argv)},
			.version{version}
		};

		if (!dst.empty())
			rec.prev_id = dst.back().id;

		dst.emplace_back(std::move(rec));
	}
	
	
	void prepare_thread_pool_with_args(dispatch::thread_pool &thread_pool, long threads) // Type from gengetopt
	{
		thread_pool.set_min_workers(2); // Needed for the BAM reader.
	
		if (threads < 0)
		{
			std::cerr << "ERROR: Number of threads must be non-negative.\n";
			std::exit(EXIT_FAILURE);
		}
		else if (0 != threads)
		{
			std::size_t thread_pool_max_size{};
			if (threads < 3)
			{
				std::cerr << "INFO: Using three threads.\n";
				thread_pool_max_size = 2;
			}
			else
			{
				thread_pool_max_size = threads;
			}
		
			thread_pool.set_max_workers(thread_pool_max_size);
		}
	}
}
