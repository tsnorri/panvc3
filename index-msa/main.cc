/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <filesystem>
#include <iostream>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/subprocess.hh>
#include <panvc3/dispatch.hh>
#include <panvc3/dispatch/event.hh>
#include <panvc3/msa_index.hh>
#include <string>
#include <string_view>
#include <vector>
#include "cmdline.h"
#include "index_handling.hh"
#include "input_processor.hh"

namespace lb		= libbio;
namespace dispatch	= panvc3::dispatch;
namespace events	= dispatch::events;
namespace fs		= std::filesystem;
namespace mi		= panvc3::msa_indices;


namespace {
	
	struct sigchld_handler final : public events::sigchld_handler
	{
		void child_did_exit_with_nonzero_status(pid_t const pid, int const exit_status, char const *reason) override
		{
			std::cerr << "ERROR: Child process " << pid << " exited with status " << exit_status;
			if (reason)
				std::cerr << " (" << reason << ')';
			std::cerr << '.' << std::endl;
		}
		
		void child_received_signal(pid_t const pid, int const signal_number) override
		{
			std::cerr << "ERROR: Child process " << pid << " received signal " << signal_number << '.' << std::endl;
		}
		
		void finish_handling(bool const did_report_error) override
		{
			if (did_report_error)
				std::exit(EXIT_FAILURE);
		}
	};
	
	
	void install_sigchld_handler(events::manager &mgr)
	{
		static sigchld_handler handler;
		events::install_sigchld_handler(mgr, dispatch::parallel_queue::shared_queue(), handler);
	}
	
	
	void list_index_contents(char const *path)
	{
		panvc3::msa_index msa_index;
		mi::load_msa_index(path, msa_index);
		
		for (auto const &chr_entry : msa_index.chr_entries)
		{
			std::size_t seq_len{};
			auto const &seq_entries(chr_entry.sequence_entries);
			if (!seq_entries.empty())
				seq_len = seq_entries.front().gap_positions.size();

			std::cout << chr_entry.chr_id << " (" << chr_entry.sequence_entries.size() << " sequences of " << seq_len << " characters)\n";
			for (auto const &seq_entry : seq_entries)
			{
				if (seq_entry.gap_positions.size() != seq_len)
				{
					std::cerr << "ERROR: Sequence " << seq_entry.seq_id << " has " << seq_entry.gap_positions.size() << " characters, expected " << seq_len << '.' << std::endl;
					std::exit(EXIT_FAILURE);
				}

				std::cout << '\t' << seq_entry.seq_id << '\n';
			}
		}
	}


	void query_index(char const *path, char const *chr_id)
	{
		auto const msa_index{[&](){
			panvc3::msa_index msa_index;
			mi::load_msa_index(path, msa_index);
			return msa_index;
		}()};

		panvc3::msa_index::chr_entry_cmp const chr_cmp;
		panvc3::msa_index::sequence_entry_cmp const seq_cmp;

		auto const &chr_entries(msa_index.chr_entries);
		auto const chr_rng(std::equal_range(chr_entries.begin(), chr_entries.end(), chr_id, chr_cmp));
		if (chr_rng.first == chr_rng.second)
		{
			std::cerr << "ERROR: No entry for chromosome ‘" <<  chr_id << "’.\n";
			std::exit(EXIT_FAILURE);
		}

		auto const &chr_entry(*chr_rng.first);
		auto const &seq_entries(chr_entry.sequence_entries);
		auto find_seq_entry{[&](auto const &seq_id) -> panvc3::msa_index::sequence_entry const * {
			auto const seq_rng(std::equal_range(seq_entries.begin(), seq_entries.end(), seq_id, seq_cmp));
			if (seq_rng.first == seq_rng.second)
			{
				std::cerr << "No entry for sequence ‘" << seq_id << "’.\n";
				return nullptr;
			}

			return &(*seq_rng.first);
		}};

		std::string buffer;
		auto read_seq_identifier_and_find{[&](char const *msg) -> panvc3::msa_index::sequence_entry const * {
			while (true)
			{
				std::cout << msg << std::flush;
				std::cin >> buffer;
				if (std::cin.eof())
					return nullptr;

				auto const * const entry_ptr(find_seq_entry(buffer));
				if (entry_ptr)
					return entry_ptr;
			}
		}};

		std::uint64_t aln_limit{};
		std::uint64_t pos_limit{};
		auto read_src_seq_identifier_and_find{[&](){
			auto const retval{read_seq_identifier_and_find("Source sequence identifier? ")};
			if (!retval)
				return retval; // Avoid declaring return type.

			aln_limit = retval->gap_positions.size();
			pos_limit = retval->gap_positions_rank0_support(aln_limit);
			return retval;
		}};
		auto read_dst_seq_identifier_and_find{[&](){
			return read_seq_identifier_and_find("Destination sequence identifier? ");
		}};

		auto src_seq_entry{read_src_seq_identifier_and_find()};
		if (!src_seq_entry)
			return;
		auto dst_seq_entry{read_dst_seq_identifier_and_find()};
		if (!dst_seq_entry)
			return;

		while (true)
		{
			std::cout << "[" << src_seq_entry->seq_id << " → " << dst_seq_entry->seq_id << "] Source co-ordinate or ‘s’ or ‘d’ to switch sequence? ([0, " << pos_limit << ")) " << std::flush;

			std::cin >> buffer;
			if (std::cin.eof())
				return;

			if ("s" == buffer)
			{
				src_seq_entry = read_src_seq_identifier_and_find();
				if (!src_seq_entry)
					return;
			}
			else if ("d" == buffer)
			{
				dst_seq_entry = read_dst_seq_identifier_and_find();
				if (!dst_seq_entry)
					return;
			}
			else
			{
				std::uint64_t pos{};
				auto const res(std::from_chars(buffer.data(), buffer.data() + buffer.size(), pos));
				if (std::errc{} != res.ec)
					continue;

				if (! (pos < pos_limit))
					continue;

				std::cout << src_seq_entry->project_position(pos, *dst_seq_entry) << '\n';
			}
		}
	}
	
	
	extern void process(gengetopt_args_info &args_info)
	{
		if (args_info.list_contents_given)
		{
			list_index_contents(args_info.msa_index_input_arg);
			std::exit(EXIT_SUCCESS);
		}
		else if (args_info.query_given)
		{
			query_index(args_info.msa_index_input_arg, args_info.chr_id_arg);
			std::exit(EXIT_SUCCESS);
		}
		else if (args_info.build_index_given)
		{
			if (! (args_info.sequence_inputs_given || args_info.sequences_given))
			{
				std::cerr << "ERROR: Either --sequence-inputs or --sequences has to be specified.\n";
				std::exit(EXIT_FAILURE);
			}
			
			if (args_info.sequence_inputs_given)
			{
				lb::log_time(std::cerr) << "Loading sequence list from " << args_info.sequence_inputs_arg << "…\n";
				static mi::sequence_list_input_processor processor(
					args_info.sequence_inputs_arg,
					args_info.msa_index_input_arg,
					args_info.msa_index_output_arg,
					args_info.pipe_input_arg,
					args_info.output_fasta_flag,
					args_info.fasta_line_width_arg
				);
				
				dispatch::main_queue().async(&processor);
			}
			else if (args_info.sequences_given)
			{
				lb::log_time(std::cerr) << "Loading sequences from A2M file at " << args_info.sequences_arg << "…\n";
				static mi::a2m_input_processor processor(
					args_info.sequences_arg,
					args_info.msa_index_input_arg,
					args_info.msa_index_output_arg,
					args_info.pipe_input_arg,
					args_info.output_fasta_flag,
					args_info.fasta_line_width_arg
				);
				
				dispatch::main_queue().async(&processor);
			}
		}
		else
		{
			std::cerr << "ERROR: Unknown mode.\n";
			std::exit(EXIT_FAILURE);
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	if (args_info.fasta_line_width_arg < 0)
	{
		std::cerr << "ERROR: FASTA line width must be non-negative.\n";
		std::exit(EXIT_FAILURE);
	}
	
	events::manager event_manager;
	event_manager.setup();
	auto manager_thread(event_manager.start_thread_and_run());
	
	install_sigchld_handler(event_manager);
	process(args_info);
	
	dispatch::main_queue().run();
	
	event_manager.stop();
	// manager_thread gets joined here.
	return EXIT_SUCCESS;
}
