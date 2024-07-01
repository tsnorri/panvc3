/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>					// std::partition
#include <iostream>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <panvc3/alignment_input.hh>
#include <panvc3/utility.hh>
#include "cmdline.h"

namespace lb	= libbio;
namespace sam	= libbio::sam;


namespace {

	void process_alignment_group(
		std::vector <sam::record> &alignments,
		std::ostream &os,
		sam::header const &header,
		std::vector <std::size_t> &buffer
	)
	{
		if (alignments.empty())
			return;
		
		// Sort by primary.
		auto it_(std::partition(alignments.begin(), alignments.end(), [](sam::record const &aln_rec){
			return aln_rec.is_primary();
		}));
		
		// Take the SEQ value from some primary alignment.
		auto const seq_([&]() -> sam::record::sequence_type const * {
			for (auto it(alignments.begin()); it != it_; ++it)
			{
				if (!it->seq.empty())
					return &(it->seq);
			}
			
			return nullptr;
		}());
		
		if (!seq_)
		{
			std::cerr << "ERROR: No primary alignment found for read ‘" << alignments.front().qname << "’; skipping.\n";
			return;
		}
		
		// Rewrite SEQ if empty.
		auto const &seq(*seq_);
		for (; it_ != alignments.end(); ++it_)
		{
			if (it_->seq.empty())
				it_->seq = seq;
		}
		
		// Output.
		for (auto const &aln_rec : alignments)
			sam::output_record_in_parsed_order(os, header, aln_rec, buffer);
	}
	
	
	void process_(panvc3::alignment_input &aln_input, std::ostream &os, gengetopt_args_info const &args_info)
	{
		os << aln_input.header;
		
		std::string prev_ref_id;
		std::vector <sam::record> rec_buffer;
		std::vector <std::size_t> buffer;
		
		std::size_t rec_idx{SIZE_MAX};
		aln_input.read_records(
			[
				&aln_input,
				&buffer,
				&os,
				&rec_buffer,
				&rec_idx
			](auto const &aln_rec){
				++rec_idx;
				if (rec_idx && 0 == rec_idx % 10'000'000)
					lb::log_time(std::cerr) << "Processed " << rec_idx << " alignments…\n";
				
				if (rec_buffer.empty())
				{
					rec_buffer.emplace_back(aln_rec);
					return;
				}
				
				auto const &eq_class_id(rec_buffer.front().qname);
				if (aln_rec.qname != eq_class_id)
				{
					process_alignment_group(rec_buffer, os, aln_input.header, buffer);
					rec_buffer.clear();
				}
				
				rec_buffer.emplace_back(aln_rec);
			}
		);
		
		if (!rec_buffer.empty())
			process_alignment_group(rec_buffer, os, aln_input.header, buffer);
	}
	
	
	void append_program_info(sam::header &output_header, int const argc, char const * const * const argv)
	{
		panvc3::append_sam_program_info(
			"panvc3.process-alignments.",
			"PanVC 3 process_alignments",
			argc,
			argv,
			CMDLINE_PARSER_VERSION,
			output_header.programs
		);
	}
	
	
	void process(gengetopt_args_info const &args_info, int const argc, char const * const * const argv)
	{
		// Open the SAM input and output.
		auto aln_input(panvc3::alignment_input::open_path_or_stdin(args_info.alignments_arg));
		
		{
			auto aln_output_fh{[&] -> lb::file_handle {
				if (args_info.output_path_arg)
					return lb::file_handle(lb::open_file_for_writing(args_info.output_path_arg, lb::writing_open_mode::CREATE));
				else
					return lb::file_handle(STDOUT_FILENO, false);
			}()};
			
			lb::file_ostream os;
			lb::open_stream_with_file_handle(os, aln_output_fh);

			auto &header(aln_input.header);
			if (sam::sort_order_type::queryname != header.sort_order)
			{
				std::cerr << "ERROR: Expected input to be sorted by QNAME, got " << header.sort_order << " instead.\n";
				std::exit(EXIT_FAILURE);
			}
			
			append_program_info(header, argc, argv); // Adding our program info does not affect parsing.
			process_(aln_input, os, args_info);
		}
		
		lb::log_time(std::cerr) << "Done.\n";
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	//lb::setup_allocated_memory_logging();
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.

	if (args_info.print_invocation_given)
	{
		std::cerr << "Invocation:";
		for (int i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << '\n';
	}

	if (args_info.print_pid_given)
		std::cerr << "PID: " << getpid() << '\n';
	
	process(args_info, argc, argv);
	
	return EXIT_SUCCESS;
}
