/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/sort/sort.hpp>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/utility/misc.hh>
#include <panvc3/alignment_input.hh>
#include <panvc3/utility.hh>
#include <range/v3/view/enumerate.hpp>
#include <string>
#include "cmdline.h"

namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace sam	= libbio::sam;


namespace {
	
	typedef std::uint8_t						mapping_quality_type;
	
	
	struct match_count
	{
		std::size_t matches{};
		std::size_t mismatches{};
	};
	
	
	void output_best_mate(
		sam::reference_id_type const mate_ref_id, 
		sam::position_type const mate_pos,
		std::vector <sam::record> &alignments,
		std::ostream &os,
		sam::header const &header,
		std::vector <std::size_t> &buffer
	)
	{
		if (alignments.empty())
			return;
		
		// Determine again the best mapping quality.
		mapping_quality_type best_mapq{};
		for (auto const &aln_rec : alignments)
		{
			if (aln_rec.rname_id != mate_ref_id)
				continue;
			
			if (aln_rec.pos != mate_pos)
				continue;
			
			auto const mapq(aln_rec.mapq);
			if (255 == mapq)
				continue;
			best_mapq = std::max(mapq, best_mapq);
		}
		
		for (auto &aln_rec : alignments)
		{
			if (aln_rec.mapq == best_mapq && aln_rec.rname_id == mate_ref_id && aln_rec.pos == mate_pos)
			{
				sam::output_record_in_parsed_order(os, header, aln_rec, buffer);
				os << '\n';
				return;
			}
		}
	}
	

	void process_alignment_group(
		std::vector <sam::record> &alignments,
		std::ostream &os,
		sam::header const &header,
		std::vector <std::size_t> &buffer
	)
	{
		if (alignments.empty())
			return;
		
		// Determine the best mapping quality.
		mapping_quality_type best_mapq{};
		for (auto const &aln_rec : alignments)
		{
			auto const mapq(aln_rec.mapq);
			if (255 == mapq)
				continue;
			
			best_mapq = std::max(mapq, best_mapq);
		}
		
		// Find again.
		for (auto &aln_rec : alignments)
		{
			if (aln_rec.mapq == best_mapq)
			{
				sam::output_record_in_parsed_order(os, header, aln_rec, buffer);
				os << '\n';
				
				auto const mate_ref_id(aln_rec.rnext_id);
				auto const mate_pos(aln_rec.pnext);
				if (sam::INVALID_REFERENCE_ID == mate_ref_id || sam::INVALID_POSITION == mate_pos)
					return;
				
				output_best_mate(mate_ref_id, mate_pos, alignments, os, header, buffer);
				return;
			}
		}
		
		// Output the first by default.
		auto &aln_rec(alignments.front());
		sam::output_record_in_parsed_order(os, header, aln_rec, buffer);
		os << '\n';
		
		auto const mate_ref_id(aln_rec.rnext_id);
		auto const mate_pos(aln_rec.pnext);
		if (sam::INVALID_REFERENCE_ID == mate_ref_id || sam::INVALID_POSITION == mate_pos)
			return;
		
		output_best_mate(mate_ref_id, mate_pos, alignments, os, header, buffer);
	}
	
	
	match_count process_(panvc3::alignment_input &aln_input, std::ostream &os, gengetopt_args_info const &args_info)
	{
		os << aln_input.header;
		
		bool const should_subset_by_read_id(args_info.read_id_flag);
		bool const should_subset_by_best_mapq(args_info.best_mapq_flag);
		match_count mc{};

		std::deque <std::string> read_names;
		
		// Read the read names.
		if (should_subset_by_read_id)
		{
			lb::log_time(std::cerr) << "Reading the read names…\n";
			std::string buffer;
			std::size_t lineno{};
			while (std::getline(std::cin, buffer))
			{
				++lineno;
				read_names.emplace_back(buffer);
				if (0 == lineno % 100000000)
					lb::log_time(std::cerr) << "Handled " << lineno << " lines…\n";
			}
			
			lb::log_time(std::cerr) << "Sorting…\n";
			boost::sort::block_indirect_sort(read_names.begin(), read_names.end());
			//std::sort(read_names.begin(), read_names.end());
		}
		
		// Process the records.
		lb::log_time(std::cerr) << "Processing the alignment records…\n";
		auto const &ref_seqs(aln_input.header.reference_sequences);
		std::string prev_ref_id;
		std::vector <sam::record> rec_buffer;
		std::vector <std::size_t> buffer;
		
		auto const expected_ref_id([&](){
			if (!args_info.chr_arg)
				return sam::INVALID_REFERENCE_ID;
			
			auto const retval(aln_input.header.find_reference(args_info.chr_arg));
			if (sam::INVALID_REFERENCE_ID == retval)
			{
				std::cerr << "ERROR: Reference ID " << args_info.chr_arg << " not found in SAM header.\n";
				std::exit(EXIT_FAILURE);
			}
			
			return retval;
		}());
		
		std::size_t rec_idx{};
		aln_input.read_records(
			[
				&aln_input,
				&buffer,
				&mc,
				&os,
				&read_names,
				&rec_buffer,
				&rec_idx,
				expected_ref_id,
				should_subset_by_best_mapq,
				should_subset_by_read_id
			](auto const &aln_rec){
				if (rec_idx && 0 == rec_idx % 10'000'000)
					lb::log_time(std::cerr) << "Processed " << rec_idx << " alignments…\n";
				++rec_idx;
				
				if (sam::INVALID_REFERENCE_ID != expected_ref_id)
				{
					auto const ref_id(aln_rec.rname_id);
					if (sam::INVALID_REFERENCE_ID == ref_id)
					{
						++mc.mismatches;
						return;
					}
					
					if (ref_id != expected_ref_id)
					{
						++mc.mismatches;
						return;
					}
				}
				
				auto const &qname(aln_rec.qname);
				
				if (should_subset_by_read_id && !std::binary_search(read_names.begin(), read_names.end(), qname))
				{
					++mc.mismatches;
					return;
				}
				
				++mc.matches;
				if (should_subset_by_best_mapq)
				{
					if (rec_buffer.empty())
					{
						rec_buffer.emplace_back(aln_rec);
						return;
					}
					
					auto const &eq_class_id(rec_buffer.front().qname);
				
					if (qname != eq_class_id)
					{
						process_alignment_group(rec_buffer, os, aln_input.header, buffer);
						rec_buffer.clear();
					}
					
					rec_buffer.emplace_back(aln_rec);
				}
				else
				{
					sam::output_record_in_parsed_order(os, aln_input.header, aln_rec, buffer);
					os << '\n';
				}
			}
		);
		
		if (should_subset_by_best_mapq && !rec_buffer.empty())
			process_alignment_group(rec_buffer, os, aln_input.header, buffer);
		
		return mc;
	}
	
	
	void append_program_info(sam::header &output_header, int const argc, char const * const * const argv)
	{
		panvc3::append_sam_program_info(
			"panvc3.subset-alignments.",
			"PanVC 3 subset_alignments",
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
		aln_input.read_header();
		
		auto const mc([&](){
			auto aln_output_fh{[&] -> lb::file_handle {
				if (args_info.output_path_arg)
					return lb::file_handle(lb::open_file_for_writing(args_info.output_path_arg, lb::writing_open_mode::CREATE));
				else
					return lb::file_handle(STDOUT_FILENO, false);
			}()};
			
			lb::file_ostream os;
			lb::open_stream_with_file_handle(os, aln_output_fh);

			auto &header(aln_input.header);
			append_program_info(header, argc, argv); // Adding our program info does not affect parsing.
		
			return process_(aln_input, os, args_info);
		}());
		
		lb::log_time(std::cerr) << "Done. Matches: " << mc.matches << ", mismatches: " << mc.mismatches << ".\n";
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
