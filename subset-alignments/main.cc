/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/sort/sort.hpp>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <libbio/utility/misc.hh>
#include <range/v3/view/enumerate.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <string>
#include "cmdline.h"

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace rsv	= ranges::views;


namespace {

	struct match_count
	{
		std::size_t matches{};
		std::size_t mismatches{};
	};
	
	template <typename t_aln_input, typename t_aln_output>
	match_count process_(t_aln_input &aln_input, t_aln_output &aln_output, gengetopt_args_info const &args_info)
	{
		char const *expected_ref_id(args_info.chr_arg);
		bool const should_subset_by_read_id(args_info.read_id_flag);
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
		auto const &ref_ids(aln_input.header().ref_ids());
		std::string prev_ref_id;
		for (auto const &[rec_idx, aln_rec] : rsv::enumerate(aln_input))
		{
			if (0 == (1 + rec_idx) % 10000000)
				lb::log_time(std::cerr) << "Processed " << (1 + rec_idx) << " alignments…\n";

			if (expected_ref_id)
			{
				auto const &ref_id_(aln_rec.reference_id());
				if (!ref_id_.has_value())
				{
					++mc.mismatches;
					continue;
				}

				auto const &ref_id(ref_ids[*ref_id_]);
				bool const does_match(ref_id == expected_ref_id);

				if (args_info.verbose_flag && ref_id != prev_ref_id)
				{
					std::cerr << "Reference ID: " << ref_id << " matches: " << does_match << '\n';
					prev_ref_id = ref_id;
				}

				if (!does_match)
				{
					++mc.mismatches;
					continue;
				}
			}

			auto const &id(aln_rec.id());
			if (should_subset_by_read_id && !std::binary_search(read_names.begin(), read_names.end(), id))
			{
				++mc.mismatches;
				continue;
			}
			
			++mc.matches;
			aln_output.push_back(aln_rec);
		}

		return mc;
	}


	void process(gengetopt_args_info const &args_info)
	{
		char const *aln_path(args_info.alignments_arg);
		bool const should_output_sam(args_info.output_sam_flag);
		match_count mc{};

		// Prepare the input and the output.
		fs::path const alignments_path(aln_path);
		seqan3::sam_file_input <> aln_input(alignments_path); // Reads everything by default.

		if (should_output_sam)
		{
			seqan3::sam_file_output <> aln_output(std::cout, seqan3::format_sam{});
			mc = process_(aln_input, aln_output, args_info);
		}
		else
		{
			seqan3::sam_file_output <> aln_output(std::cout, seqan3::format_bam{});
			mc = process_(aln_input, aln_output, args_info);
		}
			
		std::cout << std::flush;
		lb::log_time(std::cerr) << "Done. Matches: " << mc.matches << ", mismatches: " << mc.mismatches << ".\n";
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	process(args_info);
	
	return EXIT_SUCCESS;
}
