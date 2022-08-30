/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <libbio/file_handling.hh>
#include <libbio/utility.hh>
#include <range/v3/view/enumerate.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <string>
#include <vector>
#include "cmdline.h"

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace rsv	= ranges::views;


namespace {
	
	std::vector <std::string> read_prefixes(char const *path)
	{
		std::vector <std::string> retval;
		std::string buffer;
		
		lb::file_istream stream;
		lb::open_file_for_reading(path, stream);
		
		while (std::getline(stream, buffer))
			retval.emplace_back(buffer);
		
		if (retval.empty())
		{
			std::cerr << "ERROR: The prefix list was empty.\n";
			std::exit(EXIT_FAILURE);
		}
		
		// There is at least one prefix.
		std::sort(retval.begin(), retval.end());
		
		auto it(retval.begin());
		auto const end(retval.end() - 1);
		while (it != end)
		{
			if ((it + 1)->starts_with(*it))
			{
				std::cerr << "ERROR: The contig prefixes must be prefix-free but “" << *it << "” is a prefix of “" << *(it + 1) << "”.\n";
				std::exit(EXIT_FAILURE);
			}
			
			++it;
		}
		
		return retval;
	}
	
	
	inline void report_unmatched(std::string const &ref_id)
	{
		std::cerr << "WARNING: No prefix found for reference “" << ref_id << "”.\n";
	}
	
	
	void process(char const *aln_path, char const *prefixes_path, bool const should_report_unmatched)
	{
		auto const prefixes(read_prefixes(prefixes_path));
		
		fs::path const alignments_path(aln_path);
		seqan3::sam_file_input <> aln_input(alignments_path); // Reads everything by default.
		
		// Prepare the output.
		std::vector <seqan3::sam_file_output <>> aln_outputs;
		aln_outputs.reserve(prefixes.size());
		for (auto const &prefix : prefixes)
		{
			std::string path_str(prefix);
			path_str += ".bam";
			fs::path const path(path_str);
			aln_outputs.emplace_back(path);
		}
		
		// Process the records.
		lb::log_time(std::cerr) << "Processing the alignment records…\n";
		auto const &ref_ids(aln_input.header().ref_ids());
		std::size_t ref_id_missing{};
		std::size_t no_prefix_found{};
		for (auto const &[rec_idx, aln_rec] : rsv::enumerate(aln_input))
		{
			if (0 == (1 + rec_idx) % 100000)
				lb::log_time(std::cerr) << "Processed " << (1 + rec_idx) << " alignments…\n";
			
			auto const ref_id_(aln_rec.reference_id());
			if (!ref_id_.has_value())
			{
				++ref_id_missing;
				continue;
			}
			
			// If we have a prefix, it is bound to be lexicographically smaller than the reference contig name.
			auto const &ref_id(ref_ids[*ref_id_]);
			auto const it(std::upper_bound(prefixes.begin(), prefixes.end(), ref_id));
			
			if (prefixes.begin() == it)
			{
				++no_prefix_found;
				if (should_report_unmatched)
					report_unmatched(ref_id);
					
				continue;
			}
			
			auto const pit(it - 1);
			auto const &prefix(*pit);
			if (!ref_id.starts_with(prefix))
			{
				++no_prefix_found;
				if (should_report_unmatched)
					report_unmatched(ref_id);
				
				continue;
			}
			
			// Found a matchig prefix.
			auto const idx(std::distance(prefixes.begin(), pit));
			auto &output_file(aln_outputs[idx]);
			output_file.push_back(aln_rec);
		}
		
		std::cout << "Reference ID missing\t" << ref_id_missing << '\n';
		std::cout << "No matching prefix\t" << no_prefix_found << '\n';
		
		lb::log_time(std::cerr) << "Done.\n";
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	process(
		args_info.alignments_arg,
		args_info.prefixes_arg,
		args_info.report_unmatched_flag
	);
	
	return EXIT_SUCCESS;
}
