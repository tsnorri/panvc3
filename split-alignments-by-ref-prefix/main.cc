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
	
	struct reference_name_record
	{
		std::string	reference_name;
		std::size_t	matches{};
		
		reference_name_record(std::string	const &reference_name_):
			reference_name(reference_name_)
		{
		}
		
		bool operator<(reference_name_record const &other) const
		{
			return reference_name < other.reference_name;
		}
	};
	
	
	// Transparent comparator for reference_name_records and strings.
	struct reference_name_record_cmp
	{
		using is_transparent = std::true_type;
		
		bool operator()(reference_name_record const lhs, reference_name_record const &rhs) const
		{
			return lhs < rhs;
		}
		
		template <typename t_string>
		bool operator()(reference_name_record const &lhs, t_string const &rhs) const
		{
			lb::compare_strings_transparent cmp;
			return cmp(lhs.reference_name, rhs);
		}
		
		template <typename t_string>
		bool operator()(t_string const &lhs, reference_name_record const &rhs) const
		{
			lb::compare_strings_transparent cmp;
			return cmp(lhs, rhs.reference_name);
		}
	};
	
	
	std::vector <reference_name_record> read_reference_names(char const *path, bool const should_treat_reference_names_as_prefixes)
	{
		std::vector <reference_name_record> retval;
		std::string buffer;
		
		lb::file_istream stream;
		lb::open_file_for_reading(path, stream);
		
		while (std::getline(stream, buffer))
			retval.emplace_back(buffer);
		
		if (retval.empty())
		{
			std::cerr << "ERROR: The reference name list was empty.\n";
			std::exit(EXIT_FAILURE);
		}
		
		// There is at least one reference name.
		std::sort(retval.begin(), retval.end());
		
		auto it(retval.begin());
		auto const end(retval.end() - 1);
		while (it != end)
		{
			if (should_treat_reference_names_as_prefixes)
			{
				if ((it + 1)->reference_name.starts_with(it->reference_name))
				{
					std::cerr << "ERROR: The contig prefixes must be prefix-free but “" << it->reference_name << "” is a prefix of “" << (it + 1)->reference_name << "”.\n";
					std::exit(EXIT_FAILURE);
				}
			}
			else
			{
				if ((it + 1)->reference_name == it->reference_name)
				{
					std::cerr << "ERROR: Found duplicate contig name: " << it->reference_name << ".\n";
					std::exit(EXIT_FAILURE);
				}
			}
		
			++it;
		}
		
		return retval;
	}
	
	
	inline void report_unmatched(std::string const &ref_id)
	{
		std::cerr << "WARNING: No reference name found that would match “" << ref_id << "”.\n";
	}
	
	
	void process(
		char const *aln_path,
		char const *reference_names_path,
		bool const should_treat_reference_names_as_prefixes,
		bool const should_report_unmatched
	)
	{
		// Non-const needed for statistics.
		auto reference_names(read_reference_names(reference_names_path, should_treat_reference_names_as_prefixes));
		
		fs::path const alignments_path(aln_path);
		seqan3::sam_file_input <> aln_input(alignments_path); // Reads everything by default.
		
		// Prepare the output.
		std::vector <seqan3::sam_file_output <>> aln_outputs;
		aln_outputs.reserve(reference_names.size());
		for (auto const &rec : reference_names)
		{
			std::string path_str(rec.reference_name);
			path_str += ".bam";
			fs::path const path(path_str);
			aln_outputs.emplace_back(path);
		}
		
		// Process the records.
		lb::log_time(std::cerr) << "Processing the alignment records…\n";
		auto const &ref_ids(aln_input.header().ref_ids());
		std::size_t ref_id_missing{};
		std::size_t no_match{};
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
			
			// If we have a prefix or a match, it is bound to be lexicographically smaller
			// than the reference contig name.
			auto const &ref_id(ref_ids[*ref_id_]);
			reference_name_record_cmp cmp;
			auto const it(std::upper_bound(reference_names.begin(), reference_names.end(), ref_id, cmp));
			
			if (reference_names.begin() == it)
			{
				++no_match;
				if (should_report_unmatched)
					report_unmatched(ref_id);
				
				continue;
			}
			
			auto const rit(it - 1);
			auto const &reference_name(rit->reference_name);
			if (
				(should_treat_reference_names_as_prefixes  && !ref_id.starts_with(reference_name)) ||
				(!should_treat_reference_names_as_prefixes && ref_id != reference_name)
			)
			{
				++no_match;
				if (should_report_unmatched)
					report_unmatched(ref_id);
				
				continue;
			}
			
			// Found a matchig prefix.
			++rit->matches;
			auto const idx(std::distance(reference_names.begin(), rit));
			auto &output_file(aln_outputs[idx]);
			output_file.push_back(aln_rec);
		}
		
		// Report matches.
		for (auto const &rec : reference_names)
			std::cout << rec.reference_name << '\t' << rec.matches << '\n';
		
		std::cout << "Reference ID missing\t" << ref_id_missing << '\n';
		std::cout << "No matching reference ID\t" << no_match << '\n';
		
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
		args_info.reference_names_arg,
		args_info.prefixes_flag,
		args_info.report_unmatched_flag
	);
	
	return EXIT_SUCCESS;
}
