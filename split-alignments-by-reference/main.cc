/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <libbio/file_handling.hh>
#include <libbio/utility.hh>
#include <panvc3/alignment_input.hh>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>
#include <set>
#include <string>
#include <vector>
#include "cmdline.h"

namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace sam	= libbio::sam;


namespace {
	
	struct alignment_output
	{
		lb::file_handle					handle;
		lb::file_ostream				stream;
		
		// lb::file_ostream is not movable and so we unfortunately need this.
		explicit alignment_output(alignment_output &&other):
			handle(std::move(other.handle))
		{
			lb::open_stream_with_file_handle(stream, handle);
		}
		
		explicit alignment_output(std::string const &path):
			handle(lb::open_file_for_writing(path, lb::writing_open_mode::CREATE))
		{
			lb::open_stream_with_file_handle(stream, handle);
		}
		
		void output_record(sam::header const &header, sam::record const &rec) { sam::output_record(stream, header, rec); stream << '\n'; }
	};
	
	
	struct reference_name_record
	{
		std::string	reference_name;
		std::string	new_reference_name;
		std::size_t	matches{};
		
		explicit reference_name_record(std::string const &reference_name_):
			reference_name(reference_name_)
		{
		}

		template <typename t_string>
		reference_name_record(t_string &&reference_name_, t_string &&new_reference_name_):
			reference_name(std::forward <t_string>(reference_name_)),
			new_reference_name(std::forward <t_string>(new_reference_name_))
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
	
	
	std::vector <reference_name_record> read_reference_names(
		char const *path,
		bool const should_treat_reference_names_as_prefixes,
		bool const should_rewrite_reference_names
	)
	{
		std::vector <reference_name_record> retval;
		std::string buffer;
		
		lb::file_istream stream;
		lb::open_file_for_reading(path, stream);
		
		if (should_rewrite_reference_names)
		{
			// Read the tab-separated names.
			std::size_t lineno{};
			while (std::getline(stream, buffer))
			{
				++lineno;
				std::string_view const line(buffer);
				auto const tab_pos(line.find('\t'));
				if (std::string_view::npos == tab_pos)
				{
					std::cerr << "ERROR: Unable to parse reference name on line " << lineno << ".\n";
					std::exit(EXIT_FAILURE);
				}

				auto const rname(line.substr(0, tab_pos));
				auto const new_rname(line.substr(1 + tab_pos));
				retval.emplace_back(rname, new_rname);
			}
		}
		else
		{
			while (std::getline(stream, buffer))
				retval.emplace_back(buffer);
		}
		
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


	std::string alignment_output_path(
		char const *basename,
		std::string const &reference_name
	)
	{
		std::stringstream path;
		if (basename) path << basename;
		path << reference_name << ".sam";
		return path.str();
	}
	
	
	void process(
		panvc3::alignment_input &aln_input,
		char const *reference_names_path,
		char const *basename,
		bool const should_treat_reference_names_as_prefixes,
		bool const should_rewrite_reference_names,
		bool const should_report_unmatched
	)
	{
		// Non-const needed for statistics.
		auto reference_names(read_reference_names(
			reference_names_path,
			should_treat_reference_names_as_prefixes,
			should_rewrite_reference_names
		));
		
		reference_name_record_cmp const cmp;
		
		// Prepare the output.
		auto const &aln_input_header(aln_input.header);
		auto aln_output_header(aln_input_header);
		std::vector <alignment_output> aln_outputs;
		aln_outputs.reserve(reference_names.size());
		
		if (should_rewrite_reference_names)
		{
			for (auto &ref_entry : aln_output_header.reference_sequences)
			{
				auto const it(std::lower_bound(reference_names.begin(), reference_names.end(), ref_entry.name, cmp));
				if (reference_names.end() == it)
				{
					std::cerr << "ERROR: No entry for reference ID ‘" << ref_entry.name << "’.\n";
					std::exit(EXIT_FAILURE);
				}
				
				ref_entry.name = it->new_reference_name;
			}
		}
		
		// Open the output files.
		for (auto const &ref_entry : aln_output_header.reference_sequences)
			aln_outputs.emplace_back(alignment_output_path(basename, ref_entry.name));
		
		// Process the records.
		lb::log_time(std::cerr) << "Processing the alignment records…\n";
		std::size_t rec_idx{SIZE_MAX};
		std::size_t ref_id_missing{};
		std::size_t no_match{};
		aln_input.read_records(
			[
				&reference_names,
				&aln_outputs,
				&aln_input_header,
				&aln_output_header,
				&cmp,
				&rec_idx,
				&ref_id_missing,
				&no_match,
				should_report_unmatched,
				should_treat_reference_names_as_prefixes
			](sam::record const &aln_rec){
				++rec_idx;
				if (0 == (1 + rec_idx) % 10000000)
					lb::log_time(std::cerr) << "Processed " << (1 + rec_idx) << " alignments…\n";
				
				auto const rname_id(aln_rec.rname_id);
				if (sam::INVALID_REFERENCE_ID == rname_id)
				{
					++ref_id_missing;
					return;
				}
				
				// If we have a prefix or a match, it is bound to be lexicographically smaller
				// than or equal to the reference contig name.
				auto const &ref_entry(aln_input_header.reference_sequences[rname_id]);
				auto const it(std::upper_bound(reference_names.begin(), reference_names.end(), ref_entry.name, cmp));
				
				if (reference_names.begin() == it)
				{
					++no_match;
					if (should_report_unmatched)
						report_unmatched(ref_entry.name);
				
					return;
				}
				
				auto const rit(it - 1);
				auto const &reference_name(rit->reference_name);
				if (
					(should_treat_reference_names_as_prefixes  && !ref_entry.name.starts_with(reference_name)) ||
					(!should_treat_reference_names_as_prefixes && ref_entry.name != reference_name)
				)
				{
					++no_match;
					if (should_report_unmatched)
						report_unmatched(ref_entry.name);
				
					return;
				}
				
				// Found a matchig prefix.
				++rit->matches;
				auto const idx(std::distance(reference_names.begin(), rit));
				auto &aln_output(aln_outputs[idx]);
				aln_output.output_record(aln_output_header, aln_rec);
			}
		);
		
		// Report the matches.
		for (auto const &rec : reference_names)
			std::cout << rec.reference_name << '\t' << rec.matches << '\n';
		
		std::cout << "Reference ID missing\t" << ref_id_missing << '\n';
		std::cout << "No matching reference ID\t" << no_match << '\n';
		
		lb::log_time(std::cerr) << "Done.\n";
	}
	
	
	void read_reference_names(panvc3::alignment_input &aln_input, bool const should_check_alignments)
	{
		auto const &refs(aln_input.header.reference_sequences);
		
		if (should_check_alignments)
		{
			// Check if any of the alignments actually refer to the reference name.
			std::set <sam::header::reference_sequence_identifier_type> seen_ref_ids;
			aln_input.read_records(
				[
					&seen_ref_ids
				](sam::record const &aln_rec){
					if (sam::INVALID_REFERENCE_ID == aln_rec.rname_id)
						return;
					seen_ref_ids.emplace(aln_rec.rname_id);
				}
			);
			
			for (auto const ref_id : seen_ref_ids)
			{
				auto const &ref(refs[ref_id]);
				std::cout << ref.name << '\n';
			}
		}
		else
		{
			for (auto const &ref : refs)
				std::cout << ref.name << '\n';
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.

	auto aln_input(panvc3::alignment_input::open_path_or_stdin(args_info.alignments_arg));
	
	if (args_info.read_reference_names_given)
		read_reference_names(aln_input, args_info.only_used_given);
	else
	{
		process(
			aln_input,
			args_info.reference_names_arg,
			args_info.basename_arg,
			args_info.prefixes_given,
			args_info.rewrite_reference_names_given,
			args_info.report_unmatched_given
		);
	}
	
	return EXIT_SUCCESS;
}
