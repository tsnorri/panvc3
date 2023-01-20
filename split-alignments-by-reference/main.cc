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
#include <range/v3/view/transform.hpp>
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


	fs::path alignment_output_path(
		char const *basename,
		std::string const &reference_name
	)
	{
		std::stringstream path;
		if (basename) path << basename;
		path << reference_name << ".bam";
		return {path.view()};
	}
	
	
	template <typename t_aln_input>
	void process_(
		t_aln_input &aln_input,
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
		auto &aln_input_header(aln_input.header()); // ref_ids() not const.
		auto const &ref_ids(aln_input_header.ref_ids());
		typedef seqan3::sam_file_output <
			typename t_aln_input::selected_field_ids,
			seqan3::type_list <seqan3::format_sam, seqan3::format_bam>,
			std::remove_cvref_t <decltype(ref_ids)> // sam_file_output does not seem to accept a reference for this.
		> sam_file_output_type;
		std::vector <sam_file_output_type> aln_outputs;
		aln_outputs.reserve(reference_names.size());
		
		// seqan3::sam_file_output makes a non-owning reference to output reference ids
		// under some conditions, so try to make it not deallocate before the outputs.
		std::remove_cvref_t <decltype(ref_ids)> output_ref_ids;
		if (should_rewrite_reference_names)
		{
			output_ref_ids = ref_ids; // Copy.
			for (auto &ref_id : output_ref_ids)
			{
				auto const it(std::lower_bound(reference_names.begin(), reference_names.end(), ref_id, cmp));
				if (reference_names.end() == it)
				{
					std::cerr << "ERROR: No entry for reference ID '" << ref_id << "'.\n";
					std::exit(EXIT_FAILURE);
				}

				ref_id = it->new_reference_name;
			}

			for (auto const &rec : reference_names)
			{
				aln_outputs.emplace_back(
					alignment_output_path(basename, rec.reference_name),
					output_ref_ids,
					aln_input_header.ref_id_info | rsv::transform([](auto const &tup){ return std::get <0>(tup); })
				);
			}
		}
		else
		{
			for (auto const &rec : reference_names)
				aln_outputs.emplace_back(alignment_output_path(basename, rec.reference_name));
		}
		
		// Process the records.
		lb::log_time(std::cerr) << "Processing the alignment records…\n";
		std::size_t ref_id_missing{};
		std::size_t no_match{};
		for (auto const &[rec_idx, aln_rec] : rsv::enumerate(aln_input))
		{
			if (0 == (1 + rec_idx) % 10000000)
				lb::log_time(std::cerr) << "Processed " << (1 + rec_idx) << " alignments…\n";
			
			auto const &ref_id_(aln_rec.reference_id());
			if (!ref_id_.has_value())
			{
				++ref_id_missing;
				continue;
			}
			
			// If we have a prefix or a match, it is bound to be lexicographically smaller
			// than the reference contig name.
			auto const &ref_id(ref_ids[*ref_id_]);
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
		
		// Report the matches.
		for (auto const &rec : reference_names)
			std::cout << rec.reference_name << '\t' << rec.matches << '\n';
		
		std::cout << "Reference ID missing\t" << ref_id_missing << '\n';
		std::cout << "No matching reference ID\t" << no_match << '\n';
		
		lb::log_time(std::cerr) << "Done.\n";
	}
	

	void process(
		char const *aln_path,
		char const *reference_names_path,
		char const *basename,
		bool const should_treat_reference_names_as_prefixes,
		bool const should_rewrite_reference_names,
		bool const should_report_unmatched
	)
	{
		typedef seqan3::sam_file_input <
			seqan3::sam_file_input_default_traits <>,
			seqan3::fields <
				seqan3::field::seq,
				seqan3::field::id,
				seqan3::field::offset,
				seqan3::field::ref_id,
				seqan3::field::ref_offset,
				seqan3::field::alignment,
				seqan3::field::cigar,
				seqan3::field::mapq,
				seqan3::field::qual,
				seqan3::field::flag,
				seqan3::field::mate,
				seqan3::field::tags
				//seqan3::field::header_ptr // Skip this so that sam_file_output does not copy the reference name from the record.
			>
		> sam_file_input;

		if (aln_path)
		{
			sam_file_input aln_input(fs::path{aln_path});
			process_(
				aln_input,
				reference_names_path,
				basename,
				should_treat_reference_names_as_prefixes,
				should_rewrite_reference_names,
				should_report_unmatched
			);
		}
		else
		{
			std::cerr << "Reading alignments from stdin.\n";
			sam_file_input aln_input(std::cin, seqan3::format_sam{});
			process_(
				aln_input,
				reference_names_path,
				basename,
				should_treat_reference_names_as_prefixes,
				should_rewrite_reference_names,
				should_report_unmatched
			);
		}
	}


	void read_reference_names(char const *aln_path)
	{
		fs::path const alignments_path(aln_path);
		seqan3::fields <seqan3::field::ref_id> fields{};
		seqan3::sam_file_input aln_input(alignments_path, fields);

		std::set <std::string> refs;
		auto const &ref_ids(aln_input.header().ref_ids());
		for (auto const &aln_rec : aln_input)
		{
			auto const &ref_id_(aln_rec.reference_id());
			if (!ref_id_.has_value())
				continue;
			auto const &ref_id(ref_ids[*ref_id_]);
			refs.emplace(ref_id);
		}

		for (auto const &rr : refs)
			std::cout << rr << '\n';
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.

	if (args_info.read_reference_names_given)
		read_reference_names(args_info.alignments_arg);
	else
	{
		process(
			args_info.alignments_arg,
			args_info.reference_names_arg,
			args_info.basename_arg,
			args_info.prefixes_given,
			args_info.rewrite_reference_names_given,
			args_info.report_unmatched_given
		);
	}
	
	return EXIT_SUCCESS;
}
