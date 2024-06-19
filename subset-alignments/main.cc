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
#include <panvc3/utility.hh>
#include <range/v3/view/enumerate.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <string>
#include "cmdline.h"

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace rsv	= ranges::views;


namespace {
	
	typedef std::uint8_t						mapping_quality_type;
	
	
	struct match_count
	{
		std::size_t matches{};
		std::size_t mismatches{};
	};
	
	
	template <typename t_ref_id, typename t_pos, typename t_aln_record, typename t_aln_output>
	void output_best_mate(
		t_ref_id const mate_ref_id, 
		t_pos const mate_pos,
		std::vector <t_aln_record> &alignments,
		t_aln_output &&aln_output
	)
	{
		if (alignments.empty())
			return;
		
		// Determine again the best mapping quality.
		mapping_quality_type best_mapq{};
		for (auto const &aln_rec : alignments)
		{
			if (aln_rec.reference_id() != mate_ref_id)
				continue;
			
			if (aln_rec.reference_position() != mate_pos)
				continue;
			
			auto const mapq(aln_rec.mapping_quality());
			if (255 == mapq)
				continue;
			best_mapq = std::max(mapq, best_mapq);
		}
		
		for (auto &aln_rec : alignments)
		{
			if (aln_rec.mapping_quality() == best_mapq && aln_rec.reference_id() == mate_ref_id && aln_rec.reference_position() == mate_pos)
			{
				aln_rec.header_ptr() = &aln_output.header();
				aln_output.push_back(aln_rec);
				return;
			}
		}
	}
	
	
	template <typename t_aln_record, typename t_aln_output>
	void process_alignment_group(
		std::vector <t_aln_record> &alignments,
		t_aln_output &&aln_output
	)
	{
		if (alignments.empty())
			return;
		
		// Determine the best mapping quality.
		mapping_quality_type best_mapq{};
		for (auto const &aln_rec : alignments)
		{
			auto const mapq(aln_rec.mapping_quality());
			if (255 == mapq)
				continue;
			
			best_mapq = std::max(mapq, best_mapq);
		}
		
		// Find again.
		for (auto &aln_rec : alignments)
		{
			if (aln_rec.mapping_quality() == best_mapq)
			{
				auto const mate_ref_id_(aln_rec.mate_reference_id());
				auto const mate_pos_(aln_rec.mate_position());
				
				aln_rec.header_ptr() = &aln_output.header();
				aln_output.push_back(aln_rec); // Needs non-const alignment record.
				
				if (! (mate_ref_id_ && mate_pos_))
					return;
				
				output_best_mate(*mate_ref_id_, *mate_pos_, alignments, aln_output);
				return;
			}
		}
		
		// Output the first by default.
		auto &aln_rec(alignments.front());
		auto const mate_ref_id_(aln_rec.mate_reference_id());
		auto const mate_pos_(aln_rec.mate_position());
		
		aln_rec.header_ptr() = &aln_output.header();
		aln_output.push_back(aln_rec);
		
		if (! (mate_ref_id_ && mate_pos_))
			return;
		
		output_best_mate(*mate_ref_id_, *mate_pos_, alignments, aln_output);
	}
	
	
	template <typename t_aln_input, typename t_aln_output>
	match_count process_(t_aln_input &aln_input, t_aln_output &aln_output, gengetopt_args_info const &args_info)
	{
		typedef std::remove_cvref_t <t_aln_input>	input_type;
		typedef typename input_type::record_type	record_type;
		
		char const *expected_ref_id(args_info.chr_arg);
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
		auto const &ref_ids(aln_input.header().ref_ids());
		std::string prev_ref_id;
		std::vector <record_type> rec_buffer;
		for (auto const &[rec_idx, aln_rec] : rsv::enumerate(aln_input))
		{
			if (rec_idx && 0 == rec_idx % 10'000'000)
				lb::log_time(std::cerr) << "Processed " << rec_idx << " alignments…\n";
			
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
			if (should_subset_by_best_mapq)
			{
				if (rec_buffer.empty())
				{
					rec_buffer.emplace_back(std::move(aln_rec));
					continue;
				}
				
				auto const &eq_class_id(rec_buffer.front().id());
			
				if (id != eq_class_id)
				{
					process_alignment_group(rec_buffer, aln_output);
					rec_buffer.clear();
				}
				
				rec_buffer.emplace_back(std::move(aln_rec));
			}
			else
			{
				aln_rec.header_ptr() = &aln_output.header();
				aln_output.push_back(aln_rec);
			}
		}
		
		if (should_subset_by_best_mapq && !rec_buffer.empty())
			process_alignment_group(rec_buffer, aln_output);
		
		return mc;
	}
	
	
	template <typename t_header>
	void append_program_info(t_header &output_header, int const argc, char const * const * const argv)
	{
		typedef typename t_header::program_info_t program_info_type;
		panvc3::append_sam_program_info_seqan3(
			"panvc3.subset-alignments.",
			"PanVC 3 subset_alignments",
			argc,
			argv,
			CMDLINE_PARSER_VERSION,
			output_header.program_infos
		);
	}
	
	
	void process(gengetopt_args_info const &args_info, int const argc, char const * const * const argv)
	{
		// Open the SAM input and output.
		typedef seqan3::sam_file_input <>								input_type;
		auto aln_input{[&](){
			auto const make_input_type([&]<typename t_fmt>(t_fmt &&fmt){
				if (args_info.alignments_arg)
					return input_type(fs::path{args_info.alignments_arg});
				else
					return input_type(std::cin, std::forward <t_fmt>(fmt));
			});
			
			if (args_info.bam_input_flag)
				return make_input_type(seqan3::format_bam{});
			else
				return make_input_type(seqan3::format_sam{});
		}()};

		auto const &aln_input_header(aln_input.header());
		auto const &input_ref_ids(aln_input.header().ref_ids()); // ref_ids() not const.
		typedef std::remove_cvref_t <decltype(input_ref_ids)>			ref_ids_type;
		ref_ids_type output_ref_ids;
		typedef seqan3::sam_file_output <
			typename input_type::selected_field_ids,
			seqan3::type_list <seqan3::format_sam, seqan3::format_bam>,
			ref_ids_type
		>																output_type;
		
		auto aln_output{[&](){
			// Make sure that aln_output has some header information.
			auto const make_output_type([&]<typename t_fmt>(t_fmt &&fmt){
				ref_ids_type empty_ref_ids;
				if (args_info.output_path_arg)
				{
					return output_type(
						fs::path{args_info.output_path_arg},
						output_ref_ids,
						rsv::empty <std::size_t>()	// Reference lengths; the constructor expects a forward range.
					);
				}
				else
				{
					return output_type(
						std::cout,
						output_ref_ids,
						rsv::empty <std::size_t>(),	// Reference lengths; the constructor expects a forward range.
						std::forward <t_fmt>(fmt)
					);
				}
			});
			
			if (args_info.output_bam_flag)
				return make_output_type(seqan3::format_bam{});
			else
				return make_output_type(seqan3::format_sam{});
		}()};

		// Now that we have valid header information in aln_output, we can access it and copy the headers from aln_input.
		{
			output_ref_ids = input_ref_ids;

			auto &aln_output_header(aln_output.header());
			aln_output_header.sorting = aln_input_header.sorting;
			aln_output_header.subsorting = aln_input_header.subsorting;
			aln_output_header.grouping = aln_input_header.subsorting;
			aln_output_header.program_infos = aln_input_header.program_infos;
			aln_output_header.comments = aln_input_header.comments;
			aln_output_header.ref_id_info = aln_input_header.ref_id_info;
			aln_output_header.ref_dict = aln_input_header.ref_dict;
			aln_output_header.read_groups = aln_input_header.read_groups;
			
			append_program_info(aln_output_header, argc, argv);
		}
		
		auto const mc(process_(aln_input, aln_output, args_info));
		
		std::cout << std::flush;
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
