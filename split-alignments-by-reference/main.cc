/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>					// std::sort
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <libbio/dispatch.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/sam.hh>
#include <libbio/utility.hh>			// lb::compare_strings_transparent
#include <optional>
#include <panvc3/alignment_input.hh>
#include <panvc3/utility.hh>
#include <range/v3/view/zip.hpp>
#include <string>						// std::getline
#include <string_view>
#include <sstream>
#include <type_traits>
#include <utility>						// std::move
#include <vector>
#include "cmdline.h"

namespace dispatch	= libbio::dispatch;
namespace lb		= libbio;
namespace rsv		= ranges::views;
namespace sam		= libbio::sam;


namespace {

	void append_program_info(sam::header &output_header, std::string const &call)
	{
		panvc3::append_sam_program_info(
			"panvc3.split-alignments-by-reference.",
			"PanVC 3 split_alignments_by_reference",
			call,
			CMDLINE_PARSER_VERSION,
			output_header.programs
		);
	}


	struct alignment_output
	{
		lb::file_handle					handle;
		lb::file_ostream				stream;
		std::size_t						records{};

		alignment_output() = default;

		bool is_open() const { return stream.is_open(); }
		operator bool() const { return is_open(); }
		void output_record(sam::header const &header, sam::record const &rec) { sam::output_record(stream, header, rec); stream << '\n'; }
	};


	struct reference_name_record
	{
		std::string	reference_name;
		std::string	new_reference_name;
		std::size_t matches{};
		std::size_t	dst_id{SIZE_MAX};

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

	typedef std::vector <reference_name_record> reference_name_record_vector;


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


	reference_name_record_vector read_reference_names(
		char const *path,
		bool const should_treat_reference_names_as_prefixes,
		bool const should_rewrite_reference_names
	)
	{
		reference_name_record_vector retval;
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


	std::string alignment_output_path(
		std::string const &basename,
		std::string const &reference_name
	)
	{
		std::stringstream path;
		path << basename;
		path << reference_name << ".sam";
		return path.str();
	}


	class split_alignments_task : public panvc3::alignment_input_delegate
	{
	private:
		sam::header										m_aln_output_header;
		std::optional <reference_name_record_vector>	m_ref_names;
		std::vector <alignment_output>					m_aln_outputs;
		std::vector <std::size_t>						m_ref_mapping;
		std::string										m_basename;
		std::string										m_command_line_call;
		std::uint64_t									m_rec_idx{};
		std::uint64_t									m_ref_id_missing{};
		std::uint64_t									m_no_ref_name_match{};
		bool											m_should_rewrite_reference_names;
		bool											m_should_treat_reference_names_as_prefixes;

	public:
		split_alignments_task(
			std::optional <reference_name_record_vector> &&ref_names,
			std::string &&command_line_call,
			bool should_rewrite_reference_names,
			bool should_treat_reference_names_as_prefixes
		):
			m_ref_names(std::move(ref_names)),
			m_command_line_call(std::move(command_line_call)),
			m_should_rewrite_reference_names(should_rewrite_reference_names),
			m_should_treat_reference_names_as_prefixes(should_treat_reference_names_as_prefixes)
		{
		}

		void finish();

		void handle_header(sam::header &header) override;
		void handle_alignment(sam::record &aln_rec) override;
	};


	void split_alignments_task::handle_header(sam::header &header)
	{
		m_aln_output_header = header; // Copy.
		append_program_info(m_aln_output_header, m_command_line_call);


		m_ref_mapping.clear();
		m_ref_mapping.resize(header.reference_sequences.size(), SIZE_MAX);

		if (m_ref_names)
		{
			auto &ref_names(*m_ref_names);

			auto it(m_aln_output_header.reference_sequences.begin());
			reference_name_record_cmp const cmp;
			std::size_t src_ref_idx{};
			std::size_t dst_ref_idx{};

			while (it != m_aln_output_header.reference_sequences.end())
			{
				auto &ref_entry(*it);
				auto const handle_not_found([&](){
					std::cerr << "WARNING: No entry that matches reference ID ‘" << ref_entry.name << "’ in the input header.\n";
					it = m_aln_output_header.reference_sequences.erase(it);
					++src_ref_idx;
				});

				// If we have a prefix or a match, it is bound to be lexicographically smaller
				// than or equal to the reference contig name.
				auto it_(std::upper_bound(ref_names.begin(), ref_names.end(), ref_entry.name, cmp));
				if (ref_names.begin() == it_)
				{
					handle_not_found();
					continue;
				}

				auto rit(it_ - 1);
				if (
					(m_should_treat_reference_names_as_prefixes && !ref_entry.name.starts_with(rit->reference_name)) ||
					(!m_should_treat_reference_names_as_prefixes && ref_entry.name != rit->reference_name)
				)
				{
					handle_not_found();
					continue;
				}

				// Found a match.
				if (m_should_rewrite_reference_names)
				{
					if (rit->matches)
						m_ref_mapping[src_ref_idx] = rit->dst_id;
					else
					{
						m_ref_mapping[src_ref_idx] = dst_ref_idx;
						rit->dst_id = dst_ref_idx;
						++dst_ref_idx;
					}

					ref_entry.name = it_->new_reference_name;
				}
				else
				{
					m_ref_mapping[src_ref_idx] = dst_ref_idx;
					++dst_ref_idx;
					ref_entry.name = it_->reference_name;
				}

				++rit->matches;
				++src_ref_idx;
				++it;
			}

			// resize() requires that the elements are MoveInsertable, which they are not.
			m_aln_outputs = std::vector <alignment_output>(dst_ref_idx);
		}
		else
		{
			m_aln_outputs = std::vector <alignment_output>(header.reference_sequences.size());
		}
	}


	void split_alignments_task::handle_alignment(sam::record &aln_rec)
	{
		++m_rec_idx;
		if (m_rec_idx && 0 == m_rec_idx % 10'000'000)
			lb::log_time(std::cerr) << "Processed " << m_rec_idx << " alignments…\n";

		if (sam::INVALID_REFERENCE_ID == aln_rec.rname_id)
		{
			++m_ref_id_missing;
			return;
		}

		auto const open_ref([this](std::size_t idx, alignment_output &output){
			auto const &ref(m_aln_output_header.reference_sequences[idx]);
			auto const path(alignment_output_path(m_basename, ref.name));
			output.handle = lb::file_handle(lb::open_file_for_writing(path, lb::writing_open_mode::CREATE));
			lb::open_stream_with_file_handle(output.stream, output.handle);

			output.stream << m_aln_output_header;
		});

		if (m_ref_names)
		{
			auto const idx(m_ref_mapping[aln_rec.rname_id]);
			if (SIZE_MAX == idx)
			{
				++m_no_ref_name_match;
				return;
			}

			auto &output(m_aln_outputs[idx]);
			if (!output)
				open_ref(idx, output);

			++output.records;
			sam::output_record(output.stream, m_aln_output_header, aln_rec);
		}
		else
		{
			auto &output(m_aln_outputs[aln_rec.rname_id]);
			if (!output)
				open_ref(aln_rec.rname_id, output);

			++output.records;
			sam::output_record(output.stream, m_aln_output_header, aln_rec);
			output.stream << '\n';
		}
	}


	void split_alignments_task::finish()
	{
		// Report the matches.
		for (auto const &[rec, output] : rsv::zip(m_aln_output_header.reference_sequences, m_aln_outputs))
		{
			if (output)
				std::cout << rec.name << '\t' << output.records << '\n';
		}

		std::cout << "Reference ID missing\t" << m_ref_id_missing << '\n';
		std::cout << "No matching reference ID\t" << m_no_ref_name_match << '\n';

		lb::log_time(std::cerr) << "Done.\n";
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);

	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.

	dispatch::thread_pool thread_pool;

	panvc3::prepare_thread_pool_with_args(thread_pool, args_info.threads_arg);
	auto const task_count(thread_pool.max_workers() - 1); // Reader needs one thread while reading.

	dispatch::parallel_queue parallel_queue(thread_pool);

	dispatch::group group;
	auto &main_queue(dispatch::main_queue());

	std::optional <reference_name_record_vector> ref_names;
	if (args_info.reference_names_given)
	{
		ref_names = read_reference_names(
			args_info.reference_names_arg,
			args_info.prefixes_flag,
			args_info.rewrite_reference_names_flag
		);
	}

	split_alignments_task task(
		std::move(ref_names),
		panvc3::command_line_call(argc, argv),
		args_info.prefixes_flag,
		args_info.rewrite_reference_names_flag
	);
	auto aln_input(panvc3::alignment_input::open_path_or_stdin(
		args_info.alignments_arg,
		task_count,
		parallel_queue,
		main_queue,
		group,
		task
	));

	lb::log_time(std::cerr) << "Processing the alignment records…\n";
	parallel_queue.group_async(group, [&aln_input]{
		aln_input.run(); // Does not block.
	});

	group.notify(main_queue, [&task, &main_queue]{
		task.finish();
		main_queue.stop();
	});

	main_queue.run();
	return EXIT_SUCCESS;
}
