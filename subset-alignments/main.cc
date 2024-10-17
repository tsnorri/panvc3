/*
 * Copyright (c) 2022-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>												// std::max
#include <boost/sort/block_indirect_sort/block_indirect_sort.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <libbio/dispatch.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/sam.hh>
#include <libbio/utility/misc.hh>
#include <panvc3/alignment_input.hh>
#include <panvc3/utility.hh>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>
#include "cmdline.h"

namespace dispatch	= libbio::dispatch;
namespace lb		= libbio;
namespace sam		= libbio::sam;


namespace {

	void append_program_info(sam::header &output_header, std::string const &call)
	{
		panvc3::append_sam_program_info(
			"panvc3.subset-alignments.",
			"PanVC 3 subset_alignments",
			call,
			CMDLINE_PARSER_VERSION,
			output_header.programs
		);
	}


	class subset_task final : public panvc3::alignment_input_delegate
	{
	private:
		typedef sam::mapping_quality_type						mapping_quality_type;
		typedef sam::header::reference_sequence_identifier_type	reference_id_type;

	private:
		std::deque <std::string>	m_read_names;
		std::vector <sam::record>	m_rec_buffer;
		std::vector <std::size_t>	m_sorting_buffer;
		std::string					m_command_line_call;
		sam::header					*m_header{};
		std::ostream				*m_os{};
		char const					*m_expected_contig{};
		std::uint64_t				m_rec_idx{};
		std::uint64_t				m_matches{};
		std::uint64_t				m_mismatches{};
		reference_id_type			m_expected_ref_id{};
		bool						m_should_subset_by_read_id;
		bool						m_should_subset_by_best_mapq;

	private:
		void output_best_mate(
			sam::reference_id_type const mate_ref_id,
			sam::position_type const mate_pos
		);

		void process_alignment_group();

	public:
		subset_task(std::ostream &os, gengetopt_args_info const &args_info, std::string &&command_line_call):
			m_command_line_call(std::move(command_line_call)),
			m_os(&os),
			m_expected_contig(args_info.chr_arg),
			m_should_subset_by_read_id(args_info.read_id_flag),
			m_should_subset_by_best_mapq(args_info.best_mapq_flag)
		{
		}

		std::uint64_t matches() const { return m_matches; }
		std::uint64_t mismatches() const { return m_mismatches; }

		void prepare();
		void finish();

		void handle_header(sam::header &header) override;
		void handle_alignment(sam::record &aln_rec) override;
	};


	void subset_task::prepare()
	{
		// Read the read names.
		if (m_should_subset_by_read_id)
		{
			lb::log_time(std::cerr) << "Reading the read names…\n";
			std::string buffer;
			std::size_t lineno{};
			while (std::getline(std::cin, buffer))
			{
				++lineno;
				m_read_names.emplace_back(buffer);
				if (0 == lineno % 100'000'000)
					lb::log_time(std::cerr) << "Handled " << lineno << " lines…\n";
			}

			lb::log_time(std::cerr) << "Sorting…\n";
			boost::sort::block_indirect_sort(m_read_names.begin(), m_read_names.end());
			//std::sort(read_names.begin(), read_names.end());
		}
	}


	void subset_task::output_best_mate(
		sam::reference_id_type const mate_ref_id,
		sam::position_type const mate_pos
	)
	{
		if (m_rec_buffer.empty())
			return;

		// Determine again the best mapping quality.
		mapping_quality_type best_mapq{};
		for (auto const &aln_rec : m_rec_buffer)
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

		for (auto &aln_rec : m_rec_buffer)
		{
			if (aln_rec.mapq == best_mapq && aln_rec.rname_id == mate_ref_id && aln_rec.pos == mate_pos)
			{
				sam::output_record_in_parsed_order(*m_os, *m_header, aln_rec, m_sorting_buffer);
				(*m_os) << '\n';
				return;
			}
		}
	}


	void subset_task::process_alignment_group()
	{
		if (m_rec_buffer.empty())
			return;

		// Determine the best mapping quality.
		mapping_quality_type best_mapq{};
		for (auto const &aln_rec : m_rec_buffer)
		{
			auto const mapq(aln_rec.mapq);
			if (255 == mapq)
				continue;

			best_mapq = std::max(mapq, best_mapq);
		}

		// Find again.
		for (auto &aln_rec : m_rec_buffer)
		{
			if (aln_rec.mapq == best_mapq)
			{
				sam::output_record_in_parsed_order(*m_os, *m_header, aln_rec, m_sorting_buffer);
				(*m_os) << '\n';

				auto const mate_ref_id(aln_rec.rnext_id);
				auto const mate_pos(aln_rec.pnext);
				if (sam::INVALID_REFERENCE_ID == mate_ref_id || sam::INVALID_POSITION == mate_pos)
					return;

				output_best_mate(mate_ref_id, mate_pos);
				return;
			}
		}

		// Output the first by default.
		auto &aln_rec(m_rec_buffer.front());
		sam::output_record_in_parsed_order(*m_os, *m_header, aln_rec, m_sorting_buffer);
		(*m_os) << '\n';

		auto const mate_ref_id(aln_rec.rnext_id);
		auto const mate_pos(aln_rec.pnext);
		if (sam::INVALID_REFERENCE_ID == mate_ref_id || sam::INVALID_POSITION == mate_pos)
			return;

		output_best_mate(mate_ref_id, mate_pos);
	}


	void subset_task::handle_header(sam::header &header)
	{
		if (m_expected_contig)
		{
			m_expected_ref_id = m_header->find_reference(m_expected_contig);
			if (sam::INVALID_REFERENCE_ID == m_expected_ref_id)
			{
				std::cerr << "ERROR: Reference ID ‘" << m_expected_contig << "’ not found in SAM header.\n";
				std::exit(EXIT_FAILURE);
			}
		}
		else
		{
			m_expected_ref_id = sam::INVALID_REFERENCE_ID;
		}

		append_program_info(header, m_command_line_call);
		(*m_os) << header;
		m_header = &header;
	}


	void subset_task::handle_alignment(sam::record &aln_rec)
	{
		panvc3::increment_guard const guard(m_rec_idx);
		if (m_rec_idx && 0 == m_rec_idx % 10'000'000)
			lb::log_time(std::cerr) << "Processed " << m_rec_idx << " alignments…\n";

		if (! (sam::INVALID_REFERENCE_ID == m_expected_ref_id || aln_rec.rname_id == m_expected_ref_id))
		{
			++m_mismatches;
			return;
		}

		auto const &qname(aln_rec.qname);
		if (m_should_subset_by_read_id && !std::binary_search(m_read_names.begin(), m_read_names.end(), qname))
		{
			++m_mismatches;
			return;
		}

		++m_matches;

		if (m_should_subset_by_best_mapq)
		{
			if (m_rec_buffer.empty())
			{
				m_rec_buffer.emplace_back(std::move(aln_rec));
				return;
			}

			auto const &eq_class_id(m_rec_buffer.front().qname);

			if (qname != eq_class_id)
			{
				process_alignment_group();
				m_rec_buffer.clear();
			}

			m_rec_buffer.emplace_back(std::move(aln_rec));
		}
		else
		{
			sam::output_record_in_parsed_order(*m_os, *m_header, aln_rec, m_sorting_buffer);
			(*m_os) << '\n';
		}
	}


	void subset_task::finish()
	{
		if (m_should_subset_by_best_mapq && !m_rec_buffer.empty())
			process_alignment_group();
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

	dispatch::thread_pool thread_pool;

	panvc3::prepare_thread_pool_with_args(thread_pool, args_info.threads_arg);
	auto const task_count(thread_pool.max_workers() - 1); // Reader needs one thread while reading.

	dispatch::parallel_queue parallel_queue(thread_pool);

	dispatch::group group;
	auto &main_queue(dispatch::main_queue());

	auto aln_output_fh{[&] -> lb::file_handle {
		if (args_info.output_path_arg)
			return lb::file_handle(lb::open_file_for_writing(args_info.output_path_arg, lb::writing_open_mode::CREATE));
		else
			return lb::file_handle(STDOUT_FILENO, false);
	}()};

	lb::file_ostream os;
	lb::open_stream_with_file_handle(os, aln_output_fh);

	subset_task task(os, args_info, panvc3::command_line_call(argc, argv));
	auto aln_input(panvc3::alignment_input::open_path_or_stdin(
		args_info.alignments_arg,
		task_count,
		parallel_queue,
		main_queue,
		group,
		task
	));

	lb::log_time(std::cerr) << "Processing the alignment records…\n";
	parallel_queue.group_async(group, [&aln_input, &task]{
		task.prepare();
		aln_input.run(); // Does not block.
	});

	group.notify(main_queue, [&task, &main_queue]{
		task.finish();
		lb::log_time(std::cerr) << "Done. Matches: " << task.matches() << ", mismatches: " << task.mismatches() << ".\n";

		main_queue.stop();
	});

	main_queue.run();
	return EXIT_SUCCESS;
}
