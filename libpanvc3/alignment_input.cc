/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/algorithm/string/predicate.hpp>	// boost::iequals
#include <libbio/file_handling.hh>
#include <panvc3/alignment_input.hh>
#include <unistd.h>								// STDIN_FILENO

namespace bam		= libbio::bam;
namespace dispatch	= libbio::dispatch;
namespace lb		= libbio;
namespace sam		= libbio::sam;


namespace panvc3 {
	
	alignment_input alignment_input::open_file_handle(
		libbio::file_handle &&fh,
		std::size_t task_count,								// for BAM
		dispatch::queue &processing_queue,					// for BAM
		dispatch::serial_queue_base &output_queue,
		dispatch::group &group,
		alignment_input_delegate &delegate,
		bool input_is_bam
	)
	{
		if (input_is_bam)
			return alignment_input{std::make_unique <detail::bam_in_order_alignment_input>(std::move(fh), task_count, processing_queue, output_queue, group, delegate)};
		else
			return alignment_input{std::make_unique <detail::sam_alignment_input>(std::move(fh), output_queue, group, delegate)};
	}
	
	
	alignment_input alignment_input::open_path_or_stdin(
		char const *path,
		std::size_t task_count,								// for BAM
		dispatch::queue &processing_queue,					// for BAM
		dispatch::serial_queue_base &output_queue,
		dispatch::group &group,
		alignment_input_delegate &delegate
	)
	{
		if (!path)
			return open_file_handle(lb::file_handle(STDIN_FILENO, false), task_count, processing_queue, output_queue, group, delegate, false);
		
		// Determine from the file name extension if the input is BAM.
		lb::file_handle fh(lb::open_file_for_reading(path));
		auto const path_length(strlen(path));
		if (4 <= path_length)
		{
			std::string_view const path_suffix(path + path_length - 4, 4);
			if (boost::iequals(path_suffix, std::string_view{".bam"}))
				return open_file_handle(std::move(fh), task_count, processing_queue, output_queue, group, delegate, true);
		}
		
		return open_file_handle(std::move(fh), task_count, processing_queue, output_queue, group, delegate, false);
	}
}


namespace panvc3::detail {
	
	void sam_alignment_input::run()
	{
		m_output_queue->group_async(*m_group, [this](){
			m_reader.read_header(m_header, m_input_range);
			m_delegate->handle_header(m_header);
			m_reader.read_records(m_header, m_input_range, [this](sam::record &aln_rec){ m_delegate->handle_alignment(aln_rec); });
		});
	}
	
	
	void bam_in_order_alignment_input::run()
	{
		m_bgzf_reader.run(*m_processing_queue);
	}
	
	
	void bam_in_order_alignment_input::streaming_reader_did_parse_header(
		bam::in_order_streaming_reader &reader,
		bam::header &&hh,
		sam::header &&hh_
	)
	{
		m_header = std::move(hh_);
		m_delegate->handle_header(m_header);
	}
	
	
	void bam_in_order_alignment_input::streaming_reader_did_parse_records(
		bam::in_order_streaming_reader &reader,
		bam::record_buffer &records
	)
	{
		for (auto &aln_rec : records)
			m_delegate->handle_alignment(aln_rec);
		
		m_semaphore.release();
	}
}
