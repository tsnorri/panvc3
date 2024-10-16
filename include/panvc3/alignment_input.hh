/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_ALIGNMENT_INPUT_HH
#define PANVC3_ALIGNMENT_INPUT_HH

#include <libbio/bam/in_order_streaming_reader.hh>
#include <libbio/dispatch.hh>
#include <libbio/file_handle.hh>
#include <libbio/sam.hh>
#include <semaphore>
#include <utility>					// std::forward

namespace panvc3 {
	
	struct alignment_input_delegate
	{
		virtual ~alignment_input_delegate() {}
		virtual void handle_header(libbio::sam::header &header) = 0;
		virtual void handle_alignment(libbio::sam::record &rec) = 0;
	};
}


namespace panvc3::detail {
	
	struct alignment_input
	{
		virtual ~alignment_input() {}
		virtual void run() = 0;
	};
	
	
	class sam_alignment_input final : public alignment_input
	{
	private:
		libbio::sam::reader						m_reader;
		libbio::sam::header						m_header;
		libbio::sam::file_handle_input_range_	m_input_range;
		libbio::dispatch::serial_queue_base		*m_output_queue{};
		libbio::dispatch::group					*m_group{};
		alignment_input_delegate				*m_delegate{};
		
	public:
		sam_alignment_input(
			libbio::file_handle &&fh,
			libbio::dispatch::serial_queue_base &output_queue,
			libbio::dispatch::group &group,
			alignment_input_delegate &delegate
		):
			m_input_range(std::move(fh)),
			m_output_queue(&output_queue),
			m_group(&group),
			m_delegate(&delegate)
		{
			m_input_range.prepare();
		}
		
		void run() override;
	};
	
	
	class bam_in_order_alignment_input final :	public alignment_input,
												public libbio::bam::in_order_streaming_reader_delegate
	{
	private:
		typedef std::counting_semaphore <UINT16_MAX>	semaphore_type;

	private:
		semaphore_type							m_semaphore;
		libbio::file_handle						m_input_handle;
		libbio::dispatch::queue					*m_processing_queue{};
		libbio::dispatch::group					*m_group{};
		alignment_input_delegate				*m_delegate{};
		libbio::bam::in_order_streaming_reader	m_bam_reader;
		libbio::bgzf::streaming_reader			m_bgzf_reader;
		libbio::sam::header						m_header;
		
	public:
		bam_in_order_alignment_input(
			libbio::file_handle &&fh,
			std::size_t task_count,
			libbio::dispatch::queue &processing_queue,
			libbio::dispatch::serial_queue_base &output_queue,
			libbio::dispatch::group &group,
			alignment_input_delegate &delegate
		):
			m_semaphore(task_count),
			m_input_handle(std::move(fh)),
			m_processing_queue(&processing_queue),
			m_group(&group),
			m_delegate(&delegate),
			m_bam_reader(task_count, output_queue, *m_group, *this),
			m_bgzf_reader(m_input_handle, task_count, *m_group, &m_semaphore, m_bam_reader)
		{
		}
		
		libbio::sam::header const &header() const { return m_header; }
		
		void run() override;
		
		void streaming_reader_did_parse_header(libbio::bam::in_order_streaming_reader &reader, libbio::bam::header &&hh, libbio::sam::header &&hh_) override;
		void streaming_reader_did_parse_records(libbio::bam::in_order_streaming_reader &reader, libbio::bam::record_buffer &records) override;
	};
}


namespace panvc3 {
	
	class alignment_input
	{
	private:
		std::unique_ptr <detail::alignment_input>	m_input;
		
	private:
		explicit alignment_input(std::unique_ptr <detail::alignment_input> &&input):
			m_input(std::move(input))
		{
		}
		
	public:
		static alignment_input open_file_handle(
			libbio::file_handle &&fh,
			std::size_t task_count,								// for BAM
			libbio::dispatch::queue &reading_queue,				// for BAM
			libbio::dispatch::serial_queue_base &output_queue,
			libbio::dispatch::group &group,
			alignment_input_delegate &delegate,
			bool input_is_bam
		);
		
		static alignment_input open_path_or_stdin(
			char const *path,
			std::size_t task_count,								// for BAM
			libbio::dispatch::queue &reading_queue,				// for BAM
			libbio::dispatch::serial_queue_base &output_queue,
			libbio::dispatch::group &group,
			alignment_input_delegate &delegate
		);
			
		void run() { m_input->run(); }
	};
}

#endif
