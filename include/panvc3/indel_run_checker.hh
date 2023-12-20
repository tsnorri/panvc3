/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_INDEL_RUN_CHECKER_HH
#define PANVC3_INDEL_RUN_CHECKER_HH

#include <ostream>
#include <panvc3/cigar_adapter.hh>
#include <panvc3/range.hh>
#include <range/v3/view/slice.hpp>
#include <vector>


namespace panvc3 {
	
	template <typename t_adapter>
	class indel_run_checker_tpl
	{
	public:
		typedef std::vector <char>				sequence_vector;
		typedef t_adapter::vector_type			cigar_vector;
		typedef cigar_vector::const_iterator	cigar_vector_const_iterator;
		typedef std::pair <
			cigar_vector_const_iterator,
			cigar_vector_const_iterator
		>										cigar_vector_const_iterator_pair;
		
	private:
		enum
		{
			RUN_HAS_INSERTIONS	= 0x1,
			RUN_HAS_DELETIONS	= 0x2,
			RUN_HAS_BOTH		= RUN_HAS_INSERTIONS | RUN_HAS_DELETIONS,
			RUN_HAS_NEITHER		= 0x0
		};
	
	private:
		cigar_vector_const_iterator_pair	m_cigar_range{};
		cigar_vector_const_iterator_pair	m_cigar_realigned_range{};
		range								m_ref_range{};
		range								m_query_range{};
		std::size_t							m_ref_pos{};
		std::size_t							m_query_pos{};
		std::uint8_t						m_run_type{};

	public:
		void reset(cigar_vector const &cv, std::size_t const ref_pos);
		bool find_next_range_for_realigning();
		range reference_range() const { return m_ref_range; }
		range query_range() const { return m_query_range; }
		std::size_t reference_position() const { return m_ref_pos; }
		std::size_t query_position() const { return m_query_pos; }
		cigar_vector_const_iterator_pair cigar_realigned_range() const { return m_cigar_realigned_range; }
		
	private:
		void update_ranges(std::size_t const ref_pos, std::size_t const query_pos);
	};
	
	
	typedef indel_run_checker_tpl <cigar_adapter_libbio>	indel_run_checker_libbio;
	typedef indel_run_checker_tpl <cigar_adapter_seqan3>	indel_run_checker_seqan3;
	
	
	template <typename t_adapter>
	void indel_run_checker_tpl <t_adapter>::reset(cigar_vector const &cv, std::size_t const ref_pos)
	{
		m_cigar_range = {cv.begin(), cv.end()};
		m_ref_pos = ref_pos;
		m_query_pos = 0;
		m_ref_range = range(ref_pos, 0);
		m_query_range = range(0, 0);
		m_run_type = 0x0;
	}
	
	
	template <typename t_adapter>
	void indel_run_checker_tpl <t_adapter>::update_ranges(std::size_t const ref_pos, std::size_t const query_pos)
	{
		m_ref_range.update_length(ref_pos);
		m_query_range.update_length(query_pos);
	}
	
	
	template <typename t_adapter>
	bool indel_run_checker_tpl <t_adapter>::find_next_range_for_realigning()
	{
		// Process the CIGAR string while maintaining state in m_run_type and in the position variables.
		// This might be somewhat simpler with coroutines.
		t_adapter adapter{};
		while (m_cigar_range.first != m_cigar_range.second)
		{
			// Copies of the existing values for use in report_run.
			auto const ref_pos(m_ref_pos);
			auto const query_pos(m_query_pos);
			auto const run_type(m_run_type);
			
			// If the run has not started yet, update its location.
			if (RUN_HAS_NEITHER == m_run_type)
			{
				m_ref_range.location = m_ref_pos;
				m_query_range.location = m_query_pos;
				m_cigar_realigned_range.first = m_cigar_range.first;
			}
			
			auto const operation(adapter.operation_(*m_cigar_range.first));
			switch (operation)
			{
				case 'H':	// Hard clipping, consumes nothing.
				case 'P':	// Padding (silent deletion from padded reference), consumes nothing.
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
				
				case 'I':	// Insertion, consumes query.
					m_query_pos += adapter.count(*m_cigar_range.first);
					m_run_type |= RUN_HAS_INSERTIONS;
					break;
				
				case 'D':	// Deletion, consumes reference.
					m_ref_pos += adapter.count(*m_cigar_range.first);
					m_run_type |= RUN_HAS_DELETIONS;
					break;
				
				case 'S':	// Soft clipping, consumes query.
					m_query_pos += adapter.count(*m_cigar_range.first);
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
					
				case 'N':	// Skipped region, consumes reference. (In SAMv1, this is only relevant in mRNA-to-genome alignments.)
					m_ref_pos += adapter.count(*m_cigar_range.first);
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
					
				case 'M':	// Match or mismatch, consumes both.
				case '=':	// Match, consumes both.
				case 'X':	// Mismatch, consumes both.
					m_ref_pos += adapter.count(*m_cigar_range.first);
					m_query_pos += adapter.count(*m_cigar_range.first);
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
				
				default:
					libbio_fail("Unexpected CIGAR operation “", operation, "”");
					break;
			}
			
			++m_cigar_range.first;
			continue;
			
		report_run:
			update_ranges(ref_pos, query_pos);
			m_cigar_realigned_range.second = m_cigar_range.first;
			++m_cigar_range.first;
			return true;
		}
		
		// Handle the final run.
		if (RUN_HAS_BOTH == m_run_type)
		{
			m_run_type = RUN_HAS_NEITHER;
			update_ranges(m_ref_pos, m_query_pos);
			m_cigar_realigned_range.second = m_cigar_range.first;
			return true;
		}
		
		return false;
	}
}

#endif
