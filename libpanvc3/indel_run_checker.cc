/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/indel_run_checker.hh>

namespace lb	= libbio;


namespace panvc3 {
	
	void indel_run_checker::reset(panvc3::cigar_vector const &cigar_vector, std::size_t const ref_pos)
	{
		m_cigar_it = cigar_vector.begin();
		m_cigar_end = cigar_vector.end();
		m_ref_pos = ref_pos;
		m_query_pos = 0;
		m_ref_range = range(ref_pos, 0);
		m_query_range = range(0, 0);
		m_run_type = 0x0;
	}
	
	
	void indel_run_checker::update_ranges(std::size_t const ref_pos, std::size_t const query_pos)
	{
		m_ref_range.update_length(ref_pos);
		m_query_range.update_length(query_pos);
	}
	
	
	bool indel_run_checker::find_next_range_for_realigning()
	{
		// Process the CIGAR string while maintaining state in m_run_type and in the position variables.
		// This might be somewhat simpler with coroutines.
		using seqan3::get;
		
		while (m_cigar_it != m_cigar_end)
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
				m_cigar_realigned_range_begin = m_cigar_it;
			}
			
			auto const operation(get <1>(*m_cigar_it).to_char());
			switch (operation)
			{
				case 'H':	// Hard clipping, consumes nothing.
				case 'P':	// Padding (silent deletion from padded reference), consumes nothing.
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
				
				case 'I':	// Insertion, consumes query.
					m_query_pos += get <0>(*m_cigar_it);
					m_run_type |= RUN_HAS_INSERTIONS;
					break;
				
				case 'D':	// Deletion, consumes reference.
					m_ref_pos += get <0>(*m_cigar_it);
					m_run_type |= RUN_HAS_DELETIONS;
					break;
				
				case 'S':	// Soft clipping, consumes query.
					m_query_pos += get <0>(*m_cigar_it);
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
					
				case 'N':	// Skipped region, consumes reference. (In SAMv1, this is only relevant in mRNA-to-genome alignments.)
					m_ref_pos += get <0>(*m_cigar_it);
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
					
				case 'M':	// Match or mismatch, consumes both.
				case '=':	// Match, consumes both.
				case 'X':	// Mismatch, consumes both.
					m_ref_pos += get <0>(*m_cigar_it);
					m_query_pos += get <0>(*m_cigar_it);
					m_run_type = RUN_HAS_NEITHER;
					if (RUN_HAS_BOTH == run_type)
						goto report_run;
					break;
				
				default:
					libbio_fail("Unexpected CIGAR operation “", operation, "”");
					break;
			}
			
			++m_cigar_it;
			continue;
			
		report_run:
			update_ranges(ref_pos, query_pos);
			++m_cigar_it;
			return true;
		}
		
		// Handle the final run.
		if (RUN_HAS_BOTH == m_run_type)
		{
			m_run_type = RUN_HAS_NEITHER;
			update_ranges(m_ref_pos, m_query_pos);
			return true;
		}
		
		return false;
	}
}
