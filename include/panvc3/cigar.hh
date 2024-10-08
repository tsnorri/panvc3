/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_CIGAR_HH
#define PANVC3_CIGAR_HH

#include <panvc3/cigar_adapter.hh>
#include <utility>					// std::swap


namespace panvc3 {
	
	// Combine adjacent matching CIGAR operations.
	void collapse_cigar_operations(/* inout */ cigar_adapter_libbio::vector_type &ops);
	
	
	template <typename t_adapter>
	class cigar_buffer_tpl
	{
	public:
		typedef	t_adapter::vector_type		vector_type;
		typedef	t_adapter::run_type			run_type;
		typedef	t_adapter::operation_type	operation_type;
		typedef	t_adapter::count_type		count_type;
		
	protected:
		vector_type	m_operations;
		run_type	m_current_op;
		bool		m_is_first{true};
		
	public:
		void push_back(operation_type const op, count_type const count = 1);
		void finish();
		void clear();
		void swap_buffer(vector_type &ops) { using std::swap; swap(ops, m_operations); }
		
		vector_type &operations() { return m_operations; }
		vector_type const &operations() const { return m_operations; }
	};
	
	typedef cigar_buffer_tpl <cigar_adapter_libbio> cigar_buffer_libbio;
	
	
	template <typename t_adapter>
	void cigar_buffer_tpl <t_adapter>::clear()
	{
		m_is_first = true;
		m_operations.clear();
	}
	
	
	template <typename t_adapter>
	void cigar_buffer_tpl <t_adapter>::push_back(operation_type const op, count_type count)
	{
		if (0 == count)
			return;
		
		t_adapter adapter{};
		
		if (m_is_first)
		{
			m_is_first = false;
			adapter.assign(m_current_op, op);
			adapter.assign(m_current_op, count);
		}
		else if (adapter.operation(m_current_op) == op)
		{
			count += adapter.count(m_current_op);
			adapter.assign(m_current_op, count);
		}
		else
		{
			if (adapter.count(m_current_op))
				m_operations.emplace_back(m_current_op);
			adapter.assign(m_current_op, op);
			adapter.assign(m_current_op, count);
		}
	}
	
	
	template <typename t_adapter>
	void cigar_buffer_tpl <t_adapter>::finish()
	{
		if (!m_is_first && t_adapter{}.count(m_current_op))
			m_operations.emplace_back(m_current_op);
	}
}


#if defined(PANVC3_USE_SEQAN3) && PANVC3_USE_SEQAN3
namespace panvc3 {
	void collapse_cigar_operations(/* inout */ cigar_adapter_seqan3::vector_type &ops);
	typedef cigar_buffer_tpl <cigar_adapter_seqan3> cigar_buffer_seqan3;
}
#endif

#endif
