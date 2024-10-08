/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_CIGAR_ADAPTER_HH
#define PANVC3_CIGAR_ADAPTER_HH

#include <libbio/sam/cigar.hh>
#include <panvc3/utility.hh>				// type_list_to_tuple_t
#include <vector>


#if defined(PANVC3_USE_SEQAN3) && PANVC3_USE_SEQAN3
#include <seqan3/alphabet/cigar/cigar.hpp>

namespace panvc3::cigar_adapters::seqan3::detail {
	
	using ::seqan3::operator""_cigar_operation;
	
	struct cigar_operations
	{
		constexpr static inline auto const alignment_match_op{'M'_cigar_operation};
		constexpr static inline auto const insertion_op{'I'_cigar_operation};
		constexpr static inline auto const deletion_op{'D'_cigar_operation};
		constexpr static inline auto const skipped_region_op{'N'_cigar_operation};
		constexpr static inline auto const soft_clipping_op{'S'_cigar_operation};
		constexpr static inline auto const hard_clipping_op{'H'_cigar_operation};
		constexpr static inline auto const padding_op{'P'_cigar_operation};
		constexpr static inline auto const sequence_match_op{'='_cigar_operation};
		constexpr static inline auto const sequence_mismatch_op{'X'_cigar_operation};
	};
}

#endif


namespace panvc3 {
	
	struct cigar_adapter_libbio
	{
		typedef libbio::sam::cigar_run				run_type;
		typedef std::vector <run_type>				vector_type;
		typedef libbio::sam::cigar_operation		operation_type;
		typedef libbio::sam::cigar_run::count_type	count_type;
		
		constexpr static inline auto const alignment_match_op	{libbio::sam::cigar_operation::alignment_match};
		constexpr static inline auto const insertion_op			{libbio::sam::cigar_operation::insertion};
		constexpr static inline auto const deletion_op			{libbio::sam::cigar_operation::deletion};
		constexpr static inline auto const skipped_region_op	{libbio::sam::cigar_operation::skipped_region};
		constexpr static inline auto const soft_clipping_op		{libbio::sam::cigar_operation::soft_clipping};
		constexpr static inline auto const hard_clipping_op		{libbio::sam::cigar_operation::hard_clipping};
		constexpr static inline auto const padding_op			{libbio::sam::cigar_operation::padding};
		constexpr static inline auto const sequence_match_op	{libbio::sam::cigar_operation::sequence_match};
		constexpr static inline auto const sequence_mismatch_op	{libbio::sam::cigar_operation::sequence_mismatch};
		
		constexpr count_type count(run_type const &ci) const { return ci.count(); }
		constexpr operation_type operation(run_type const &ci) const { return ci.operation(); }
		constexpr char operation_(run_type const &ci) const { return to_char(operation(ci)); }
		constexpr auto unpack(run_type const &ci) const { return std::make_tuple(count(ci), operation(ci), operation_(ci)); }
		
		constexpr auto assign(run_type &ci, operation_type const op) const { ci.assign(op); }
		constexpr auto assign(run_type &ci, count_type const count) const { ci.assign(count); }
	};
	
	
#if defined(PANVC3_USE_SEQAN3) && PANVC3_USE_SEQAN3
	struct cigar_adapter_seqan3 : public cigar_adapters::seqan3::detail::cigar_operations
	{
		typedef seqan3::cigar											run_type;
		typedef std::vector <run_type>									vector_type;
		typedef seqan3::cigar::operation								operation_type;
		
		typedef type_list_to_tuple_t <run_type::seqan3_required_types>	cigar_component_types;
		typedef std::tuple_element_t <0, cigar_component_types>			count_type;
		
		constexpr count_type count(run_type const &ci) const { using seqan3::get; return get <0>(ci); }
		constexpr operation_type operation(run_type const &ci) const { using seqan3::get; return get <1>(ci); }
		constexpr char operation_(run_type const &ci) const { using seqan3::get; return get <1>(ci).to_char(); }
		constexpr auto unpack(run_type const &ci) const { return std::make_tuple(count(ci), operation(ci), operation_(ci)); }
		
		constexpr auto assign(run_type &ci, operation_type const op) const { ci = op; }
		constexpr auto assign(run_type &ci, count_type const count) const { ci = count; }
	};
#endif
}

#endif
