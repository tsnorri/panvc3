/*
 * Copyright (c) 2022-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <panvc3/cigar.hh>

namespace lb	= libbio;


namespace {
	
	template <typename t_vector, typename t_adapter>
	void collapse_cigar_operations_(/* inout */ t_vector &ops, t_adapter &&adapter)
	{
		// Essentially a conditional copy.
		auto const begin(ops.begin());
		auto it(begin);
		auto it_(it + 1);
		auto const end(ops.end());

		// If the operation of the item at the reading position (*it_) matches that of
		// the current item at the writing position (*it), add to the count. Otherwise,
		// copy to the next place (i.e. *(it + 1) = *it_).
		while (it_ != end)
		{
			auto &curr_item(*it);
			auto const curr_op(adapter.operation(curr_item));
			auto const next_item(*it_);
			auto const next_op(adapter.operation(next_item));

			if (curr_op == next_op)
			{
				auto const curr_count(adapter.count(curr_item));
				auto const next_count(adapter.count(next_item));
				adapter.assign(curr_item, curr_count + next_count);
				++it_;
			}
			else
			{
				auto &dst_item(*(it + 1));
				dst_item = next_item;
				++it;
				++it_;
			}
		}

		// Remove the excess items.
		ops.resize(std::distance(begin, it + 1));
	}
}


namespace panvc3 {
	
	void collapse_cigar_operations(/* inout */ cigar_adapter_libbio::vector_type &ops)
	{
		collapse_cigar_operations_(ops, cigar_adapter_libbio{});
	}
	
	
#if defined(PANVC3_USE_SEQAN3) && PANVC3_USE_SEQAN3
	void collapse_cigar_operations(/* inout */ cigar_adapter_seqan3::vector_type &ops)
	{
		collapse_cigar_operations_(ops, cigar_adapter_seqan3{});
	}
#endif
}
