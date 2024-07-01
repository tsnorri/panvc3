/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_ALIGNMENT_INPUT_HH
#define PANVC3_ALIGNMENT_INPUT_HH

#include <libbio/file_handle.hh>
#include <libbio/sam.hh>
#include <utility>					// std::forward


namespace panvc3 {
	
	struct alignment_input
	{
		libbio::sam::file_handle_input_range_	input_range;
		libbio::sam::header						header;
		libbio::sam::reader						reader;
		
		explicit alignment_input(libbio::file_handle &&fh):
			input_range(std::move(fh))
		{
			input_range.prepare();
		}
		
		static alignment_input open_path_or_stdin(char const *path);
		
		template <typename t_cb>
		void read_records(t_cb &&cb) { reader.read_records(header, input_range, std::forward <t_cb>(cb)); }
	};
}

#endif
