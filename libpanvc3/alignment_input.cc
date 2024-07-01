/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/file_handling.hh>
#include <panvc3/alignment_input.hh>
#include <unistd.h>						// STDIN_FILENO

namespace lb	= libbio;


namespace panvc3 {
	
	alignment_input alignment_input::open_path_or_stdin(char const *path)
	{
		if (path)
			return alignment_input(lb::file_handle(lb::open_file_for_reading(path)));
		else
			return alignment_input(lb::file_handle(STDIN_FILENO, false));
	}
}
