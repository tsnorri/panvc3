/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/file_handling.hh>
#include <libbio/utility.hh>
#include "index_handling.hh"

namespace lb	= libbio;


namespace panvc3::msa_indices {

	void load_msa_index(char const *path, panvc3::msa_index &msa_index)
	{
		lb::log_time(std::cerr) << "Loading the input MSA index…\n";
		lb::file_istream stream;
		lb::open_file_for_reading(path, stream);
		cereal::PortableBinaryInputArchive archive(stream);
		archive(msa_index);
	}
}
