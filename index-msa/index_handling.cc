/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <iostream>
#include <libbio/file_handling.hh>
#include <libbio/utility.hh>
#include <panvc3/msa_index.hh>
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
