/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef PANVC3_MSA_INDICES_INDEX_HANDLING_HH
#define PANVC3_MSA_INDICES_INDEX_HANDLING_HH

#include <panvc3/msa_index.hh>


namespace panvc3::msa_indices {
	void load_msa_index(char const *path, panvc3::msa_index &msa_index);
}

#endif
