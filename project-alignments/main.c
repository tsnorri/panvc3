/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <stdlib.h> // EXIT_SUCCESS

// FIXME: This file is no longer needed; we can use the C++ compiler to compile main now that GCD is not used.

extern void panvc3_project_alignments(int argc, char **argv);


int main(int argc, char **argv)
{
	panvc3_project_alignments(argc, argv);
	return 0;
}
