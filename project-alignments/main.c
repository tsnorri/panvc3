/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/dispatch/dispatch_compat.h>
#include <stdlib.h> // EXIT_SUCCESS

// The compiler (as of GCC 12) seems to optimise process() in project_alignments.cc
// in some way that interferes with the call to dispatch_main(). This results in
// "terminate called without an active exception" on Linux but not on macOS at least
// with some inputs. The solution is to have process() in another compilation unit;
// for extra safety main() is now in pure C.

extern void panvc3_project_alignments(int argc, char **argv);


int main(int argc, char **argv)
{
	panvc3_project_alignments(argc, argv);
	dispatch_main();
	
	// Not reached.
	return EXIT_SUCCESS;
}
