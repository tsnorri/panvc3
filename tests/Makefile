MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../make/local.$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

CPPFLAGS += -DCATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS -I../lib/libbio/lib/Catch2/single_include

OBJECTS	=	main.o \
			alignment_projector.o \
			cigar_eq.o \
			indel_run_checker.o \
			indel_run_checker_arbitrary.o \
			rewrite_cigar.o \
			rewrite_cigar_arbitrary.o

all: tests


tests: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../libpanvc3/libpanvc3.a ../lib/rapidcheck/build/librapidcheck.a ../lib/libbio/build-gcc/libbio.a $(BOOST_LIBS)


clean:
	$(RM) $(OBJECTS) tests
