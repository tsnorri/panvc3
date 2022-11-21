MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../make/local.$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

CPPFLAGS += -DCATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS -I../lib/libbio/lib/Catch2/single_include

OBJECTS	=	main.o \
			cigar_eq.o \
			rewrite_cigar.o \
			rewrite_cigar_arbitrary.o

CATCH2_HEADERS	= $(shell find ../lib/libbio/lib/Catch2/include)

all: tests


tests: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../libpanvc3/libpanvc3.a ../lib/rapidcheck/build/librapidcheck.a $(BOOST_LIBS)


clean:
	$(RM) $(OBJECTS) tests


../lib/Catch2/single_include/catch2/catch.hpp: $(CATCH2_HEADERS)
	cd ../lib/libbio/lib/Catch2 && $(PYTHON) scripts/generateSingleHeader.py
