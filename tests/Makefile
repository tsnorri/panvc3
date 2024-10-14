MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

CPPFLAGS +=	-DCATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS \
			-I../lib/libbio/lib/Catch2/src \
			-I../lib/libbio/lib/Catch2/build/generated-includes \
			-I../lib/libbio/lib/rapidcheck/include \
			-I../lib/libbio/lib/rapidcheck/extras/catch/include

TEST_DEPENDENCIES =	../libpanvc3/libpanvc3.a \
					../lib/libbio/lib/rapidcheck/build/librapidcheck.a \
					../lib/libbio/lib/Catch2/build/src/libCatch2Main.a \
					../lib/libbio/lib/Catch2/build/src/libCatch2.a \
					../$(LIBBIO_LIB)

OBJECTS	=	alignment_projector.o \
			cigar_eq.o \
			indel_run_checker.o \
			indel_run_checker_arbitrary.o \
			rewrite_cigar.o \
			rewrite_cigar_arbitrary.o


all: tests


tests: $(OBJECTS) $(TEST_DEPENDENCIES)
	$(CXX) -o $@ $^ $(LDFLAGS) $(BOOST_LIBS)


clean:
	$(RM) $(OBJECTS) tests
