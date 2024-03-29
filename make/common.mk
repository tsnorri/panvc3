# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

# Default values.
WARNING_FLAGS_		?=
WARNING_CXXFLAGS_	?=
WARNING_FLAGS		?= -Wall -Werror -Wno-deprecated-declarations -Wno-unused $(WARNING_FLAGS_)
WARNING_CXXFLAGS	?= $(WARNING_CXXFLAGS_)
OPT_FLAGS			?= -O2 -ggdb

CMAKE			?= cmake
CP				?= cp
DOT				?= dot
GENGETOPT		?= gengetopt
MKDIR			?= mkdir
NINJA			?= ninja
PATCH			?= patch
RAGEL			?= ragel
TAR				?= tar
WGET			?= wget

CFLAGS			?=
CXXFLAGS		?=
CPPFLAGS		?=
LDFLAGS			?=
SYSTEM_CFLAGS	?=
SYSTEM_CXXFLAGS	?=
SYSTEM_CPPFLAGS	?=
SYSTEM_LDFLAGS	?=

LIBDISPATCH_LIBRARIES =

# Target type description, used currently in the .tar.gz name.
TARGET_TYPE		?=

CFLAGS			+= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CFLAGS)
CXXFLAGS		+= -std=c++2b $(OPT_FLAGS) $(WARNING_FLAGS) $(WARNING_CXXFLAGS) $(SYSTEM_CXXFLAGS)
CPPFLAGS		+=	-DHAVE_CONFIG_H \
					-I../include \
					-I../lib/cereal/include \
					-I../lib/libbio/include \
					-I../lib/libbio/lib/GSL/include \
					-I../lib/libbio/lib/range-v3/include \
					-I../lib/rapidcheck/include \
					-I../lib/rapidcheck/extras/catch/include \
					-I../lib/sdsl-lite/include \
					$(BOOST_INCLUDE) \
					$(SYSTEM_CPPFLAGS)
LDFLAGS			:= $(BOOST_LIBS) $(LDFLAGS) $(SYSTEM_LDFLAGS)

ifeq ($(shell uname -s),Linux)
	CPPFLAGS				+= -I../lib/swift-corelibs-libdispatch
	LIBDISPATCH_LIBRARIES	:= ../lib/swift-corelibs-libdispatch/build/src/libdispatch.a ../lib/swift-corelibs-libdispatch/build/src/BlocksRuntime/libBlocksRuntime.a
endif


# FIXME: the first two likely only work with Clang; I think GCC uses something else than -coverage.
%.cov.o: %.c
	$(CC) -c -coverage $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.cov.o: %.cc
	$(CXX) -c -coverage $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"

%.cc: %.rl
	$(RAGEL) -L -C -G2 -o $@ $<

%.dot: %.rl
	$(RAGEL) -V -p -o $@ $<

%.pdf: %.dot
	$(DOT) -Tpdf $< > $@
