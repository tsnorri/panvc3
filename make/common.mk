MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKE_SCRIPT_DIR)/os-name.mk
-include $(MAKE_SCRIPT_DIR)/../local.mk
include  $(MAKE_SCRIPT_DIR)/$(OS_NAME).mk
include  $(MAKE_SCRIPT_DIR)/gcc.mk

# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

LIBBIO_LIB =	lib/libbio/src/libbio.a

# Default values.
WARNING_FLAGS_		?=
WARNING_CXXFLAGS_	?=
WARNING_FLAGS		?= -Wall -Werror -Wno-deprecated-declarations -Wno-unknown-pragmas -Wno-unused $(WARNING_FLAGS_)
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

# Target type description, used currently in the .tar.gz name.
TARGET_TYPE		?=

CFLAGS			+= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS)
CXXFLAGS		+= -std=c++2b $(OPT_FLAGS) $(WARNING_FLAGS) $(WARNING_CXXFLAGS)
CPPFLAGS		+=	-DHAVE_CONFIG_H \
					-I../include \
					-I../lib/cereal/include \
					-I../lib/libbio/include \
					-I../lib/libbio/lib/GSL/include \
					-I../lib/libbio/lib/range-v3/include \
					-I../lib/rapidcheck/include \
					-I../lib/rapidcheck/extras/catch/include \
					-I../lib/sdsl-lite/include \
					-I../lib/seqan3/include \
					$(BOOST_INCLUDE)
LDFLAGS			:=	$(BOOST_LIBS) -lbacktrace -lz $(LDFLAGS)

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
