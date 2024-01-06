GCC_VERSION = 12
CLANG_ROOT = /usr/lib/llvm-16
LLVM_ROOT = /usr/lib/llvm-16
OPT_FLAGS = -O2 -g

BOOST_ROOT		= /usr
BOOST_INCLUDE	=
BOOST_LIBDIR	= /usr/lib/x86_64-linux-gnu

LIBDISPATCH_CFLAGS =
LIBDISPATCH_CXXFLAGS =
LIBDISPATCH_LDFLAGS = -pthread

LLVM_CFLAGS = -fblocks
LLVM_CXXFLAGS = -fblocks
LLVM_LIBGCC_LDFLAGS = -lgcc -lgcc_s
LLVM_LDFLAGS = -stdlib=libstdc++ -pthread -lbsd -lz -ldl -lm -lc -lrt $(LLVM_LIBGCC_LDFLAGS)

SYSTEM_CXXFLAGS =

LDFLAGS = -static -static-libgcc
