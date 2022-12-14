# Build a static binary.

LLVM_ROOT				= /usr/lib/llvm-15
CLANG_INCLUDE_DIR		= $(LLVM_ROOT)/lib/clang/15/include

CC							= clang-15
CXX							= clang++-15
LIBDISPATCH_CFLAGS			= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) 
LIBDISPATCH_CXXFLAGS		= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) -stdlib=libc++

CPPFLAGS				= -DBOOST_NO_CXX98_FUNCTION_BASE
CFLAGS					= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR)
CXXFLAGS				= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) -stdlib=libc++
WARNING_FLAGS_		= -Wno-deprecated-builtins
LDFLAGS					= -static -static-libgcc -stdlib=libc++ -lc++ -lpthread -lbsd -lz -ldl

BOOST_ROOT				= /home/tnorri/local/boost-1.80.0-clang-15
BOOST_LIBS				= -L$(BOOST_ROOT)/lib -lboost_iostreams

# Used in .tar.gz name.
TARGET_TYPE				= static
