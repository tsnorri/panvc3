CLANG_ROOT			= /opt/homebrew/opt/llvm
GCC_ROOT			= /opt/homebrew

CC					= $(CLANG_ROOT)/bin/clang
CXX					= $(CLANG_ROOT)/bin/clang++
OPT_FLAGS			= -O0 -g
#OPT_FLAGS			= -Og -g
#OPT_FLAGS			= -O2 -g
DOT					= /opt/homebrew/bin/dot
WGET				= /usr/local/homebrew/bin/wget

ISYSROOT			= -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk

# Make Boost not use std::unary_function and std::binary_function with BOOST_NO_CXX98_FUNCTION_BASE. (These have been deprecated.)
# Boost (as of 1.80.0) correctly defines the macro for libstdc++ but for some reason does not for libc++.
CPPFLAGS			= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_NO_CXX98_FUNCTION_BASE -DBOOST_STACKTRACE_USE_NOOP -I/usr/local/homebrew/include
LDFLAGS				= -L/usr/local/homebrew/lib -static -static-libgcc -stdlib=libc++ -lc++ -lpthread -lbsd -lz -ldl

SYSTEM_CFLAGS		= -mmacosx-version-min=10.11 $(ISYSROOT)
SYSTEM_CXXFLAGS		= -mmacosx-version-min=10.11 -faligned-allocation -stdlib=libc++ -nostdinc++ -I$(CLANG_ROOT)/include/c++/v1 $(ISYSROOT)
SYSTEM_LDFLAGS		= -mmacosx-version-min=10.11 -faligned-allocation -stdlib=libc++ -L$(CLANG_ROOT)/lib -Wl,-rpath,$(CLANG_ROOT)/lib

BOOST_ROOT			= /opt/homebrew/opt/boost
BOOST_LIBS			= -L$(BOOST_ROOT)/lib -lboost_iostreams-mt
