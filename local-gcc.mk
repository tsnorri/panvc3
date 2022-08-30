GCC_ROOT			= /opt/homebrew/opt/gcc

CC					= $(GCC_ROOT)/bin/gcc-12
CXX					= $(GCC_ROOT)/bin/g++-12
OPT_FLAGS			= -O0 -g
#OPT_FLAGS			= -Og -g
#OPT_FLAGS			= -O2 -g
DOT					= /opt/homebrew/bin/dot
WGET				= /usr/local/homebrew/bin/wget

ISYSROOT			= -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk

CPPFLAGS			= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH -I/usr/local/homebrew/include
LDFLAGS				= -L/usr/local/homebrew/lib
SYSTEM_CFLAGS		= -mmacosx-version-min=10.11 $(ISYSROOT)
SYSTEM_CXXFLAGS		= -mmacosx-version-min=10.11 -nostdinc++ -I$(GCC_ROOT)/include/c++/12 -I$(GCC_ROOT)/include/c++/12/aarch64-apple-darwin21 $(ISYSROOT)
SYSTEM_LDFLAGS		= -mmacosx-version-min=10.11 -L$(GCC_ROOT)/lib -Wl,-rpath,$(GCC_ROOT)/lib

#BOOST_ROOT			= /opt/homebrew/opt/boost
BOOST_ROOT			= /Users/tsnorri/boost-1.80.0-gcc-12
BOOST_LIBS			= -L$(BOOST_ROOT)/lib -lboost_iostreams
