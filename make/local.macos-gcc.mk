MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/local-macos.mk

GCC_ROOT			= /opt/homebrew/opt/gcc

CC					= $(GCC_ROOT)/bin/gcc-12
CXX					= $(GCC_ROOT)/bin/g++-12
OPT_FLAGS			= -O0 -gdwarf-3
#OPT_FLAGS			= -Og -gdwarf-5
#OPT_FLAGS			= -O2 -gdwarf-5

ISYSROOT			= -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk

CPPFLAGS			= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH -I/usr/local/homebrew/include
LDFLAGS				= -L/opt/homebrew/lib -lz
WARNING_FLAGS_		= -Wno-deprecated-builtins
WARNING_CXXFLAGS_	= -Wno-interference-size

SYSTEM_CFLAGS		= $(ISYSROOT)
SYSTEM_CXXFLAGS		= -nostdinc++ -I$(GCC_ROOT)/include/c++/12 -I$(GCC_ROOT)/include/c++/12/aarch64-apple-darwin21 $(ISYSROOT)
SYSTEM_LDFLAGS		= -L$(GCC_ROOT)/lib -Wl,-rpath,$(GCC_ROOT)/lib

BOOST_ROOT			= /Users/tsnorri/boost-1.80.0-gcc-12
BOOST_LIBS			= -L$(BOOST_ROOT)/lib -lboost_iostreams
