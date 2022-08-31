GCC_ROOT			= /home/tnorri/.linuxbrew

CC					= $(GCC_ROOT)/bin/gcc-12
CXX					= $(GCC_ROOT)/bin/g++-12
OPT_FLAGS			= -O0 -g
#OPT_FLAGS			= -Og -g
#OPT_FLAGS			= -O2 -g
DOT					= /opt/homebrew/bin/dot
WGET				= /usr/local/homebrew/bin/wget

ISYSROOT			=

CPPFLAGS			= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH -I/usr/local/homebrew/include
LDFLAGS				= -L/home/tnorri/.linuxbrew/lib -static -static-libgcc -lpthread -lz -ldl

SYSTEM_CFLAGS		= $(ISYSROOT)
SYSTEM_CXXFLAGS		= -nostdinc++ -I$(GCC_ROOT)/include/c++/12 -I$(GCC_ROOT)/include/c++/12/x86_64-pc-linux-gnu $(ISYSROOT)
SYSTEM_LDFLAGS		= -L$(GCC_ROOT)/lib -Wl,-rpath,$(GCC_ROOT)/lib

#BOOST_ROOT			= /opt/homebrew/opt/boost
BOOST_ROOT			= /home/tnorri/local/boost-1.80.0-gcc-12
BOOST_LIBS			= -L$(BOOST_ROOT)/lib -lboost_iostreams
