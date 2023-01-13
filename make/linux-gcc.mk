MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/gcc.mk

OPT_FLAGS			?= -O2 -g

CPPFLAGS			?= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH
WARNING_FLAGS_		= -Wno-deprecated-builtins
WARNING_CXXFLAGS_	= -Wno-interference-size

LDFLAGS				?= -lpthread -lrt
SYSTEM_CFLAGS		?= $(GCC_INCLUDES)
SYSTEM_CXXFLAGS		?= -nostdinc++ $(GCC_INCLUDES)
SYSTEM_LDFLAGS		?= -L$(GCC_ROOT)/lib -Wl,-rpath,$(GCC_ROOT)/lib

BOOST_LIBS			?= $(BOOST_ROOT)/lib/libboost_iostreams.a
