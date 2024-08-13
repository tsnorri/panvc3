MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/gcc.mk

OPT_FLAGS			?= -O2 -ggdb
#OPT_FLAGS			?= -O0 -ggdb

GCC_CPPFLAGS		?= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH -I/usr/include/kqueue
CPPFLAGS			?= $(GCC_CPPFLAGS)
WARNING_FLAGS_		= -Wno-deprecated-builtins
WARNING_CXXFLAGS_	= -Wno-interference-size

LDFLAGS				?= -lpthread -lrt $(GCC_LDFLAGS)
SYSTEM_CFLAGS		?= $(GCC_INCLUDES)
SYSTEM_CXXFLAGS		?= -nostdinc++ $(GCC_INCLUDES)
SYSTEM_LDFLAGS		?= -L$(GCC_ROOT)/lib -Wl,-rpath,$(GCC_ROOT)/lib

#BOOST_LIBS			?= $(BOOST_LIBDIR)/libboost_iostreams.a
BOOST_LIBS			?= -lboost_iostreams

#CFLAGS += -fsanitize=thread -fsanitize=undefined
#CXXFLAGS += -fsanitize=thread -fsanitize=undefined
#LDFLAGS += -fsanitize=thread -fsanitize=undefined
#SYSTEM_CFLAGS += -fsanitize=thread
#SYSTEM_CXXFLAGS += -fsanitize=thread
#SYSTEM_CXXFLAGS += -fsanitize=thread
#SYSTEM_LDFLAGS += -fsanitize=thread
