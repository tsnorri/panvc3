ifeq ($(GCC_VERSION),)
GCC_SUFFIX :=
else
GCC_SUFFIX := -$(GCC_VERSION)
endif

GCC_ROOT			?= /usr
GCC					?= $(GCC_ROOT)/bin/gcc$(GCC_SUFFIX)
GXX					?= $(GCC_ROOT)/bin/g++$(GCC_SUFFIX)

ifeq "$(origin CC)" "default"
CC					= $(GCC)
endif

ifeq "$(origin CXX)" "default"
CXX					= $(GXX)
endif

BOOST_ROOT_GCC		?= /usr
BOOST_ROOT			?= $(BOOST_ROOT_GCC)
BOOST_INCLUDE_GCC	?= -I$(BOOST_ROOT)/include
BOOST_INCLUDE		?= $(BOOST_INCLUDE_GCC)
BOOST_LIBDIR_GCC	?= $(BOOST_ROOT)/lib
BOOST_LIBDIR		?= $(BOOST_LIBDIR_GCC)
