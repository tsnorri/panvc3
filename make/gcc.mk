ifeq ($(GCC_VERSION),)
GCC_SUFFIX :=
else
GCC_SUFFIX := -$(GCC_VERSION)
endif

GCC					?= $(GCC_ROOT)/bin/gcc$(GCC_SUFFIX)
GXX					?= $(GCC_ROOT)/bin/g++$(GCC_SUFFIX)
CC					= $(GCC)
CXX					= $(GXX)

BOOST_ROOT_GCC		?= /usr
BOOST_ROOT			?= $(BOOST_ROOT_GCC)
