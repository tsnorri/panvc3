ifeq ($(CLANG_VERSION),)
CLANG_SUFFIX :=
else
CLANG_SUFFIX := -$(CLANG_VERSION)
endif

CLANG				?= $(CLANG_ROOT)/bin/clang$(CLANG_SUFFIX)
CLANGXX				?= $(CLANG_ROOT)/bin/clang++$(CLANG_SUFFIX)
CC					= $(CLANG)
CXX					= $(CLANGXX)
WARNING_CXXFLAGS_	+= -Wno-redundant-consteval-if

BOOST_ROOT_LLVM		?= /usr
BOOST_ROOT			?= $(BOOST_ROOT_LLVM)
BOOST_LIBDIR		?= $(BOOST_ROOT)/lib
