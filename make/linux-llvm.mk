# Link libc++ and libc++abi statically.
# Requires that CLANG_INCLUDE_DIR (e.g. LLVM_ROOT/lib/clang/15/include), LIBCXX_ROOT and BOOST_ROOT are defined in local.mk.

MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/llvm.mk

LIBDISPATCH_CFLAGS		?= -U__STDC_HOSTED__ $(CLANG_INCLUDES)
LIBDISPATCH_CXXFLAGS	?= -U__STDC_HOSTED__ $(CLANG_INCLUDES) -stdlib=libc++
LIBDISPATCH_LDFLAGS		?= -lpthread -lrt

CPPFLAGS				?= -DBOOST_NO_CXX98_FUNCTION_BASE
CFLAGS					?= -fblocks -U__STDC_HOSTED__ $(CLANG_INCLUDES)
CXXFLAGS				?= -fblocks -nostdinc++ -U__STDC_HOSTED__ $(CLANG_INCLUDES)
WARNING_FLAGS_			= -Wno-deprecated-builtins
LDFLAGS					?= -stdlib=libc++ -nodefaultlibs $(LIBCXX_ROOT)/lib/libc++.a $(LIBCXX_ROOT)/lib/libc++abi.a -lpthread -lbsd -lz -ldl -lm -lc -lgcc_s -lgcc -lrt

BOOST_LIBS				?= $(BOOST_ROOT)/lib/libboost_iostreams.a
