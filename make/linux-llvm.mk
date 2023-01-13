# Link libc++ statically.
# Requires that CLANG_INCLUDE_DIR (e.g. LLVM_ROOT/lib/clang/15/include), LIBCXX_ROOT and BOOST_ROOT are defined in local.mk.

MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/llvm.mk

LIBDISPATCH_CFLAGS		?= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) 
LIBDISPATCH_CXXFLAGS	?= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) -stdlib=libc++

CPPFLAGS				?= -DBOOST_NO_CXX98_FUNCTION_BASE
CFLAGS					?= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR)
CXXFLAGS				?= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) -stdlib=libc++
WARNING_FLAGS_			= -Wno-deprecated-builtins
LDFLAGS					?= -stdlib=libc++ -nodefaultlibs $(LIBCXX_ROOT)/libc++.a $(LIBCXX_ROOT)/libc++abi.a -lpthread -lbsd -lz -ldl -lm -lc -lgcc_s -lgcc -lrt

BOOST_LIBS				?= -L$(BOOST_ROOT)/lib -lboost_iostreams
