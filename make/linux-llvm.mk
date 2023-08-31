# Link libc++ and libc++abi statically.
# Requires that CLANG_INCLUDE_DIR (e.g. LLVM_ROOT/lib/clang/15/include), LIBCXX_ROOT and BOOST_ROOT are defined in local.mk.

MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/llvm.mk

LIBDISPATCH_CFLAGS		?= -U__STDC_HOSTED__ $(CLANG_INCLUDES)
LIBDISPATCH_CXXFLAGS	?= -U__STDC_HOSTED__ $(CLANG_INCLUDES) -stdlib=libc++
LIBDISPATCH_LDFLAGS		?= -lpthread -lrt -lc++abi

WARNING_FLAGS_			= -Wno-deprecated-builtins
LLVM_CPPFLAGS			?=
LLVM_CFLAGS				?= -fblocks -U__STDC_HOSTED__ $(CLANG_INCLUDES)
LLVM_CXXFLAGS			?= -fblocks -nostdinc++ -U__STDC_HOSTED__ $(CLANG_INCLUDES)
LLVM_LIBGCC_LDFLAGS		?= -lgcc_s -lgcc
LLVM_LDFLAGS			?= -stdlib=libc++ -nodefaultlibs $(LIBCXX_ROOT)/lib/libc++.a $(LIBCXX_ROOT)/lib/libc++abi.a -lpthread -lbsd -lz -ldl -lm -lc -lrt $(LLVM_LIBGCC_LDFLAGS)
CPPFLAGS				?= $(LLVM_CPPFLAGS)
CFLAGS					?= $(LLVM_CFLAGS)
CXXFLAGS				?= $(LLVM_CXXFLAGS)
LDFLAGS					?= $(LLVM_LDFLAGS)

BOOST_LIBS				?= $(BOOST_ROOT)/lib/libboost_iostreams.a

#SYSTEM_CFLAGS += -fsanitize=thread
#SYSTEM_CXXFLAGS += -fsanitize=thread
#SYSTEM_LDFLAGS += -fsanitize=thread
#LIBDISPATCH_CFLAGS += -fsanitize=thread
#LIBDISPATCH_CXXFLAGS += -fsanitize=thread
#LIBDISPATCH_LDFLAGS += -fsanitize=thread
