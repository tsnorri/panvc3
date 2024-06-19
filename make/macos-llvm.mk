# Requires that MACOS_VERSION_MIN, CLANG_ROOT, LIBCXX_ROOT and BOOST_ROOT are defined in local.mk.

MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/llvm.mk

ISYSROOT			?= -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk

# Make Boost not use std::unary_function and std::binary_function with BOOST_NO_CXX98_FUNCTION_BASE. (These have been deprecated.)
# Boost (as of 1.80.0) correctly defines the macro for libstdc++ but for some reason does not for libc++.
CPPFLAGS			?= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP
LDFLAGS				?= -stdlib=libc++ -lz

SYSTEM_CFLAGS		?= -mmacosx-version-min=$(MACOS_VERSION_MIN) $(ISYSROOT)
SYSTEM_CXXFLAGS		?= -mmacosx-version-min=$(MACOS_VERSION_MIN) -faligned-allocation -stdlib=libc++ -nostdinc++ -I$(CLANG_ROOT)/include/c++/v1 $(ISYSROOT) -fexperimental-library
SYSTEM_LDFLAGS		?= -mmacosx-version-min=$(MACOS_VERSION_MIN) -faligned-allocation -stdlib=libc++ -L$(CLANG_ROOT)/lib -L$(CLANG_ROOT)/lib/c++ -Wl,-rpath,$(CLANG_ROOT)/lib -Wl,-rpath,$(CLANG_ROOT)/lib/c++

BOOST_LIBS			?= -L$(BOOST_LIBDIR) -lboost_iostreams-mt
