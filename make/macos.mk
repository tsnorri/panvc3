# Requires that MACOS_VERSION_MIN, GCC_INCLUDES, GCC_ROOT and BOOST_ROOT are defined in local.mk

OPT_FLAGS			?= -O2 -gdwarf-3
ISYSROOT			?= -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk

CPPFLAGS			?= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH $(ISYSROOT)
CFLAGS				?= -mmacosx-version-min=$(MACOS_VERSION_MIN) $(CPPFLAGS)
CXXFLAGS			?= -mmacosx-version-min=$(MACOS_VERSION_MIN) $(CPPFLAGS)
LDFLAGS				?= -mmacosx-version-min=$(MACOS_VERSION_MIN) -L$(GCC_ROOT)/lib -Wl,-rpath,$(GCC_ROOT)/lib
WARNING_FLAGS_		+= -Wno-deprecated-builtins
WARNING_CXXFLAGS_	+= -Wno-interference-size

BOOST_LIBS			?= -L$(BOOST_LIBDIR) -lboost_iostreams
