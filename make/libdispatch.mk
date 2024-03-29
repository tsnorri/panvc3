MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(MAKE_SCRIPT_DIR)/../local.mk
include  $(MAKE_SCRIPT_DIR)/linux-llvm.mk
include  $(MAKE_SCRIPT_DIR)/common.mk


all: lib/swift-corelibs-libdispatch/build/src/libdispatch.a


lib/swift-corelibs-libdispatch/CMakeLists.txt.original:
	$(CP) lib/swift-corelibs-libdispatch/CMakeLists.txt lib/swift-corelibs-libdispatch/CMakeLists.txt.original
	$(PATCH) lib/swift-corelibs-libdispatch/CMakeLists.txt.original \
		tools/swift-cmakelists.patch \
		-o lib/swift-corelibs-libdispatch/CMakeLists.txt


# Linker flags below needed for building a test program.
lib/swift-corelibs-libdispatch/build/src/libdispatch.a: lib/swift-corelibs-libdispatch/CMakeLists.txt.original
	$(RM) -rf lib/swift-corelibs-libdispatch/build && \
	cd lib/swift-corelibs-libdispatch && \
	$(MKDIR) build && \
	cd build && \
	$(CP) ../../../tools/disable_warnings.cmake ../cmake/modules/V2MDisableCompilerWarnings.cmake && \
	$(CMAKE) \
		-G Ninja \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DCMAKE_C_FLAGS="$(LIBDISPATCH_CFLAGS)" \
		-DCMAKE_CXX_FLAGS="$(LIBDISPATCH_CXXFLAGS)" \
		-DCMAKE_EXE_LINKER_FLAGS="$(LIBDISPATCH_LDFLAGS)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(NINJA) -v
