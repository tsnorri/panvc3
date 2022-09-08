include local.llvm.mk
include common.mk


all: lib/swift-corelibs-libdispatch/build/src/libdispatch.a


lib/swift-corelibs-libdispatch/CMakeLists.txt.original:
	$(CP) lib/swift-corelibs-libdispatch/CMakeLists.txt lib/swift-corelibs-libdispatch/CMakeLists.txt.original
	$(PATCH) lib/swift-corelibs-libdispatch/CMakeLists.txt.original \
		tools/swift-cmakelists.patch \
		-o lib/swift-corelibs-libdispatch/CMakeLists.txt


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
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(NINJA) -v
