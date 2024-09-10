MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKE_SCRIPT_DIR)/common.mk


all: lib/libbio/lib/Catch2/build/src/libCatch2.a


clean:
	$(RM) -rf lib/libbio/lib/Catch2/build


lib/libbio/lib/Catch2/build/src/libCatch2.a:
	$(RM) -rf lib/libbio/lib/Catch2/build && \
	cd lib/libbio/lib/Catch2 && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(MAKE)
