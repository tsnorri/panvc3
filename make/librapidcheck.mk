MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKE_SCRIPT_DIR)/common.mk


all: lib/libbio/lib/rapidcheck/build/librapidcheck.a


clean:
	$(RM) -rf lib/libbio/lib/rapidcheck/build


lib/libbio/lib/rapidcheck/build/librapidcheck.a:
	$(RM) -rf lib/libbio/lib/rapidcheck/build && \
	cd lib/libbio/lib/rapidcheck && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(MAKE)
