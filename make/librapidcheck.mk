MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKE_SCRIPT_DIR)/common.mk


all: lib/rapidcheck/build/librapidcheck.a


clean:
	$(RM) -rf lib/rapidcheck/build


lib/rapidcheck/build/librapidcheck.a:
	$(RM) -rf lib/rapidcheck/build && \
	cd lib/rapidcheck && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DCMAKE_CXX_FLAGS="$(CXXFLAGS)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	VERBOSE=1 $(MAKE)
