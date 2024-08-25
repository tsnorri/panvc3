MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKE_SCRIPT_DIR)/common.mk


all: lib/libkqueue/build/libkqueue.a


clean:
	$(RM) -rf lib/libkqueue/build


lib/libkqueue/build/libkqueue.a:
	$(RM) -rf lib/libkqueue/build && \
	cd lib/libkqueue && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(MAKE)
