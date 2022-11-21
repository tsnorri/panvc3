MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKE_SCRIPT_DIR)/os-name.mk
-include $(MAKE_SCRIPT_DIR)/local.$(OS_NAME)-gcc.mk
include  $(MAKE_SCRIPT_DIR)/common.mk


all: lib/rapidcheck/build/librapidcheck.a


lib/rapidcheck/build/librapidcheck.a:
	$(RM) -rf lib/rapidcheck/build && \
	cd lib/rapidcheck && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(MAKE)
