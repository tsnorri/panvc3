MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../make/local.$(OS_NAME)-llvm.mk
include  $(MAKEFILE_DIR)/../make/common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: convert_bed_positions

clean:
	$(RM) $(OBJECTS) convert_bed_positions cmdline.h cmdline.c

convert_bed_positions: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/build-llvm/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
