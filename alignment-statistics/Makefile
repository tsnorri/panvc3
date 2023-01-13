MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../local.mk
include  $(MAKEFILE_DIR)/../make/$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

.PRECIOUS: 

.PHONY: all clean

OBJECTS		=	cmdline.o \
				main.o

all: alignment_statistics

clean:
	$(RM) $(OBJECTS) alignment_statistics cmdline.h cmdline.c

alignment_statistics: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/build-gcc/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
