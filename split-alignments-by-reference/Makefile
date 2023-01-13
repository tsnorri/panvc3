MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../local.mk
include  $(MAKEFILE_DIR)/../make/$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: split_alignments_by_reference

clean:
	$(RM) $(OBJECTS) split_alignments_by_reference cmdline.h cmdline.c

split_alignments_by_reference: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/build-gcc/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
