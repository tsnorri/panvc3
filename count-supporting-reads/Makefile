MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../local.mk
include  $(MAKEFILE_DIR)/../make/$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: count_supporting_reads

clean:
	$(RM) $(OBJECTS) count_supporting_reads cmdline.h cmdline.c

count_supporting_reads: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libpanvc3/libpanvc3.a ../lib/libbio/build-gcc/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
