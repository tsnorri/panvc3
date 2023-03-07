MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../local.mk
include  $(MAKEFILE_DIR)/../make/$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

.PRECIOUS: 

.PHONY: all clean

OBJECTS		=	cmdline.o \
				log_memory_usage.o \
				main.o

CPPFLAGS += -DLIBBIO_LOG_ALLOCATED_MEMORY

all: project_alignments

clean:
	$(RM) $(OBJECTS) project_alignments cmdline.h cmdline.c

project_alignments: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/build-gcc/libbio.a ../libpanvc3/libpanvc3.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
