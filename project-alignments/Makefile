MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

.PRECIOUS: 

.PHONY: all clean

OBJECTS		=	cmdline.o \
				log_memory_usage.o \
				project_alignments.o \
				main.o

#CPPFLAGS += -DLIBBIO_LOG_ALLOCATED_MEMORY

all: project_alignments

clean:
	$(RM) $(OBJECTS) project_alignments cmdline.h cmdline.c

project_alignments: $(OBJECTS) ../libpanvc3/libpanvc3.a ../$(LIBBIO_LIB)
	$(CXX) -o $@ $^ $(LDFLAGS) -pthread

main.c: cmdline.c
project_alignments.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
