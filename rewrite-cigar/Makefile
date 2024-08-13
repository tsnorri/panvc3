MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

.PRECIOUS: 

.PHONY: all clean

OBJECTS		=	cmdline.o \
				main.o

#CPPFLAGS += -DLIBBIO_LOG_ALLOCATED_MEMORY

all: rewrite_cigar

clean:
	$(RM) $(OBJECTS) rewrite_cigar cmdline.h cmdline.c

rewrite_cigar: $(OBJECTS) ../libpanvc3/libpanvc3.a ../$(LIBBIO_LIB)
	$(CXX) -o $@ $^ $(LDFLAGS) -pthread

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
