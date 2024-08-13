MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: process_alignments

clean:
	$(RM) $(OBJECTS) process_alignments cmdline.h cmdline.c

process_alignments: $(OBJECTS) ../libpanvc3/libpanvc3.a ../$(LIBBIO_LIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
