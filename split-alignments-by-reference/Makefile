MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: split_alignments_by_reference

clean:
	$(RM) $(OBJECTS) split_alignments_by_reference cmdline.h cmdline.c

split_alignments_by_reference: $(OBJECTS) ../libpanvc3/libpanvc3.a ../$(LIBBIO_LIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
