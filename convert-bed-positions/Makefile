MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: convert_bed_positions

clean:
	$(RM) $(OBJECTS) convert_bed_positions cmdline.h cmdline.c

convert_bed_positions: $(OBJECTS) ../$(LIBBIO_LIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
