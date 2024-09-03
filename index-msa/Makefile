MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				index_handling.o \
				index_vector_builder.o \
				input_processor.o \
				main.o

all: index_msa

clean:
	$(RM) $(OBJECTS) index_msa cmdline.h cmdline.c

index_msa: $(OBJECTS) ../libpanvc3/libpanvc3.a ../$(LIBBIO_LIB)
	$(CXX) -o $@ $^ $(LDFLAGS) -pthread

main.cc: cmdline.c
cmdline.c: config.h

include $(MAKEFILE_DIR)/../make/config.mk
