MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../local.mk
include  $(MAKEFILE_DIR)/../make/$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

OBJECTS		=	alignment_projector.o \
				cigar.o \
				compressed_fasta_reader.o \
				indel_run_checker.o \
				sequence_buffer_store.o

CFLAGS		+= -fPIC
CXXFLAGS	+= -fPIC

all: libpanvc3.a

clean:
	$(RM) $(OBJECTS) libpanvc3.a

libpanvc3.a: $(OBJECTS)
	$(AR) rcs $@ $^
