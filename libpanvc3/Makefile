MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../make/local.$(OS_NAME)-gcc.mk
include  $(MAKEFILE_DIR)/../make/common.mk
include  $(MAKEFILE_DIR)/../make/seqan3.mk

OBJECTS		=	cigar.o \
				indel_run_checker.o

CFLAGS		+= -fPIC
CXXFLAGS	+= -fPIC

all: libpanvc3.a

clean:
	$(RM) $(OBJECTS) libpanvc3.a

libpanvc3.a: $(OBJECTS)
	$(AR) rcs $@ $^
