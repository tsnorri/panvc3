MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/common.mk

OBJECTS		=	alignment_input.o \
				cigar.o \
				sequence_buffer_store.o \
				utility.o

ifeq ($(shell uname -s),Linux)
	OBJECTS	+=	dispatch_event_linux.o
else
	OBJECTS	+=	dispatch_event_kqueue.o
endif

CFLAGS		+= -fPIC
CXXFLAGS	+= -fPIC

all: libpanvc3.a

clean:
	$(RM) $(OBJECTS) libpanvc3.a

libpanvc3.a: $(OBJECTS)
	$(AR) rcs $@ $^
