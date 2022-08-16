-include ../local.mk
include ../common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: make_indexable_fasta

clean:
	$(RM) $(OBJECTS) make_indexable_fasta cmdline.h cmdline.c

make_indexable_fasta: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/libbio/src/libbio.a

main.cc: cmdline.c
cmdline.c: config.h

include ../config.mk
