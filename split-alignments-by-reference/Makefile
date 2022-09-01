-include ../local.gcc.mk
include ../common.mk
include ../seqan3.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: split_alignments_by_reference

clean:
	$(RM) $(OBJECTS) split_alignments_by_reference cmdline.h cmdline.c

split_alignments_by_reference: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/src/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include ../config.mk
