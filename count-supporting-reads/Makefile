-include ../local.gcc.mk
include ../common.mk
include ../seqan3.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: count_supporting_reads

clean:
	$(RM) $(OBJECTS) count_supporting_reads cmdline.h cmdline.c

count_supporting_reads: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/src/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include ../config.mk
