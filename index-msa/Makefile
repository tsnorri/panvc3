-include ../local.llvm.mk
include ../common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: make_indexable_fasta

clean:
	$(RM) $(OBJECTS) make_indexable_fasta cmdline.h cmdline.c

make_indexable_fasta: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/build-llvm/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include ../config.mk
