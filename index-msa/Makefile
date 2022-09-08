-include ../local.llvm.mk
include ../common.mk

.PRECIOUS: 

OBJECTS		=	cmdline.o \
				main.o

all: index_msa

clean:
	$(RM) $(OBJECTS) index_msa cmdline.h cmdline.c

index_msa: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../lib/libbio/build-llvm/libbio.a $(LDFLAGS)

main.cc: cmdline.c
cmdline.c: config.h

include ../config.mk
