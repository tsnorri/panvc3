MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKEFILE_DIR)/../make/os-name.mk
-include $(MAKEFILE_DIR)/../make/local.$(OS_NAME)-llvm.mk
include  $(MAKEFILE_DIR)/../make/common.mk

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

include $(MAKEFILE_DIR)/../make/config.mk
