include local.mk
include common.mk

# Some of the tools require Clang for building while others require GCC.
# Unfortunately using either libstdc++ or libc++ with all of the tools
# was not possible, so we build two versions of libbio.
DEPENDENCIES =	lib/libbio/build-gcc/libbio.a \
				lib/libbio/build-llvm/libbio.a

OS_NAME = $(shell ./tools/os_name.sh)

ifeq ($(OS_NAME),Linux)
	DEPENDENCIES += lib/swift-corelibs-libdispatch/build/src/libdispatch.a
endif

# “$() $()” is a literal space.
VERSION = $(subst $() $(),-,$(shell tools/git_version.sh))
DIST_TARGET_DIR = panvc3-$(VERSION)
DIST_NAME_SUFFIX = $(if $(TARGET_TYPE),-$(TARGET_TYPE),)
DIST_TAR_GZ = panvc3-$(VERSION)-$(OS_NAME)$(DIST_NAME_SUFFIX).tar.gz


.PHONY: all clean-all clean clean-dependencies dependencies

all:	index-msa/index_msa \
		count-supporting-reads/count_supporting_reads \
		split-alignments-by-reference/split_alignments_by_reference

clean-all: clean clean-dependencies clean-dist

clean:
	$(MAKE) -C index-msa clean
	$(MAKE) -C count-supporting-reads clean
	$(MAKE) -C split-alignments-by-reference clean

clean-dependencies:
	$(RM) -r lib/libbio/build-gcc lib/libbio/build-llvm
	$(RM) -r lib/swift-corelibs-libdispatch/build

clean-dist:
	$(RM) -rf $(DIST_TARGET_DIR)

dependencies: $(DEPENDENCIES)

dist: $(DIST_TAR_GZ)

test:
	$(MAKE) -C tests

index-msa/index_msa: lib/libbio/build-llvm/libbio.a
	$(MAKE) -C index-msa

count-supporting-reads/count_supporting_reads: lib/libbio/build-gcc/libbio.a
	$(MAKE) -C count-supporting-reads

split-alignments-by-reference/split_alignments_by_reference: lib/libbio/build-gcc/libbio.a
	$(MAKE) -C split-alignments-by-reference

$(DIST_TAR_GZ):	index-msa/index_msa
	$(MKDIR) -p $(DIST_TARGET_DIR)
	$(CP) index-msa/index_msa $(DIST_TARGET_DIR)
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/swift-corelibs-libdispatch/LICENSE $(DIST_TARGET_DIR)/swift-corelibs-libdispatch-license.txt
	$(CP) lib/cereal/LICENSE $(DIST_TARGET_DIR)/cereal-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

lib/libbio/build-%/libbio.a:
	$(MKDIR) -p lib/libbio/build-$*
	VPATH=../src $(MAKE) -C lib/libbio/build-$* -f ../../../local.$(OS_NAME)-$*.mk -f ../src/Makefile

lib/swift-corelibs-libdispatch/build/src/libdispatch.a:
	$(MAKE) -f libdispatch.mk
