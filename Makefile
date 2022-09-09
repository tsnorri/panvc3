PANVC3_PROJECT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(PANVC3_PROJECT_DIR)/make/os-name.mk
-include $(PANVC3_PROJECT_DIR)/make/local.mk
include  $(PANVC3_PROJECT_DIR)/make/common.mk

# Some of the tools require Clang for building while others require GCC.
# Unfortunately using either libstdc++ or libc++ with all of the tools
# was not possible, so we build two versions of libbio.
DEPENDENCIES =	lib/libbio/build-gcc/libbio.a \
				lib/libbio/build-llvm/libbio.a

ifeq ($(OS_NAME),linux)
	DEPENDENCIES += lib/swift-corelibs-libdispatch/build/src/libdispatch.a
endif

BUILD_PRODUCTS	=	index-msa/index_msa \
					count-supporting-reads/count_supporting_reads \
					convert-bed-positions/convert_bed_positions \
					split-alignments-by-reference/split_alignments_by_reference

# “$() $()” is a literal space.
VERSION = $(subst $() $(),-,$(shell tools/git_version.sh))
DIST_TARGET_DIR = panvc3-$(VERSION)
DIST_NAME_SUFFIX = $(if $(TARGET_TYPE),-$(TARGET_TYPE),)
DIST_TAR_GZ = panvc3-$(VERSION)-$(OS_NAME)$(DIST_NAME_SUFFIX).tar.gz


.PHONY: all clean-all clean clean-dependencies dependencies

all: $(BUILD_PRODUCTS)

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

convert-bed-positions/convert_bed_positions: lib/libbio/build-llvm/libbio.a
	$(MAKE) -C convert-bed-positions

split-alignments-by-reference/split_alignments_by_reference: lib/libbio/build-gcc/libbio.a
	$(MAKE) -C split-alignments-by-reference

$(DIST_TAR_GZ):	$(BUILD_PRODUCTS)
	$(MKDIR) -p $(DIST_TARGET_DIR)
	for f in $(BUILD_PRODUCTS); do $(CP) $$f $(DIST_TARGET_DIR); done
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/swift-corelibs-libdispatch/LICENSE $(DIST_TARGET_DIR)/swift-corelibs-libdispatch-license.txt
	$(CP) lib/cereal/LICENSE $(DIST_TARGET_DIR)/cereal-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

lib/libbio/build-%/libbio.a:
	$(MKDIR) -p lib/libbio/build-$*
	VPATH=../src $(MAKE) -C lib/libbio/build-$* -f ../../../make/local.$(OS_NAME)-$*.mk -f ../src/Makefile

lib/swift-corelibs-libdispatch/build/src/libdispatch.a:
	$(MAKE) -f make/libdispatch.mk
