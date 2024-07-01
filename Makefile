PANVC3_PROJECT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(PANVC3_PROJECT_DIR)/make/os-name.mk
-include $(PANVC3_PROJECT_DIR)/local.mk
include  $(PANVC3_PROJECT_DIR)/make/common.mk

# Some of the tools require Clang for building while others require GCC.
# Unfortunately using either libstdc++ or libc++ with all of the tools
# was not possible, so we build two versions of libbio.
DEPENDENCIES =	lib/libbio/build-gcc/libbio.a \
				lib/libbio/build-llvm/libbio.a \
				lib/rapidcheck/build/librapidcheck.a

LIBDISPATCH_LINUX_BIN :=
ifeq ($(OS_NAME),linux)
	LIBDISPATCH_LINUX_BIN := lib/swift-corelibs-libdispatch/build/src/libdispatch.a
	DEPENDENCIES += $(LIBDISPATCH_LINUX_BIN)
endif

BUILD_PRODUCTS	=	alignment-statistics/alignment_statistics \
					convert-bed-positions/convert_bed_positions \
					count-supporting-reads/count_supporting_reads \
					index-msa/index_msa \
					libpanvc3/libpanvc3.a \
					process-alignments/process_alignments \
					project-alignments/project_alignments \
					recalculate-mapq/recalculate_mapq \
					rewrite-cigar/rewrite_cigar \
					split-alignments-by-reference/split_alignments_by_reference \
					subset-alignments/subset_alignments \

# “$() $()” is a literal space.
VERSION = $(subst $() $(),-,$(shell tools/git_version.sh))
DIST_TARGET_DIR = panvc3-$(VERSION)
DIST_NAME_SUFFIX = $(if $(TARGET_TYPE),-$(TARGET_TYPE),)
DIST_TAR_GZ = panvc3-$(VERSION)-$(OS_NAME)$(DIST_NAME_SUFFIX).tar.gz

CATCH2_HEADERS	= $(shell find lib/libbio/lib/Catch2/include)


.PHONY: all clean-all clean clean-dependencies dependencies libbio-tests
.PRECIOUS: lib/libbio/lib/rapidcheck/build-gcc/librapidcheck.a lib/libbio/lib/rapidcheck/build-llvm/librapidcheck.a


all: $(BUILD_PRODUCTS)

clean-all: clean clean-dependencies clean-dist

clean:
	$(MAKE) -C alignment-statistics clean
	$(MAKE) -C convert-bed-positions clean
	$(MAKE) -C count-supporting-reads clean
	$(MAKE) -C index-msa clean
	$(MAKE) -C libpanvc3 clean
	$(MAKE) -C project-alignments clean
	$(MAKE) -C recalculate-mapq clean
	$(MAKE) -C rewrite-cigar clean
	$(MAKE) -C split-alignments-by-reference clean
	$(MAKE) -C subset-alignments clean
	$(MAKE) -C tests clean

clean-dependencies:
	$(RM) -r lib/libbio/build-gcc lib/libbio/build-llvm
	$(RM) -r lib/swift-corelibs-libdispatch/build
	$(RM) -r lib/rapidcheck/build

clean-dist:
	$(RM) -rf $(DIST_TARGET_DIR)

dependencies: $(DEPENDENCIES)

dist: $(DIST_TAR_GZ)

tests: lib/rapidcheck/build/librapidcheck.a libpanvc3/libpanvc3.a lib/Catch2/single_include/catch2/catch.hpp
	$(MAKE) -C tests

libbio-tests: lib/libbio/build-tests-gcc/tests lib/libbio/build-tests-llvm/tests

alignment-statistics/alignment_statistics: lib/libbio/build-gcc/libbio.a
	$(MAKE) -C alignment-statistics

convert-bed-positions/convert_bed_positions: lib/libbio/build-llvm/libbio.a
	$(MAKE) -C convert-bed-positions
	
count-supporting-reads/count_supporting_reads: lib/libbio/build-gcc/libbio.a
	$(MAKE) -C count-supporting-reads

index-msa/index_msa: lib/libbio/build-llvm/libbio.a $(LIBDISPATCH_LINUX_BIN)
	$(MAKE) -C index-msa

libpanvc3/libpanvc3.a:
	$(MAKE) -C libpanvc3

process-alignments/process_alignments: lib/libbio/build-gcc/libbio.a libpanvc3/libpanvc3.a
	$(MAKE) -C process-alignments

project-alignments/project_alignments: lib/libbio/build-gcc/libbio.a libpanvc3/libpanvc3.a
	$(MAKE) -C project-alignments

recalculate-mapq/recalculate_mapq: lib/libbio/build-gcc/libbio.a libpanvc3/libpanvc3.a
	$(MAKE) -C recalculate-mapq

rewrite-cigar/rewrite_cigar: lib/libbio/build-gcc/libbio.a libpanvc3/libpanvc3.a
	$(MAKE) -C rewrite-cigar

split-alignments-by-reference/split_alignments_by_reference: lib/libbio/build-gcc/libbio.a
	$(MAKE) -C split-alignments-by-reference

subset-alignments/subset_alignments: lib/libbio/build-gcc/libbio.a
	$(MAKE) -C subset-alignments

$(DIST_TAR_GZ):	$(BUILD_PRODUCTS)
	$(MKDIR) -p $(DIST_TARGET_DIR)
	for f in $(BUILD_PRODUCTS); do if [ $${f##*/} = libpanvc3.a ]; then $(CP) $$f $(DIST_TARGET_DIR); else $(CP) $$f $(DIST_TARGET_DIR)/panvc3_$${f##*/}; fi; done
	$(CP) count-supporting-reads/calculate_reference_bias.py $(DIST_TARGET_DIR)/panvc3_calculate_reference_bias.py
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/swift-corelibs-libdispatch/LICENSE $(DIST_TARGET_DIR)/swift-corelibs-libdispatch-license.txt
	$(CP) lib/cereal/LICENSE $(DIST_TARGET_DIR)/cereal-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

lib/libbio/build-%/libbio.a:
	$(MKDIR) -p lib/libbio/build-$*
	VPATH=../src $(MAKE) -C lib/libbio/build-$* -f ../../../local.mk -f ../../../make/$*.mk -f ../../../make/$(OS_NAME)-$*.mk -f ../src/Makefile

lib/libbio/build-tests-%/tests: lib/libbio/build-%/libbio.a lib/libbio/lib/rapidcheck/build-%/librapidcheck.a
	$(MKDIR) -p lib/libbio/build-tests-$*
	VPATH=../tests LIBBIO_PATH=../build-$*/libbio.a RAPIDCHECK_BUILD_DIR="build-$*" $(MAKE) -C lib/libbio/build-tests-$* -f ../../../local.mk -f ../../../make/$*.mk -f ../../../make/$(OS_NAME)-$*.mk -f ../tests/Makefile tests

lib/libbio/lib/rapidcheck/build-%/librapidcheck.a:
	RAPIDCHECK_BUILD_DIR="build-$*" $(MAKE) -C lib/libbio -f ../../local.mk -f ../../make/$*.mk -f ../../make/$(OS_NAME)-$*.mk -f Makefile lib/rapidcheck/build-$*/librapidcheck.a

lib/rapidcheck/build/librapidcheck.a:
	$(MAKE) -f make/librapidcheck.mk

lib/swift-corelibs-libdispatch/build/src/libdispatch.a:
	$(MAKE) -f make/libdispatch.mk

lib/Catch2/single_include/catch2/catch.hpp: $(CATCH2_HEADERS)
	cd lib/libbio/lib/Catch2 && $(PYTHON) scripts/generateSingleHeader.py
