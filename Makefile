PANVC3_PROJECT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(PANVC3_PROJECT_DIR)/make/os-name.mk
-include $(PANVC3_PROJECT_DIR)/local.mk
include  $(PANVC3_PROJECT_DIR)/make/common.mk

DEPENDENCIES =	$(LIBBIO_LIB) \
				lib/libbio/lib/rapidcheck/build/librapidcheck.a

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


.PHONY: all clean-all clean clean-dependencies clean-dist dependencies tests libpanvc3 libbio libbio-tests
#.PRECIOUS: lib/libbio/lib/rapidcheck/build/librapidcheck.a


all: $(BUILD_PRODUCTS)

clean-all: clean clean-dependencies clean-dist

clean:
	$(MAKE) -C alignment-statistics clean
	$(MAKE) -C convert-bed-positions clean
	$(MAKE) -C count-supporting-reads clean
	$(MAKE) -C index-msa clean
	$(MAKE) -C libpanvc3 clean
	$(MAKE) -C process-alignments clean
	$(MAKE) -C project-alignments clean
	$(MAKE) -C recalculate-mapq clean
	$(MAKE) -C rewrite-cigar clean
	$(MAKE) -C split-alignments-by-reference clean
	$(MAKE) -C subset-alignments clean
	$(MAKE) -C tests clean

clean-dependencies:
	$(MAKE) -C lib/libbio clean
	$(MAKE) -f make/librapidcheck.mk clean
	$(MAKE) -f make/libcatch2.mk clean

clean-dist:
	$(RM) -rf $(DIST_TARGET_DIR)

dependencies: $(DEPENDENCIES)

dist: $(DIST_TAR_GZ)

tests: lib/libbio/lib/rapidcheck/build/librapidcheck.a lib/libbio/lib/Catch2/build/src/libCatch2.a $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C tests

libbio-tests: lib/libbio/tests/tests

libbio: lib/libbio/src/libbio.a

libpanvc3: libpanvc3/libpanvc3.a

alignment-statistics/alignment_statistics: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C alignment-statistics

convert-bed-positions/convert_bed_positions: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C convert-bed-positions
	
count-supporting-reads/count_supporting_reads: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C count-supporting-reads

index-msa/index_msa: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C index-msa

libpanvc3/libpanvc3.a:
	$(MAKE) -C libpanvc3

process-alignments/process_alignments: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C process-alignments

project-alignments/project_alignments: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C project-alignments

recalculate-mapq/recalculate_mapq: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C recalculate-mapq

rewrite-cigar/rewrite_cigar: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C rewrite-cigar

split-alignments-by-reference/split_alignments_by_reference: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C split-alignments-by-reference

subset-alignments/subset_alignments: $(LIBBIO_LIB) libpanvc3/libpanvc3.a
	$(MAKE) -C subset-alignments

$(DIST_TAR_GZ):	$(BUILD_PRODUCTS)
	$(MKDIR) -p $(DIST_TARGET_DIR)
	for f in $(BUILD_PRODUCTS); do if [ $${f##*/} = libpanvc3.a ]; then $(CP) $$f $(DIST_TARGET_DIR); else $(CP) $$f $(DIST_TARGET_DIR)/panvc3_$${f##*/}; fi; done
	$(CP) count-supporting-reads/calculate_reference_bias.py $(DIST_TARGET_DIR)/panvc3_calculate_reference_bias.py
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/cereal/LICENSE $(DIST_TARGET_DIR)/cereal-license.txt
	$(CP) lib/sdsl-lite/LICENSE $(DIST_TARGET_DIR)/sdsl-lite-license.txt
	$(CP) lib/seqan3/LICENSE.md $(DIST_TARGET_DIR)/seqan3-license.md
	$(CP) lib/libbio/lib/GSL/LICENSE $(DIST_TARGET_DIR)/GSL-license.txt
	$(CP) lib/libbio/lib/range-v3/LICENSE.txt $(DIST_TARGET_DIR)/range-v3-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

lib/libbio/local.mk:
	$(CP) make/libbio-local.mk $@

$(LIBBIO_LIB): lib/libbio/local.mk
	$(MAKE) -C lib/libbio

lib/libbio/tests/tests: lib/libbio/local.mk
	$(MAKE) -C lib/libbio tests

lib/libbio/lib/rapidcheck/build/librapidcheck.a:
	$(MAKE) -f make/librapidcheck.mk

lib/libbio/lib/Catch2/build/src/libCatch2.a:
	$(MAKE) -f make/libcatch2.mk
