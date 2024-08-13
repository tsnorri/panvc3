MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include  $(MAKE_SCRIPT_DIR)/../../make/os-name.mk
-include $(MAKE_SCRIPT_DIR)/../../local.mk
include  $(MAKE_SCRIPT_DIR)/../../make/$(OS_NAME).mk
include  $(MAKE_SCRIPT_DIR)/../../make/gcc.mk
