MAKE_SCRIPT_DIR_ := $(dir $(lastword $(MAKEFILE_LIST)))
OS_NAME := $(shell $(MAKE_SCRIPT_DIR_)/../tools/os_name.sh)
