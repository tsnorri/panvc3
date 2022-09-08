MAKE_SCRIPT_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
OS_NAME := $(shell $(MAKE_SCRIPT_DIR)/../tools/os_name.sh)
