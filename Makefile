###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR     - Root directory of the MOOSE project
# FRAMEWORK_DIR - Location of the MOOSE framework
#
###############################################################################
#//THESE VARIABLES ARE SET ONLY IF THEY ARE NOT DEFINED (OPERATOR ?=)
EXAMPLE_DIR        ?= $(shell dirname `pwd`)
MOOSE_DIR          ?= /Users/liuy2/Documents/projects/moose
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
MODULE_DIR	   ?= $(MOOSE_DIR)/modules
###############################################################################


# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk


# dependency on MODULES
ALL_MODULES := yes
include $(MODULE_DIR)/modules.mk




APPLICATION_NAME := lubrication
# dep apps
APPLICATION_DIR    := $(shell pwd)
APPLICATION_NAME   := ASFracture
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk


# Include dependency files for this example
ex_srcfiles := $(shell find $(APPLICATION_DIR) -name "*.C")
ex_deps     := $(patsubst %.C, %.$(obj-suffix).d, $(ex_srcfiles))
-include $(ex_deps)

###############################################################################
# Additional special case targets should be added here

