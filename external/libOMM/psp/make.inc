ROOTDIR = ..
include $(dir $(lastword $(MAKEFILE_LIST)))/$(ROOTDIR)/make.inc

LIBNAME = $(PSPLIBNAME)
BUILDPATH = $(PSPBUILDPATH)
