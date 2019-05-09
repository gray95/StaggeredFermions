# $Id: Makefile,v 1.1 2007/02/25 15:21:41 cmcneile Exp cmcneile $
#
# Makefile to link in chroma/qdp++ library
# functions to a main program.
#
#  This makefile uses commands specific
#  to the gnu version of make.
#
#  This allows you to build a main program
#  against an installed grid library.
#
#  1) The location of the  program has to be set.
#  2) The main program you want to compile needs to be 
#     added to this makefile.
#
#

# The grid utility to extract the compile flags
gridinstall=/home/gray/code/grid
config=$(gridinstall)/bin/grid-config
#
#  flags from the grid utility program
#
CXXFLASGS= $(shell $(config) --cxxflags)   -I$(gridinstall)/include
LDFLAGS=$(shell $(config) --ldflags)
LIBS= $(shell $(config) --libs) -L$(gridinstall)/lib  -lGrid
#COMPILER=$(shell $(config) --cxx)
COMPILER= g++
# other libraries from SciDac


pion_corr: staggered_pion.cc 
	$(COMPILER) -g -o $@ $(CXXFLASGS)   staggered_pion.cc   $(LDFLAGS) $(LIBS) 

rho_corr: rho_corr.cc 
	$(COMPILER) -g -o $@ $(CXXFLASGS)  rho_corr.cc   $(LDFLAGS) $(LIBS) 

lab: laboratory.cc
	$(COMPILER) -g -o $@ $(CXXFLASGS)   laboratory.cc   $(LDFLAGS) $(LIBS) 

clean: 
	rm  pion_corr rho_corr lab