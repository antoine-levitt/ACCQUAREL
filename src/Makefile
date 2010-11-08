### MAKEFILE PREAMBLE ###
F95 = g95
#F95 = gfortran
#F95 = ifort

# common to compilers and linker
DEBUGFFLAGS = -fbounds-check -fimplicit-none -ffpe-trap=invalid,zero,denormal
PRODFFLAGS = -fimplicit-none
COMMONFLAGS = -O3 -g
FFLAGS = $(PRODFFLAGS) $(COMMONFLAGS)

ifndef ASPICROOT
   ASPICROOT=../A.S.P.I.C
endif

ifeq ($(F95),g95)
   F95FLAGS = -c $(FLAGS) -fmod=$(MODDIR) $(FFLAGS)
   LD = g95
   FLIBS = -llapack -lblas -lm
endif
ifeq ($(F95),gfortran)
   F95FLAGS = -c $(FLAGS) -fopenmp -J $(MODDIR) $(FFLAGS)
   LD = gfortran
   FLIBS = -lgomp -llapack -lblas -lm
endif
ifeq ($(F95),ifort)
   F95FLAGS = -c $(FLAGS) -module $(MODDIR) $(FFLAGS)
   LD = ifort
   FLIBS = -llapack -lblas -lm
endif

LDFLAGS = $(COMMONFLAGS)

CPP = g++
CXXFLAGS = -c $(COMMONFLAGS)

## DIRECTORIES ##
TARGET = accquarel.exe
EXEDIR = ../
MODDIR = ../mod/
OBJDIR = ../obj/
SRCDIR = ./

## ASPIC code ##
ASPICINCPATH = -I$(ASPICROOT)/include

ASPICLIBPATH = -L$(ASPICROOT)/lib

XERCESCLIBPATH = -L$(XERCESCROOT)/lib

ifeq ($(shell uname),Darwin)
## cluster osx darwin
   ASPICLIBS = $(ASPICLIBPATH) -lgIntegrals -lchemics -lxmlParser -lgaussian -lpolynome -lcontainor -laspicUtils $(XERCESCLIBPATH) -lxerces-c -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
endif
ifeq ($(shell uname),Linux)
## linux ubuntu
   ASPICLIBS = $(ASPICLIBPATH) -lgIntegrals -lchemics -lxmlParser -lgaussian -lpolynome -lcontainor -laspicUtils $(XERCESCLIBPATH) -lxerces-c -L/usr/lib/gcc/i686-linux-gnu/4.4/ -lstdc++
endif

## ACCQUAREL code ##
OBJ =				\
 $(OBJDIR)expokit.o		\
 $(OBJDIR)optimization.o	\
 $(OBJDIR)rootfinding.o		\
 $(OBJDIR)tools.o		\
 $(OBJDIR)setup.o		\
 $(OBJDIR)common.o      	\
 $(OBJDIR)integrals_c.o		\
 $(OBJDIR)integrals_f.o		\
 $(OBJDIR)basis.o		\
 $(OBJDIR)matrices.o		\
 $(OBJDIR)scf.o			\
 $(OBJDIR)roothaan.o		\
 $(OBJDIR)levelshifting.o	\
 $(OBJDIR)diis.o		\
 $(OBJDIR)oda.o			\
 $(OBJDIR)esa.o			\
 $(OBJDIR)algorithms.o		\
 $(OBJDIR)drivers.o		\
 $(OBJDIR)main.o

# Compilation rules

$(TARGET) :  $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $(EXEDIR)$(TARGET) $(FLIBS) $(ASPICLIBS)
	@echo " ----------- $(TARGET) created ----------- "

$(OBJDIR)%.o : $(SRCDIR)%.f90
	$(F95) $(F95FLAGS) -o $@ $<

$(OBJDIR)%.o : $(SRCDIR)%.f
	$(F95) $(F95FLAGS) -o $@ $<

$(OBJDIR)%.o : $(SRCDIR)%.cpp
	$(CPP) $(CXXFLAGS) -o $@ $(ASPICINCPATH) $<
#
clean :
	rm $(EXEDIR)$(TARGET) $(OBJ) $(MODDIR)*.mod
