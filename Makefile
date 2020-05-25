#############################################################################
##
## $Id: Makefile,v 1.38 2009/08/07 14:47:08 sahughes Exp sahughes $
##
#############################################################################
#
# Makefile for GremlinEq
#
#############################################################################

.PHONY : dummy

#############################################################################
#
# All the paths we need to know.
#
#############################################################################

TOP = $(CURDIR)
INCL = $(TOP)/include
INCLGMP = /usr/local/include
INCLHDF5 = /usr/local/include
BIN = $(TOP)/bin
LIB = $(TOP)/lib
LIBGMP = /usr/local/lib
LIBHDF5 = /usr/local/lib

SRC = $(TOP)/src
CIRCEQSRC = $(SRC)/circeq
EXECSRC = $(SRC)/exec
SWSHSRC = $(SRC)/swsh
FTSRC = $(SRC)/fujtag
UTILSRC = $(SRC)/utility
ALLSRCS = $(CIRCEQSRC):$(EXECSRC):$(SWSHSRC):$(FTSRC):$(UTILSRC)

VPATH = $(BIN):$(INCL):$(LIB):$(ALLSRCS)

#############################################################################
#
# Useful variables.
#
#############################################################################

CC = g++

AR = ar rv

SYSLIBS = -lhdf5 -lhdf5_hl -lgsl -lgslcblas -lfftw3 -lm -lgmp

CFLAGS = -O2 -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated

#############################################################################
#
# -Wall to catch as many warnings as possible, except:
#
# -Wno-unused because some Numerical Recipes functions don't use all of
#             their parameters, and
#
# -Wno-uninitialized because a few Numerical Recipes functions have
#                    variables that look like they'll be used uninitialized.
#                    (They won't be, don't worry about it.)
#
# -Wno-deprecated because on newer compilers, my C++ generates complaints
#                 about the header format.
#
#############################################################################

#############################################################################
#
# All .o files.
#
#############################################################################

CIRCEQOBJS = CEDR.o CEID.o CEKG.o CEKR.o CETD.o

GKGOBJS = GKG.o

RRGWOBJS = RRGW.o

TIDALHOBJS = TidalH.o

SWSHOBJS = SWSHCGUtil.o SWSHSpherical.o SWSHSpheroid.o

FTOBJS = FT.o fsum.o gammln.o hypergeom.o radialfrac.o renangmom.o specialradial.o

.INTERMEDIATE : $(CIRCEQOBJS) $(GKGOBJS) $(RRGWOBJS) $(TIDALHOBJS) $(SWSHOBJS) $(FTOBJS)

#############################################################################
#
# All .cc files used to generate libraries.
#
#############################################################################

CIRCEQCC = CEDR.cc CEID.cc CEKG.cc CEKR.cc CETD.cc

GKGCC = GKG.cc

RRGWCC = RRGW.cc

TIDALHCC = TidalH.cc

SWSHCC = SWSHCGUtil.cc SWSHSpherical.cc SWSHSpheroid.cc

FTCC = FT.cc fsum.cc gammln.cc hypergeom.cc radialfrac.cc renangmom.cc \
	specialradial.cc

#############################################################################
#
# Begin setting up dependencies.
#
#############################################################################

all : core #horizgeom

core : Circ_Eq Circ_Eq_Seq Circ_Eq_Seq2 Circ_Eq_Traj Circ_Eq_TotFlux_r Circ_Eq_TotFlux_v Circ_Eq_Ymode_v Circ_Eq_Smode_v Circ_Eq_lmode_v Circ_Eq_mmode_v Circ_Eq_Wave Circ_Eq_Clm Circ_Eq_Clm2 Fluxes GW Wrapper

horizgeom : EmbedEq EmbedEq_lmode EmbedSurf TidalHEq TidalHSurf

#############################################################################
#
# Dependencies for executables.
#
#############################################################################

Circ_Eq : Circ_Eq.cc -lCircEq -lSWSH -lGKG -lFT -lRRGW Globals.h CEKG.h CEKR.h Tensors.h RRGW.h FT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq.cc -o $(BIN)/Circ_Eq -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lGKG -lFT -lRRGW $(SYSLIBS)

Circ_Eq_Seq : Circ_Eq_Seq.cc -lCircEq -lSWSH -lGKG -lFT -lRRGW Globals.h CEKG.h CEKR.h CETD.h Tensors.h RRGW.h FT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Seq.cc -o $(BIN)/Circ_Eq_Seq -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lGKG -lFT -lRRGW $(SYSLIBS)

 Circ_Eq_Seq2 : Circ_Eq_Seq2.cc -lCircEq -lSWSH -lGKG -lFT -lRRGW Globals.h CEKG.h CEKR.h CETD.h Tensors.h RRGW.h FT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Seq2.cc -o $(BIN)/Circ_Eq_Seq2 -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lGKG -lFT -lRRGW $(SYSLIBS)

Circ_Eq_Traj : Circ_Eq_Traj.cc -lCircEq -lSWSH -lRRGW Globals.h CEID.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Traj.cc -o $(BIN)/Circ_Eq_Traj -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_TotFlux_r : Circ_Eq_TotFlux_r.cc -lCircEq -lSWSH -lRRGW Globals.h CEID.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_TotFlux_r.cc -o $(BIN)/Circ_Eq_TotFlux_r -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_TotFlux_v : Circ_Eq_TotFlux_v.cc -lCircEq -lSWSH -lRRGW Globals.h CEID.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_TotFlux_v.cc -o $(BIN)/Circ_Eq_TotFlux_v -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_Ymode_v : Circ_Eq_Ymode_v.cc -lCircEq -lSWSH -lRRGW Globals.h CEID.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Ymode_v.cc -o $(BIN)/Circ_Eq_Ymode_v -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_Smode_v : Circ_Eq_Smode_v.cc -lCircEq -lSWSH -lRRGW Globals.h CEID.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Smode_v.cc -o $(BIN)/Circ_Eq_Smode_v -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_lmode_v : Circ_Eq_lmode_v.cc -lCircEq -lSWSH -lRRGW Globals.h CEID.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_lmode_v.cc -o $(BIN)/Circ_Eq_lmode_v -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_mmode_v : Circ_Eq_mmode_v.cc -lCircEq -lSWSH -lRRGW Globals.h CEID.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_mmode_v.cc -o $(BIN)/Circ_Eq_mmode_v -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_Wave : Circ_Eq_Wave.cc -lCircEq -lSWSH -lRRGW Globals.h Tensors.h SWSH.h RRGW.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Wave.cc -o $(BIN)/Circ_Eq_Wave -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_Clm : Circ_Eq_Clm.cc -lCircEq -lSWSH -lRRGW Globals.h Tensors.h SWSH.h RRGW.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Clm.cc -o $(BIN)/Circ_Eq_Clm -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Circ_Eq_Clm2 : Circ_Eq_Clm2.cc -lCircEq -lSWSH -lRRGW Globals.h Tensors.h SWSH.h RRGW.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Clm2.cc -o $(BIN)/Circ_Eq_Clm2 -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

Fluxes : Fluxes.cc -lCircEq -lSWSH -lRRGW Globals.h CEDR.h Tensors.h SWSH.h RRGW.h
	$(CC) $(CFLAGS) $(EXECSRC)/Fluxes.cc -o $(BIN)/Fluxes -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

GW : GW.cc -lCircEq -lSWSH -lRRGW Globals.h CEDR.h RRGW.h SWSH.h Tensors.h
	$(CC) $(CFLAGS) $(EXECSRC)/GW.cc -o $(BIN)/GW -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lRRGW $(SYSLIBS)

TidalHSurf : TidalHSurf.cc -lCircEq -lFT -lGKG -lSWSH  Globals.h CEKR.h CEKG.h FT.h SWSH.h Tensors.h fsum.h
	$(CC) $(CFLAGS) $(EXECSRC)/TidalHSurf.cc -o $(BIN)/TidalHSurf -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lFT -lGKG -lSWSH  $(SYSLIBS)

TidalHEq : TidalHEq.cc -lCircEq -lFT -lGKG -lTidalH -lSWSH  Globals.h CEKR.h CEKG.h FT.h SWSH.h Tensors.h fsum.h
	$(CC) $(CFLAGS) $(EXECSRC)/TidalHEq.cc -o $(BIN)/TidalHEq -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lFT -lGKG -lTidalH -lSWSH  $(SYSLIBS)

EmbedEq : EmbedEq.cc -lCircEq -lFT -lGKG -lTidalH -lSWSH  Globals.h TidalH.h CEKR.h CEKG.h FT.h SWSH.h Tensors.h fsum.h
	$(CC) $(CFLAGS) $(EXECSRC)/EmbedEq.cc -o $(BIN)/EmbedEq -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lFT -lGKG -lTidalH -lSWSH  $(SYSLIBS)

EmbedEq_lmode : EmbedEq_lmode.cc -lCircEq -lFT -lGKG -lTidalH -lSWSH  Globals.h TidalH.h CEKR.h CEKG.h FT.h SWSH.h Tensors.h fsum.h
	$(CC) $(CFLAGS) $(EXECSRC)/EmbedEq_lmode.cc -o $(BIN)/EmbedEq_lmode -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lFT -lGKG -lTidalH -lSWSH  $(SYSLIBS)

EmbedSurf : EmbedSurf.cc -lCircEq -lFT -lGKG -lTidalH -lSWSH  Globals.h TidalH.h CEKR.h CEKG.h FT.h SWSH.h Tensors.h fsum.h
	$(CC) $(CFLAGS) $(EXECSRC)/EmbedSurf.cc -o $(BIN)/EmbedSurf -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lFT -lGKG -lTidalH -lSWSH  $(SYSLIBS)

Wrapper : Wrapper.cc -lCircEq -lSWSH -lGKG -lFT -lRRGW Globals.h CEKG.h CEKR.h Tensors.h RRGW.h FT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Wrapper.cc -o $(BIN)/Wrapper -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -L$(LIB) -L$(LIBGMP) -L$(LIBHDF5) -lCircEq -lSWSH -lGKG -lFT -lRRGW $(SYSLIBS)

#############################################################################
#
# Dependencies for libraries.
#
#############################################################################

-lCircEq : $(CIRCEQCC) $(CIRCEQOBJS) Globals.h CEDR.h CEID.h CEKG.h CEKR.h CETD.h Tensors.h SWSH.h FT.h
	$(AR) $(LIB)/libCircEq.a $(CIRCEQOBJS)

-lGKG : $(GKGCC) $(GKGOBJS) Globals.h GKG.h Tensors.h
	$(AR) $(LIB)/libGKG.a $(GKGOBJS)

-lRRGW : $(RRGWCC) $(RRGWOBJS) Globals.h RRGW.h SWSH.h
	$(AR) $(LIB)/libRRGW.a $(RRGWOBJS)

-lTidalH : $(TIDALHCC) $(TIDALHOBJS) Globals.h TidalH.h SWSH.h Tensors.h
	$(AR) $(LIB)/libTidalH.a $(TIDALHOBJS)

-lSWSH : $(SWSHCC) $(SWSHOBJS) Globals.h SWSH.h
	$(AR) $(LIB)/libSWSH.a $(SWSHOBJS)

-lFT: $(FTCC) $(FTOBJS) Globals.h FT.h
	$(AR) $(LIB)/libFT.a $(FTOBJS)

#############################################################################
#
# Rule for implicitly generating intermediate .o from .cc (needed to make
# libaries).
#
#############################################################################

%.o : %.cc
	$(CC) $(CFLAGS) -I$(INCL) -I$(INCLGMP) -I$(INCLHDF5) -c $< -o $@

#############################################################################
#
# make clean
#
#############################################################################

clean : dummy
	$(RM) $(BIN)/*
	$(RM) $(LIB)/*
