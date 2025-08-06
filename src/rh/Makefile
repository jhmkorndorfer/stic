# ---- Optimized RH Makefile for fast development ----

CPU := $(shell uname -m)
OS := $(shell uname -s)

include ../makefiles/makefile.$(CPU).$(OS)

# Add debug symbols
CFLAGS   += -fopenmp
CXXFLAGS += -fopenmp
FFLAGS   += -fopenmp

CFLAGS += -I../rh  # <--- ADD THIS if atom.h is in ../rh/
CFLAGS += -I../rh/rh_1d  # <--- if some headers are in rh_1d/

# Include paths
RHINC = -I../ -I../eigen3
CPPFLAGS += $(RHINC)

# Object files (deduplicated)
OBS = readAtomFile.o solveLinearCXX.o abundance.o accelerate.o background.o backgropac_xdr.o barklem.o broad.o brs_xdr.o chemequil.o cocollisions.o collision.o complex.o cubeconvol.o duplicate.o error.o expint.o expspline.o fillgamma.o fixedrate.o fpehandler.o gammafunc.o gaussleg.o getcpu.o getlambda.o getline.o giigen.o h2collisions.o hunt.o humlicek.o hydrogen.o initial_xdr.o initscatter.o kurucz.o linear.o ltepops.o ludcmp.o matrix.o maxchange.o metal.o molzeeman.o nemetals.o ohchbf.o opacity.o options.o order.o parse.o paschen.o planck.o pops_xdr.o profile.o radrate_xdr.o readatom.o readb_xdr.o readj.o readvalue.o rayleigh.o readinput.o readmolecule.o scatter.o solvene.o sortlambda.o spline.o statequil.o statequil_H.o stokesopac.o stopreq.o thomson.o vacuumtoair.o voigt.o w3.o wigner.o writeatmos_xdr.o writeatom_xdr.o writecoll_xdr.o writedamp_xdr.o writeinput_xdr.o writemetal_xdr.o writemolec_xdr.o writeopac_xdr.o writespect_xdr.o zeeman.o hui_.o humlicek_.o

.SUFFIXES: .o .f90 .c .cc

# Pattern rules with automatic dependency tracking
%.o: %.cc
	$(CXX) $(CXXFLAGS) -MMD -MP $(CPPFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -MMD -MP $(CPPFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -MMD -MP -c $< -o $@

# Default target: compile all object files
all: objects

objects: $(OBS)

# Clean rule
clean:
	rm -f *.o *.d *.mod

# Automatically include dependency files
-include $(wildcard *.d)
