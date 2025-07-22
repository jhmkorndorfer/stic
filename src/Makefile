CPU     = $(shell uname -m)
OS      = $(shell uname -s)

include makefiles/makefile.$(CPU).$(OS)

RHFOLD = rh/

# Debug flags
FFLAGS += -g3
CXXFLAGS += -g3

# Paths and includes
INCLUDE = $(shell ncxx4-config --cflags) -I$(GSLPATH)/include -I./ -I$(RHFOLD) -I$(RHFOLD)/rh_1d
LIBS_FDENS = $(shell ncxx4-config --libs) -lstdc++ -L$(GSLPATH)/lib/ -ltirpc
LIBS = $(LIBS_FDENS) -L./ -L$(RHFOLD)/rh_1d/ -lrhf1d -L$(GSLPATH) -lfftw3 -ltirpc

# Machine-specific binary
MNAME = $(shell uname -n)
STMAC = STiC_$(MNAME).x

# Sources
FFILES = eos_math_special.o eos_eqns.o eos.o 
OFILES_SPARSE = input.o clm.o cop.o witt.o ceos.o comm.o depthmodel.o spectral.o fpigen.o specrebin.o specprefilter.o fpi.o atmosphere.o clte.o crh.o io.o slave.o master_sparse.o main_sparse.o

FDENS = cop.o ceos.o io.o depthmodel.o fillDensities.o

# Patterns
.SUFFIXES: .o .f90 .cc

%.o: %.cc
	$(MPCXX) $(CXXFLAGS) $(OPENMP) $(INCLUDE) -MMD -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

SUBDIRS = rh rh/rh_1d

.PHONY: default all $(SUBDIRS)

# Default target executed when you just run 'make'
default: $(SUBDIRS) STiC

# Rule to build each subdirectory
$(SUBDIRS):
	$(MAKE) -C $@
	
# Link targets
STiC: $(FFILES) $(OFILES_SPARSE)
	$(LINKER) -o STiC.x $(OPTS) $(OPENMP) $^ $(INCLUDE) $(LIBS) $(LINKEROPTS)

machine: $(FFILES) $(OFILES_SPARSE)
	$(LINKER) -o $(STMAC) $(OPTS) $(OPENMP) $^ $(INCLUDE) $(LIBS) $(LINKEROPTS)

fillDensities: $(FFILES) $(FDENS)
	$(LINKER) -o fillDensities.x $(CXXFLAGS) $(OPENMP) $(FFILES) $(FDENS) $(LIBS_FDENS) $(INCLUDE) $(LINKEROPTS)

# Clean rule
clean:
	rm -f *.o *.mod *.d STiC.x $(STMAC) fillDensities.x

# Auto-include generated dependencies
-include $(wildcard *.d)
