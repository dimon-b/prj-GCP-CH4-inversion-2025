# ============================================================
# === MIROC4-ACTM CH4 / GCP2025 Bayesian inversion (c_gcpv3/inversion)
#     Intel Fortran + NetCDF Fortran + Intel MKL
#     read_globalview.o is ALWAYS built from read_globalview.f90.
# ============================================================

FC = ifort

NETCDF_ROOT ?= /work/opt/local/apps/intel/2020.4.304/netcdf/4.7.4
FFLAGS ?= -O3 -I$(NETCDF_ROOT)/include
LIBS ?= -L$(NETCDF_ROOT)/lib -lnetcdff -mkl=sequential

# tdi.f is free-form Fortran even though its extension is .f.
FREE_FORM_FLAG = -free

OBJS = bayesinv.o bldfrech.o bldgreens.o datcor.o \
       dumpgreen.o findname.o getconsts.o getdat.o getlgreen.o \
       getnum.o getsrc.o getsub.o output_v3.o srccor.o tdi.o \
       lnblnk.o fnblnk.o defgrps.o getgrps.o dump.o recipe.o getblklst.o \
       split.o zero.o getrecipe.o getmap.o mapinterp.o influence.o \
       diffdat.o intdat.o svdcalc.o inrange.o getmats.o tunevar.o packmat.o \
       read_pseudodata.o read_sourcep.o output_tc3.o \
       read_monthlyobs.o read_sourcetdi.o read_globalview.o tc3greens.o \
       output_tc3tdi.o

.PHONY: all tdi clean check_sources import_legacy
all: tdi_ch4
# Compatible with the previous 'make tdi' invocation.
tdi: tdi_ch4

# Use explicit pattern rules instead of the old .SUFFIXES restriction.
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# The driver is free-form; other .f sources remain fixed-form.
tdi.o: tdi.f tdi.h common.h pointer.h
	$(FC) $(FFLAGS) $(FREE_FORM_FLAG) -c tdi.f -o tdi.o

# Explicit target overrides the old read_globalview.f when it exists.
read_globalview.o: read_globalview.f90
	$(FC) $(FFLAGS) -c read_globalview.f90 -o read_globalview.o

# The original routines include local header files.
$(filter-out tdi.o,$(OBJS)): tdi.h common.h pointer.h

tdi_ch4: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $@

# Check the key build tools and the new Fortran sources before compiling.
check_sources:
	@command -v $(FC) >/dev/null || { echo 'ERROR: compiler $(FC) not on PATH'; exit 1; }
	@test -f $(NETCDF_ROOT)/include/netcdf.inc || { echo 'ERROR: missing netcdf.inc under NETCDF_ROOT'; exit 1; }
	@for header in tdi.h common.h pointer.h; do \
	    test -f "$$header" || { echo "ERROR: missing $$header (copy from codeint54)"; exit 1; }; \
	  done
	@for obj in $(OBJS); do \
	    base=$${obj%.o}; \
	    if ! test -f "$$base.f90" && ! test -f "$$base.f"; then \
	      echo "ERROR: missing Fortran source $$base.f90 or $$base.f"; exit 1; \
	    fi; \
	  done
	@echo 'All 43 Fortran sources, three headers, compiler and NetCDF: OK'

import_legacy:
	bash ./copy_legacy_sources.sh ../../codeint54

clean:
	rm -f *.o *.mod tdi_ch4
