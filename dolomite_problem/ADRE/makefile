# fortran compiler
# FC90 = gfortran
FC90 = gfortran-4.9
	# for some reason gfortran 8.1.0 is causing problems with file 'open' statements
	# thus, I am using gfortran 4.9.4

# PHREEQCRM library flag
	# (if you do the default install, this is where it will be on macOS)
LIB = -L/usr/local/lib -lphreeqcrm

# optimization level
FLAGS = -O3

# for debugging
# FLAGS = -g -fcheck=all -fbacktrace -Og
# generate warnings and direct to compile.log
# WARN = -Wall > compile.log 2>&1

# add LIB to FLAGS
FLAGS += $(LIB)
# add WARN to FLAGS
FLAGS += $(WARN)

adre_files = dolo_ADRE.f90 ADRE_mod.f90 RM_interface.F90

all: adre_exe

adre_exe: $(adre_files)
	$(FC90) RM_interface.F90 -c $(FLAGS)
	$(FC90) ADRE_mod.f90 -c $(FLAGS)
	$(FC90) $(adre_files) -o $@ $(FLAGS)

clean:
	rm *_exe *.mod *.o *.out *.err dolomite_chem.in myrun* error.inp
