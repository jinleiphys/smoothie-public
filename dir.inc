# this file contains definitions of names of the subdirectories
GENMOD = $(BASE)/general_modules
POTMOD = $(BASE)/pot_modules
PWMOD = $(BASE)/pw_modules
MCMOD = $(BASE)/mc_modules
MESHMOD = $(BASE)/mesh_modules
HDFMOD = $(BASE)/hdf_modules
HOMOD = $(BASE)/ho_modules
NNMOD = $(BASE)/NNpotentiale
SMOOTHIEMOD= $(BASE)/smoothie
YOGURTMOD=$(BASE)/yogurt
SMOOTHIEHMMOD= $(BASE)/smoothieHM
CHOCOLATEMOD= $(BASE)/chocolate
OPTPOT= $(BASE)/NApotentiale
SMOOTHIEXMOD= $(BASE)/smoothieX
COOKIEMOD= $(BASE)/cookie


.PHONY : general_modules
general_modules:
	make -C$(GENMOD)

.PHONY : pot_modules
pot_modules:
	make -C$(POTMOD)

.PHONY : pw_modules
pw_modules:
	make -C$(PWMOD)

.PHONY : mesh_modules
mesh_modules:
	make -C$(MESHMOD)





.PHONY : src
src:
		make  -C$(SRCMOD)


.SUFFIXES: .F90 .f90

.F90.o          :
	$(F90) $(FFLAGS) -c $<

.f90.o          :
	$(F90) $(FFLAGS) -c $<
