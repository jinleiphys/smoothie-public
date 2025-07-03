objects = input.o lagrange_mesh.o scatt.o green.o  bound.o fuspot.o zerorange.o iavzerospin.o  iavdwba.o \
iavdwbarbx.o \
lagrange_mesh_source.o iavdwbarbx_lagrange.o iavdwbarb_lagrange.o

BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
# VERDATE := $(shell expr "$(VERDATE)" : '\(.*\)(.*).*')
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)

# Reset the default goal.
.DEFAULT_GOAL := all


#define compiler options etc.
include $(BASE)/make.inc

#define subdirectories
include $(BASE)/dir.inc

LIB=$(LIBSTD) -L $(SMOOTHIEMOD) -lsmoothie \
              -L $(GENMOD) -lgeneralmodules \
              -L $(PWMOD) -lpwmodules \
              -L $(MESHMOD) -lmeshmodules \
              -L $(POTMOD) -lpotmodules \
              -L $(SMOOTHIEMOD) -lsmoothie \
              -lblas -llapack

CFLAGS=$(COPT) -I$(GENMOD) -I$(PWMOD) -I$(MESHMOD) -I$(POTMOD) -I$(SMOOTHIEMOD)
FFLAGS=$(FOPT) -I$(GENMOD) -I$(PWMOD) -I$(MESHMOD) -I$(POTMOD) -I$(SMOOTHIEMOD)



.SUFFIXES: .F90 .f90

all: libsmoothie.a smoothie


libsmoothie.a:  general_modules pw_modules mesh_modules pot_modules $(objects)
	$(AR) $(ARFLAGS) libsmoothie.a $(objects)
	$(RANLIB) libsmoothie.a

smoothie:  libsmoothie.a smoothie.o
	$(FC) -o smoothie $(FFLAGS)  smoothie.o $(LIB)


.PHONY : clean
clean: cleanhere
	make clean -C$(GENMOD)
	make clean -C$(PWMOD)
	make clean -C$(MESHMOD)
	make clean -C$(POTMOD)


.PHONY : cleanhere
cleanhere:
	rm -f  *.o *.mod *.a smoothie


.F90.o          :
	$(F90) $(FFLAGS) -c $<

.f90.o          :
	$(F90) $(FFLAGS) -c $<










#
