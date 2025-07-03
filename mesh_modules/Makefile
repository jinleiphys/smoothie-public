# make file generates objects and modules in pw_modules
objects=gauss.o interpolation.o mesh.o angularmesh.o derivative.o 

# Reset the default goal.
.DEFAULT_GOAL := libmeshmodules.a

BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
# VERDATE := $(shell expr "$(VERDATE)" : '\(.*\)(.*).*')
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)



include $(BASE)/make.inc
include $(BASE)/dir.inc

LIB=$(LIBSTD) -L $(GENMOD) -lgeneralmodules \
              -L $(MESHMOD) -lmeshmodules
CFLAGS=$(COPT) -I$(MESHMOD) -I$(GENMOD)
FFLAGS=-I$(MESHMOD) -I$(GENMOD) $(FOPT)


libmeshmodules.a:  general_modules $(objects)
	$(AR) $(ARFLAGS) libmeshmodules.a $(objects)
	$(RANLIB) libmeshmodules.a



# clean target
.PHONY : clean
clean: cleanhere
	make clean -C$(GENMOD)

.PHONY : cleanhere
cleanhere:
	rm -f $(objects) *.mod *.a
