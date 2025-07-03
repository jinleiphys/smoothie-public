# make file generates objects and modules in pw_modules
objects= precision.o constants.o systems.o channels.o 



BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
# VERDATE := $(shell expr "$(VERDATE)" : '\(.*\)(.*).*')
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)

# Reset the default goal.
.DEFAULT_GOAL := libgeneralmodules.a

include $(BASE)/make.inc
include $(BASE)/dir.inc

LIB=$(LIBSTD) -L $(GENMOD) -lgeneralmodules
CFLAGS=$(COPT) -I$(GENMOD)
FFLAGS=$(FOPT) -I$(GENMOD)


libgeneralmodules.a: $(objects)
	$(AR) $(ARFLAGS) libgeneralmodules.a $(objects)
	$(RANLIB) libgeneralmodules.a

test_general: libgeneralmodules.a test_general.o
	$(F90)  $(FFLAGS)  -o test_general  test_general.o $(LIB)

# clean target
.PHONY : clean
clean:
	rm -f $(objects) *.mod *.a test_general.o test_general
