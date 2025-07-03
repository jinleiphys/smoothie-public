# make file generates objects and modules in pw_modules
objects=clebsch.o coul90.o spharm.o whittaker.o coulcc.o

# Reset the default goal.
.DEFAULT_GOAL := libpwmodules.a

BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
# VERDATE := $(shell expr "$(VERDATE)" : '\(.*\)(.*).*')
VERREV := $(shell git log -1 --pretty=format:"%h")


include $(BASE)/make.inc
include $(BASE)/dir.inc

LIB=$(LIBSTD) -L $(GENMOD) -lgeneralmodules \
              -L $(PWMOD)  -lpwmodules
CFLAGS=$(COPT) -I$(GENMOD) -I$(PWMOD)
FFLAGS=$(FOPT) -I$(GENMOD) -I$(PWMOD)

libpwmodules.a: general_modules $(objects)
	$(AR) $(ARFLAGS) libpwmodules.a $(objects)
	$(RANLIB) libpwmodules.a


# clean target
.PHONY : clean
clean: cleanhere
	make clean -C$(GENMOD)

.PHONY : cleanhere
cleanhere:
	rm -f $(objects) *.mod *.a *.o
