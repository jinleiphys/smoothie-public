# define positions of submodules required

# define the compilers and standard libraries for the current computer

# macbook pro

LIBSTD1 = -L /usr/lib/x86_64-linux-gnu/ -llapack 
#LIBSTD1 = -L /home/moro/cloud/THO/thox/v2.5/lapack64/ -llapack 
#LIBSTD2 = -lhdf5_hl -lhdf5hl_fortran -lhdf5 -lhdf5_fortran -lz
LIBSTD3 =
LIBSTD4 =
LIBSTD5 =
LIBSTD6 =
LIBSTD7 =
LIBSTD8 =

LIBSTD = $(LIBSTD1) $(LIBSTD2) $(LIBSTD3) $(LIBSTD4) $(LIBSTD5) $(LIBSTD6) $(LIBSTD7) $(LIBSTD8)

CC = gfortran
COPT = -O3
FC = gfortran
F90 = gfortran
COMPILE_OPT1= -fPIC -fopenmp  -O3   
COMPILE_OPT2= -ffixed-line-length-0
#COMPILE_OPT3= -I /opt/local/include/
#COMPILE_OPT4= -Wunused-variable
#COMPILE_OPT5=-ffpe-trap=invalid,zero,overflow
COMPILE_OPT6=
COMPILE_OPT7=
# COMPILE_OPT8=-fdefault-real-8  -fno-automatic

FOPT = $(COMPILE_OPT1) $(COMPILE_OPT2) $(COMPILE_OPT3) $(COMPILE_OPT4) \
       $(COMPILE_OPT5) $(COMPILE_OPT6) $(COMPILE_OPT7) $(COMPILE_OPT8) \
      -DVERDATE="'$(VERDATE)'" -DVERREV="'$(VERREV)'" \
      -DMAKEF90="'$(F90)'" -DMAKEF77="'$(FC)'" \
      -DMAKEFOPT1="'$(COMPILE_OPT1)'" \
      -DMAKEFOPT2="'$(COMPILE_OPT2)'" \
      -DMAKEFOPT3="'$(COMPILE_OPT3)'" \
      -DMAKEFOPT4="'$(COMPILE_OPT4)'" \
      -DMAKEFOPT5="'$(COMPILE_OPT5)'" \
      -DMAKEFOPT6="'$(COMPILE_OPT6)'" \
      -DMAKEFOPT7="'$(COMPILE_OPT7)'" \
      -DMAKEFOPT8="'$(COMPILE_OPT8)'" \
      -DMAKECOPT="'$(COPT)'" -DMAKECC="'$(CC)'" \
      -DMAKELIBSTD1="'$(LIBSTD1)'" \
      -DMAKELIBSTD2="'$(LIBSTD2)'" \
      -DMAKELIBSTD3="'$(LIBSTD3)'" \
      -DMAKELIBSTD4="'$(LIBSTD4)'" \
      -DMAKELIBSTD5="'$(LIBSTD5)'" \
      -DMAKELIBSTD6="'$(LIBSTD6)'" \
      -DMAKELIBSTD7="'$(LIBSTD7)'" \
      -DMAKELIBSTD8="'$(LIBSTD8)'" \
      -DCOMPDATE="'$(COMPDATE)'"
AR     = ar
ARFLAGS= rv
RANLIB   = touch
#  -ffpe-trap=invalid,zero,overflow  -fcheck=all -fno-automatic -finit-integer=-1000 -finit-real=inf -fbounds-check\
# -ffpe-trap=invalid,zero,overflow  -fcheck=all -fno-automatic -finit-integer=-1000 -finit-real=zero\
#-fdefault-real-8 -fno-automatic  -ffpe-trap=invalid,zero,overflow  -fbounds-check -fcheck=all -fcheck=bounds
#   -ffpe-trap=invalid,zero,overflow \
#      -finit-integer=-1000 -finit-real=inf\
#      -fcheck=all \
