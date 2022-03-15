# Makefile for the fireHK program

FC=gfortran
FCFLAGS=-ffree-form -ffree-line-length-none

netcdf=/home/public/easybuild/software/netCDF-Fortran/4.5.3-gmpich-2021.01
NC_LIB=$(netcdf)/lib
NC_INC=$(netcdf)/include

CPPFLAGS = -I$(NC_INC)
LDFLAGS  = -L$(NC_LIB)
LIBS     = -lnetcdff

#---------------------------------------------

OBJS = parametersmod.f90\
       statevarsmod.f90\
       errormod.f90\
			 geohashmod.f90\
			 statsmod.f90\
			 utilitiesmod.f90\
			 randomdistmod.f90\
			 getdatamod.f90\
			 newsplinemod.f90\
			 netcdfinputmod.f90\
			 weathergenmod.f90\
			 gwgenmod.f90\
			 orbitmod.f90\
			 daylengthmod.f90\
			 radiationmod.f90\
			 diurnaltempmod.f90\
			 ffdimod.f90\
			 gfdimod.f90\
			 main.f90

#---------------------------------------------

.SUFFIXES: .o .f90 .f .mod

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $<

all::	fireHK

fireHK: $(OBJS1)
	$(FC) $(FCFLAGS) -o fireHK $(OBJS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

clean::
	-rm fireHK *.o *.mod
