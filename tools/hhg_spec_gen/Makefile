FC = gfortran
#FC = ifort
#LIBS = -lfftw3
#LIBS= -O3 -L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64 -mkl -traceback -guide-opts="transform.f90, 1-105"
#LIBS= -O3 -L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64 -mkl -traceback
#LIBS= -O3 -L$(MKLROOT) -mkl -traceback
LIBS = -O3  -L${MKLROOT}/lib/intel64/ -L${MKLROOT} -Wl, -lfftw3 -lm -ldl

#LIBS= -O3 -mkl -traceback
#FLAGS = -O3 -traceback -guide-opts="transform.f90, 1-100"
FLAGS = -O3  -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -c

OBJECTS = init.o transform.o main.o
#globalvars.o globalconst.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all:	$(OBJECTS)
	$(FC) -o FT_transform $(OBJECTS) $(LIBS)

.PHONY : clean
clean :
	rm -f $(OBJECTS) *.mod FT_transform *.txt *.tmp
