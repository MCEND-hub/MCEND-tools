#F2C = intelem
C2C = intelem
#F2C = gnu95
F2C = gfortran
C2C = gcc
#FLAGS = -lblas -llapack
FC = ifort
FC = gfortran
F2FLAGS = -O3 #-traceback --no-latex-doc
F2LINK_FLAGS = -O3

FTC = $(FC)
#ifort
#FFLAGS = -pg -mkl -qopenmp -fast -O3 -xCORE-AVX2
#FFLAGS = -pg -mkl -qopenmp -fast -O3
FFLAGS = -pg -llapack -lblas  -fopenmp -O3
#--compiler=$(C2C)

writeints: integral_writer.f90
	f2py -c --fcompiler=$(F2C)  --f90flags="$(F2FLAGS)" --opt="$(F2LINK_FLAGS)" -m writeints integral_writer.f90
#        f2py -c --fcompiler=$(FC) --f90flags="$(FFLAGS)" --opt="$(LINK_FLAGS)" -m writeints integral_writer.f90

writecube: dens_io.f90
	f2py -c --fcompiler=$(F2C) --compiler=$(C2C) --f90flags="$(F2FLAGS)" --opt="$(F2LINK_FLAGS)"  -m writecube dens_io.f90
#        f2py -c --fcompiler=$(FC) --f90flags="$(FFLAGS)" --opt="$(LINK_FLAGS)"  -m writecube dens_io.f90
#        f2py -c --fcompiler=$(FC) --f90flags="$(FFLAGS)" --opt="$(LINK_FLAGS)" -m writecube dens_io.f90
generate:
	$(FTC) generate_int_v2.0.f90 $(FFLAGS) -o generate_int_v2.0.o
#        $(FTC) write2string.f90 generate_int_v2.0.f90 $(FFLAGS) -o generate_int_v2.0.o

cubemake: mcend_cube_base.f90
	f2py -c --fcompiler=$(F2C) --compiler=$(C2C) --f90flags="$(F2FLAGS)" --opt="$(F2LINK_FLAGS)" -m cubemake mcend_cube_base.f90

#OBJECTS = generate_ints_v2.o
#
#
#generate_ints_v2.0: generate_ints_v2.f90
#        ifort -c generate_ints_v2.f90
#        ifort -o generate_ints_v2.o generate_ints_v2.0
#
#f2pyclean:
#        rm -f writeints*.so
#
clean:
	rm -f *.mod *.o

cleanints:
	rm -f writeints*.so

cleancube:
	rm -f writecube*.so

cleancubemake:
	rm -f cubemake.*so

#generate_ints: generate_ints_v2.f90
	#ifort $(CFLAGS) generate_ints_v2.f90
#CFLAGS = -O3 -c -traceback -i8 -qopenmp -I${MKLROOT}/include
#LINK_FLAGS = -O3 -traceback -i8 -qopenmp -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

#write_h5: generate_ints_v2.f90
#        $(FC) -o $@ write-h5-ints.f90 main-h5.f90 $(LIBSHDF) $(LIB)


#.f90.o: $(FSOURCE)
#        $(FC) $(FLAGS) $?
#
#.SUFFIXES:.o.c.f90
