# Location of HDF5 binaries (with include/ and lib/ underneath)
HDF5 = /usr/local/hdf5-1.10.1
#HDF5 = /home/lucas/local
# cruntch3
#HDF5 = /export/apps/UTILS/hdf5/1.8.20/intel/17.0
# Location of External Libraries
#LIBZ = /usr/lib/x86_64-linux-gnu/libz.a
#LIBSZ = /usr/lib/x86_64-linux-gnu/libsz.a
#LIBZ    = /mnt/hdf/packages/hdf5/v189/Linux_2.6/hdf5-1.8.9-linux-static/lib/libz.a
#LIBSZ   = /mnt/hdf/packages/hdf5/v189/Linux_2.6/hdf5-1.8.9-linux-static/lib/libsz.a

# Compiler
#FC        = gfortran
FC        = ifort

#FSOURCE = write_h5 write_h5_debug
FSOURCE = write_h5
#OBJ = utils/simple-h5-reader.o main-h5.o
OBJ = simple-h5-reader.o main-h5.o
HDF5_LIB = ${HDF5}/lib
CFLAGS     = -O3  -c -traceback -I$(HDF5)/include -pg
#LINK_FLAGS = -O3 -traceback $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -lrt -lm -pg
LINK_FLAGS = -O3 -traceback $(HDF5_LIB)/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -lrt -lm -pg

#CFLAGS_debug     = -O0  -c -traceback -I$(HDF5)/include -pg
#LINK_FLAGS_debug = -O0 -traceback $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -lrt -lm -pg
#LINK_FLAGS_debug = -O0 -traceback $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -lrt -lm -pg

all:    $(FSOURCE)

#write_h5: main-h5.f90 simple_h5_io.f90
#        $(FC) -o $@ simple_h5_io.f90 main-h5.f90 $(LIBSHDF) $(LIB)
write_h5: $(OBJ)
	$(FC) -p -o $(FSOURCE) $(OBJ) $(LINK_FLAGS)

#write_h5_debug: $(OBJ)
#        $(FC) -p -o $(FSOURCE) $(OBJ) $(LINK_FLAGS_debug)
clean:
	rm -f *.o *.mod *.h5
#        cd utils; rm -f *.o *.mod *.h5
#        cd ..
%.o: %.f90
	$(FC) $(CFLAGS)  $<

#.SUFFIXES:.o.c.f90
