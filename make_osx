# Makefile for lensview on Mac OSX. Change the make variables according to your needs
# this file is for lensview including the fastell.f (and associated) code
# be sure to "make clean" if something goes wrong

INCLUDES1=-I/Users/cheng109/work/phosim/phosim/source/cfitsio
INCLUDES2=-I/opt/local/include
INCLUDES3=-I/opt/local/fftw-2.1.5/rfftw
INCLUDES4=-I/opt/local/fftw-2.1.5/fftw

LIBS1=-L/Users/cheng109/work/phosim/phosim/source/cfitsio/lib
LIBS2=-L/opt/local/lib
#FC=g77
FC=gfortran

# Regular compiling on both platforms with debugging symbols
#CFLAGS= -g $(CFLAG_OS) $(INCLUDES) $(LIBS)
# Optimised compiling on most platforms
CFLAGS= -O $(INCLUDES1) $(INCLUDES2) $(INCLUDES3) $(INCLUDES4) $(LIBS1) $(LIBS2)

COMMON_HDRS= log.h common.h
OBJS= log.o parseopts.o lv_common.o lv_lens.o lv_geom.o lv_image.o lv_mem.o lv_critline.o lv_paramfit.o lensview.o
lensview:   libfortranstuff.a $(OBJS) $(COMMON_HDRS) libfortranstuff.a
# with fastell
	$(CC) $(CFLAGS) -o lensview $(OBJS) -lgsl -lgslcblas -lcfitsio -ldrfftw -ldfftw -lm -L. -lfortranstuff -lgfortran
#-lg2c

libfortranstuff.a:
	$(FC) -O -c slatec/src/*.f fastell.f
	ar -r libfortranstuff.a *.o
	rm *.o

clean:
	rm -f *.o 
