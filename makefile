#COMPILERFLAGS = -Wall -DNOTGRAPH #-O3
COMPILERFLAGS = -w #-DNOTGRAPH #-O3
#COMPILERFLAGS = -DNOTGRAPH -O3
CC = gcc

CLAPACK_ROOTPATH = /home/rangan/CLAPACK
CLAPACK_F2CDIR  = $(CLAPACK_ROOTPATH)/F2CLIBS
INCLUDE = -I/usr/local/include -I/home/rangan/local_sphinx/include/ -I$(CLAPACK_ROOTPATH)/SRC -I$(CLAPACK_ROOTPATH)
LIBDIR  = -L/home/rangan/local_sphinx/lib
LIBRARIES = -lfftw3 $(CLAPACK_ROOTPATH)/lapack_LINUX.a $(CLAPACK_ROOTPATH)/blas_LINUX.a $(CLAPACK_F2CDIR)/libF77.a $(CLAPACK_F2CDIR)/libI77.a -lm -lglut -lGL -lGLU -I/usr/include/GL

CFLAGS = $(COMPILERFLAGS) $(INCLUDE)

d_v1: d_v1.o
	$(CC) $(CFLAGS) -o $@.out $(LIBDIR) $< $(LIBRARIES)
d_v1.o : d_lgn.c d_lgn.h d_cortex.c d_cortex.h d_llists.c d_llists.h d_ogletools.c d_ogletools.h d_inputhandling.c d_inputhandling.h d_datahandling.c d_datahandling.h

al: al.o
	$(CC) $(CFLAGS) -o $@.out $(LIBDIR) $< $(LIBRARIES)
al.o : matrices.c matrices.h al_cortex.c al_cortex.h d_llists.c d_llists.h al_ogletools.c al_ogletools.h al_inputhandling.c al_inputhandling.h al_datahandling.c al_datahandling.h

inputraracreate: inputraracreate.o
	$(CC) $(CFLAGS) -o $@.out $(LIBDIR) $< $(LIBRARIES)
inputraracreate.o : inputraracreate.c

clean : 
	rm d_v1.o
	rm al.o
