# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-par.o askparams.o displaymatrix-mpi.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-par

partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) *.o *~

partdiff-seq.o: partdiff-par.c Makefile

askparams.o: askparams.c Makefile

displaymatrix-mpi.o: displaymatrix-mpi.c Makefile
