# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-posix.o partdiff-seq.o askparams.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-posix partdiff-seq

partdiff-posix: partdiff-posix.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-posix.o askparams.o $(LIBS)

partdiff-seq: partdiff-seq.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-seq.o askparams.o $(LIBS)

clean:
	$(RM) *.o *~

partdiff-posix.o: partdiff-posix.c Makefile
partdiff-seq.o: partdiff-seq.c Makefile

askparams.o: askparams.c Makefile
