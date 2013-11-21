# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb 
LFLAGS = $(CFLAGS)
LIBS   = -lm

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: timempi timempi2

timempi: timempi.o Makefile
	$(CC) $(LFLAGS) -o $@ timempi.o $(LIBS)

timempi2: timempi2.o Makefile
	$(CC) $(LFLAGS) -o $@ timempi2.o $(LIBS)

clean:
	$(RM) *.o *~

timempi.o: timempi.c Makefile
timempi2.o: timempi2.c Makefile

