all: pointer hello

#define options for gcc
CFLAGS: -Wall -Wextra -Wpointer-arith 

#define objects that are made from .c files
OBJ: pointer.o hello.o

#make clean a phony target so that it will always be created
.PHONY: clean

pointer: pointer.o
	gcc $(CFLAGS) -o $@ $+

hello: hello.o
	gcc $(CFLAGS) -o $@ $+

#make objects out of the zugehörige .c files	
%.o: %.c
	gcc $(CFLAGS) -c $<

#remove the programm and all objects
clean:
	rm -f pointer
	rm -f hello
	rm -f *.o
	rm -f *~
