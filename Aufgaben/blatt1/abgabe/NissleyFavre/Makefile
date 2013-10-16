all: pointer helloworld
	
pointer: objectfiles
	gcc -o pointer.x pointer.o

helloworld: objectfiles
	gcc -o helloworld.x helloworld.o
	
objectfiles:
	gcc -c $(wildcard *.c)

clean:
	rm -f *.x
	rm -f *.o
	rm -f *~
