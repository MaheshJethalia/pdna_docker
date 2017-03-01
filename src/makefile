# defining necessary variables for compilation
CC=g++
CFLAGS=-c

main: main.o molecule.o
	$(CC) main.o molecule.o -o main

main.o: main.cpp 
	$(CC) $(CFLAGS) main.cpp -o main.o

molecule.o: ./libs/molecule.cpp ./libs/molecule.h
	$(CC) $(CFLAGS) ./libs/molecule.cpp -o molecule.o

clean:
	rm -rf *.o main