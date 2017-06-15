# available c compiler for c files
CC = gcc

# activate compilation warnings and debugging information
CFLAGS = -g -Wall

# linker flags for Unix systems (for math.h and fftw3.h)
LINKERFLAGS = -lm -lfftw3

# header files
HEADERS = config.h docker.h

# clean up command
RM = rm

all: main

main: main.o preprocessing.o docker.o scoring.o
	$(CC) $(CFLAGS) preprocessing.c docker.c scoring.c main.c -o main $(LINKERFLAGS)

main.o: main.c $(HEADERS)
	$(CC) $(CFLAGS) -c main.c -o main.o

preprocessing.o: preprocessing.c $(HEADERS)
	$(CC) $(CFLAGS) -c preprocessing.c -o preprocessing.o

docker.o: docker.c $(HEADERS)
	$(CC) $(CFLAGS) -c docker.c -o docker.o

scoring.o: scoring.c $(HEADERS)
	$(CC) $(CFLAGS) -c scoring.c -o scoring.o

clean:
	$(RM) main *.o 
