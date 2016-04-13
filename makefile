#
CC=g++

#
CFLAGS=-std=c++0x -pthread
#
LDFLAGS=-fopenmp


all: test

test: test.o
	$(CC) $(LDFLAGS) test.o -o test

test.o: test.cpp
	$(CC) $(CFLAGS) test.cpp -o test.o 

clean:
	rm -rf *.o test
