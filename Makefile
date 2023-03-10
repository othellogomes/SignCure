CC=c++
CFLAGS=-O2 -std=c++14
INC= 
L= -larmadillo

all: main clean

main: SignCure.cpp
	$(CC) $(CFLAGS) -o signcure SignCure.cpp $(L) 

clean:
	rm -f *.o
