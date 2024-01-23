CC=c++
CFLAGS=-O2 -std=c++14
INC= 
L= -larmadillo

all: main clean

main: sincore.cpp
	$(CC) $(CFLAGS) -o sincore SinCore.cpp $(L) 

clean:
	rm -f *.o