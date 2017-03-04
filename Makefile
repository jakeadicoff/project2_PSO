CC = g++ -O3
CCFLAGS = -std=c++11 -Wall

all: PSO

debug :
	g++ -std=c++11 -g main.cpp PSO.cpp -o debug.o

PSO : main.o pso.o
	$(CC) -o $@ main.o pso.o

main.o : main.cpp main.h
	$(CC) -c $(CCFLAGS) main.cpp -o $@

pso.o : PSO.cpp PSO.h
	$(CC) -c $(CCFLAGS) PSO.cpp -o $@

clean:
	rm PSO
	rm *.o
