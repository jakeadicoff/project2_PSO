CC = g++ -O3
CFLAGS = -std=c++11 -Wall

PSO : main.o pso.o
	$(CC) -o $@ main.o pso.o

main.o : main.cpp main.h
	$(CC) -c $(CFLAGS) main.cpp -o $@

pso.o : PSO.cpp PSO.h
	$(CC) -c $(CFLAGS) PSO.cpp -o $@

clean:
	rm PSO
	rm *.o