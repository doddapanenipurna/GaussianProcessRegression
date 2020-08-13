gpr: main.o generators.o gpr.o
	g++ -fopenmp -g -Wall main.o generators.o gpr.o -o gpr
main.o: main.cpp generators.h gpr.h
	g++ -fopenmp -g -Wall -c main.cpp
generators.o: generators.h generators.cpp
	g++ -fopenmp -g -Wall -c generators.cpp
gpr.o: gpr.h gpr.cpp
	g++ -fopenmp -g -Wall -c gpr.cpp
clean:
	rm *.o gpr
