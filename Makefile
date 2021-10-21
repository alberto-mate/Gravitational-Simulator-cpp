output: sim-aos.o
	g++ -o sim-aos -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors sim-aos.o
	g++ -o sim-soa -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors sim-soa.o

sim-aos.o: sim-aos.cpp
	g++ -c sim-aos.cpp

sim-soa.o: sim-soa.cpp
	g++ -c sim-soa.cpp

clean:
	rm *.o 