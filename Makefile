all:
	g++ -Wall -O3 -o MAIN main.cpp

clean:
	rm -f main *.o
