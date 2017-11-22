new:
	make clean
	make all

all:
	gcc -c -I/usr/local/include -I. -std=c++11 *.cpp
	ar -cvq libHEAAN.a *.o

clean:
	rm *.o
	rm *.a
