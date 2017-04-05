main: factor.c primes.c functions.c factorizationtable.c linalg.c main.c
	gcc -c main.c -c factor.c -c primes.c -c functions.c -c factorizationtable.c -c linalg.c -lm -lmpfr -lgmp -Wall -I/usr/include/gsl/include
	gcc -L/usr/include/gsl/lib main.o factor.o primes.o functions.o factorizationtable.o linalg.o -lm -lmpfr -lgmp -o main
