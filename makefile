testmain: testmain.c primes.c functions.c factorizationtable.c
	gcc -c testmain.c -c primes.c -c functions.c -c factorizationtable.c -lm -lmpfr -lgmp -lgsl -lgslcblas -Wall -I/usr/include/gsl/include
	gcc -L/usr/include/gsl/lib testmain.o primes.o functions.o factorizationtable.o -lm -lmpfr -lgmp -lgsl -lgslcblas -o testmain
