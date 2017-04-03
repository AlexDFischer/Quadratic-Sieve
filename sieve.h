#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

typedef struct
{
  size_t len;
  mpz_t *nums;
} BigNumList;

typedef size_t prime_t;

typedef struct
{
  size_t len;
  prime_t *primes;
} PrimeList;

typedef unsigned char exponent_t;  // assuming no prime power factors will ahve exponents > 255

typedef struct
{
  size_t len;
  size_t numPrimes;
  exponent_t *mem;
} FactorizationTable;

typedef struct
{
  size_t numRows, numCols;
  size_t *mem;
} Matrix;

/* defined in primes.c */
PrimeList primesLEq(prime_t B);
void freePrimeList(PrimeList primeList);
void printPrimes(PrimeList primeList);

/* defined in functions.c */
prime_t smoothnessBound(mpz_t N);
void lowerBoundT(mpz_t rop, mpz_t N);
void upperBoundT(mpz_t rop, mpz_t N);

BigNumList initList(size_t len);
void freeBigNumList(BigNumList list);
void printBigNumList(BigNumList list);

FactorizationTable initFactorizationTable(size_t len, size_t numPrimes);
void freeFactorizationTable(FactorizationTable table);
void factorizationTableIncrementExponent(FactorizationTable table, size_t fTindex, size_t primeIndex);
exponent_t factorizationTableExponent(FactorizationTable table, size_t fTindex, size_t primeIndex);
void printFactorization(FactorizationTable table, PrimeList primes, size_t ftIndex);

void divideAtInterval(BigNumList values, FactorizationTable table,
  size_t initialIndex, size_t offset, mpz_t divisor, size_t primeIndex);
