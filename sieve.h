#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

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

void factor(mpz_t N, double numRelationsPower, double smoothnessPower);

/* defined in primes.c */
PrimeList primesLEq(prime_t B);
void freePrimeList(PrimeList primeList);
void printPrimes(PrimeList primeList);

/* defined in functions.c */
prime_t LofNpow(mpz_t N, double pow);
void lowerBoundT(mpz_t rop, mpz_t N);

BigNumList initList(size_t len);
void freeBigNumList(BigNumList list);
void printBigNumList(BigNumList list);

FactorizationTable initFactorizationTable(size_t len, size_t numPrimes);
void freeFactorizationTable(FactorizationTable table);
void factorizationTableIncrementExponent(FactorizationTable table, size_t fTindex, size_t primeIndex);
exponent_t factorizationTableExponent(FactorizationTable table, size_t fTindex, size_t primeIndex);
void printFactorization(FactorizationTable table, PrimeList primes, size_t ftIndex);

void divideAtInterval(BigNumList values,
  size_t initialIndex, size_t offset, mpz_t divisor);

int findSolutions(size_t startIndex, size_t *sol1index,
  size_t *sol2index, size_t prime, mpz_t primeMP, size_t primePower, mpz_t primePowerMP,
  BigNumList origValues, BigNumList tValues, mpz_t dummy);

/* defined in linalg.c */
Matrix initMatrix(size_t rows, size_t cols);
void set1(Matrix matrix, size_t row, size_t col);
int get(Matrix matrix, size_t row, size_t col);
void printMatrix(Matrix m);
