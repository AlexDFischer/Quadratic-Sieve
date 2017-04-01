#include "sieve.h"
#define FLOAT_PRECISION 128

/*
 * Sets the value of rop to L(N)
 */
void L(mpfr_t rop, mpz_t N)
{
  mpfr_t Nf, lnN, lnlnN, lnN_lnlnN, sqrt_lnN_lnlnN;
  mpfr_init2(Nf, FLOAT_PRECISION);
  mpfr_init2(lnN, FLOAT_PRECISION);
  mpfr_init2(lnlnN, FLOAT_PRECISION);
  mpfr_init2(lnN_lnlnN, FLOAT_PRECISION);
  mpfr_init2(sqrt_lnN_lnlnN, FLOAT_PRECISION);
  mpfr_set_z(Nf, N, MPFR_RNDN);
  mpfr_log(lnN, Nf, MPFR_RNDN);
  mpfr_log(lnlnN, lnN, MPFR_RNDN);
  mpfr_mul(lnN_lnlnN, lnN, lnlnN, MPFR_RNDN);
  mpfr_sqrt(sqrt_lnN_lnlnN, lnN_lnlnN, MPFR_RNDN);
  mpfr_exp(rop, sqrt_lnN_lnlnN, MPFR_RNDN);
  mpfr_clear(Nf);
  mpfr_clear(lnN);
  mpfr_clear(lnlnN);
  mpfr_clear(lnN_lnlnN);
  mpfr_clear(sqrt_lnN_lnlnN);
}

/**
 * Returns L(N)^(1/sqrt(2)) as an unsigned long or prime_t
 */
prime_t smoothnessBound(mpz_t N)
{
  mpfr_t l, sqrt2inverse, b;
  mpfr_init2(l, FLOAT_PRECISION);
  mpfr_init2(sqrt2inverse, FLOAT_PRECISION);
  mpfr_init2(b, FLOAT_PRECISION);
  L(l, N);
  mpfr_set_d(sqrt2inverse, 1.0/sqrt(2.0), MPFR_RNDN);
  mpfr_pow(b, l, sqrt2inverse, MPFR_RNDN);
  prime_t result = mpfr_get_ui(b, MPFR_RNDU);
  mpfr_clear(l);
  mpfr_clear(b);
  mpfr_clear(sqrt2inverse);
  return result;
}
/**
 * floor(sqrt(N))+1
 */
void lowerBoundT(mpz_t rop, mpz_t N)
{
  mpfr_t sqrtN;
  mpfr_init_set_z(sqrtN, N, MPFR_RNDN);
  mpfr_sqrt(sqrtN, sqrtN, MPFR_RNDN);
  mpfr_get_z(rop, sqrtN, MPFR_RNDD);
  mpz_add_ui(rop, rop, 1);
  mpfr_clear(sqrtN);
}

void upperBoundT(mpz_t rop, mpz_t N)
{
  mpz_t lower;
  mpfr_t l, sqrt2, lToSqrt2;

  mpz_init(lower);
  mpfr_init2(l, FLOAT_PRECISION);
  mpfr_init2(sqrt2, FLOAT_PRECISION);
  mpfr_init2(lToSqrt2, FLOAT_PRECISION);

  lowerBoundT(lower, N);
  L(l, N);
  mpfr_set_d(sqrt2, sqrt(2.0), MPFR_RNDN);
  mpfr_pow(lToSqrt2, l, sqrt2, MPFR_RNDN);
  mpfr_get_z(rop, lToSqrt2, MPFR_RNDU);
  mpz_add(rop, rop, lower);

  mpz_clear(lower);
  mpfr_clear(l);
  mpfr_clear(sqrt2);
  mpfr_clear(lToSqrt2);
}

/**
 * Creates a list of mpz_t elements with len number elements and does the initialization of the mpz_t's
 */
BigNumList initList(size_t len)
{
  BigNumList list;
  list.len = len;
  list.nums = malloc(len * sizeof(mpz_t));
  size_t i;
  for (i = 0; i < len; i++)
  {
    mpz_init(list.nums[i]);
  }
  return list;
}

void freeBigNumList(BigNumList list)
{
  size_t i;
  for (i = 0; i < list.len; i++)
  {
    mpz_clear(list.nums[i]);
  }
  free(list.nums);
}

void printBigNumList(BigNumList list)
{
  size_t i;
  for (i = 0; i < list.len; i++)
  {
    mpz_out_str(stdout, 10, list.nums[i]);
    printf(" ");
  }
  printf("\n");
}

FactorizationTable initFactorizationTable(size_t len, size_t numPrimes)
{
  FactorizationTable result;
  result.mem = calloc(len * numPrimes, sizeof(exponent_t));
  result.len = len;
  result.numPrimes = numPrimes;
  return result;
}

void freeFactorizationTable(FactorizationTable table)
{
  free(table.mem);
}

void factorizationTableIncrementExponent(FactorizationTable table, size_t fTindex, size_t primeIndex)
{
  table.mem[fTindex * table.numPrimes + primeIndex]++;
}

exponent_t factorizationTableExponent(FactorizationTable table, size_t fTindex, size_t primeIndex)
{
  return table.mem[fTindex * table.numPrimes + primeIndex];
}

void printFactorization(FactorizationTable table, PrimeList primes, size_t ftIndex)
{
  exponent_t *mem = table.mem + ftIndex * table.numPrimes;
  size_t primeIndex, firstFactor = 1;
  for (primeIndex = 0; primeIndex < table.numPrimes; primeIndex++)
  {
    if (mem[primeIndex] > 0)
    {
      if (firstFactor)
      {
        firstFactor = 0;
      } else
      {
        printf(" * ");
      }
      printf("%lu^%d", primes.primes[primeIndex], mem[primeIndex]);
    }
  }
}
