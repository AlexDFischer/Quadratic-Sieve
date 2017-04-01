#include "sieve.h"

int main(int argc, char **argv)
{
  mpz_t N, lower, upper;
  mpz_init_set_ui(N, 1234567890);
  mpz_init(lower);
  mpz_init(upper);

  prime_t B = smoothnessBound(N);
  printf("We will generate primes up to %lu\n", B);
  PrimeList factorBase = primesLEq(B);
  size_t numPrimes = factorBase.len;
  printPrimes(factorBase);
  printf("We generated %lu primes\n", numPrimes);
  lowerBoundT(lower, N);
  upperBoundT(upper, N);
  // get number of value's we'll try as a size_t
  mpz_t numOriginalValues;
  mpz_init_set(numOriginalValues, upper);
  mpz_sub(numOriginalValues, numOriginalValues, lower);
  // len is the number of values of f(T) that we'll find
  size_t len = mpz_get_ui(numOriginalValues);
  len++;
  printf("We will try from T=");
  mpz_out_str(stdout, 10, lower);
  printf(" to T=");
  mpz_out_str(stdout, 10, upper);
  printf(".  We will try %lu different values ot T.\n", len);

  BigNumList tValues = initList(len), origValues = initList(len), currValues = initList(len);
  size_t i;

  // set T values
  mpz_set(tValues.nums[0], lower);
  for (i = 1; i < len; i++)
  {
    mpz_add_ui(tValues.nums[i], tValues.nums[i-1], 1);
  }
  // set f(T) values
  for (i = 0; i < len; i++)
  {
    mpz_mul(origValues.nums[i], tValues.nums[i], tValues.nums[i]);
    mpz_sub(origValues.nums[i], origValues.nums[i], N);
    mpz_set(currValues.nums[i], origValues.nums[i]);
  }
  FactorizationTable factorizationTable = initFactorizationTable(len, numPrimes);
  printf("checkpoint 4\n");
  // basically do trial division with p=2
  size_t primeIndex = 0;
  size_t prime = factorBase.primes[primeIndex];
  size_t primePower = prime;
  mpz_t primeMP, primePowerMP, fTmodpe;
  mpz_init_set_ui(primeMP, prime);
  mpz_init_set_ui(primePowerMP, primePower);
  mpz_init(fTmodpe);
  size_t var1 = 0, var2 = 0;

  //printBigNumList(tValues);
  //printBigNumList(origValues);

  // outer loop: loop over powers of our prime 2

  // number of solutions (unique modulo p^e) to the congruence f(T) = 0 mode p^e
  int foundSolutions = 1;
  do
  {
    foundSolutions = 0;
    // FIRST: find solutions to f(T) = 0 mod 2^e
    for (var1 = 0; var1 < primePower && var1 < currValues.len; var1++)
    {
      mpz_mod(fTmodpe, origValues.nums[var1], primePowerMP);
      if (mpz_get_ui(fTmodpe) == 0)
      {
        foundSolutions++;
        // INNER LOOP: now advance by 2^e and divide each entry by 2, incrementing
        // the prime power in the factorization table as well
        for (var2 = var1; var2 < currValues.len; var2 += primePower)
        {
          mpz_fdiv_q(currValues.nums[var2], currValues.nums[var2], primeMP);
          factorizationTableIncrementExponent(factorizationTable, var2, primeIndex);
        }
      }
    }
    if (foundSolutions > 0)
    {
      //printf("For prime power =  %lu: ", primePower);
      //printBigNumList(currValues);
    }
    // now try the next prime power
    primePower *= prime;
    mpz_mul(primePowerMP, primePowerMP, primeMP);
  } while (foundSolutions > 0);

  // now for the more general case: sieve over all of the primes in our factor base
  for (primeIndex = 1; primeIndex < factorBase.len; primeIndex++)
  {
    prime = factorBase.primes[primeIndex];
    primePower = prime;
    mpz_set_ui(primeMP, prime);
    mpz_set_ui(primePowerMP, primePower);
    // middle loop: loop over powers of our prime p
    foundSolutions = 1;
    do
    {
      foundSolutions = 0;
      // FIRST: find solutions to f(T) = 0 mod p^e
      for (var1 = 0; var1 < primePower && var1 < currValues.len; var1++)
      {
        mpz_mod(fTmodpe, origValues.nums[var1], primePowerMP);
        if (mpz_get_ui(fTmodpe) == 0)
        {
          foundSolutions++;
          // INNER LOOP: now advance by p^e and divide each entry by p, incrementing
          // the prime power in the factorization table as well
          for (var2 = var1; var2 < currValues.len; var2 += primePower)
          {
            mpz_fdiv_q(currValues.nums[var2], currValues.nums[var2], primeMP);
            factorizationTableIncrementExponent(factorizationTable, var2, primeIndex);
          }
        }
      }
      if (foundSolutions)
      {
        //printf("For prime power =  %lu: ", primePower);
        //printBigNumList(currValues);
      }
      // now try the next prime power
      primePower *= prime;
      mpz_mul(primePowerMP, primePowerMP, primeMP);
    } while (foundSolutions > 0);
  }

  // find all values in currValues that sieved down to 1
  for (var1 = 0; var1 < currValues.len; var1++)
  {
    if (mpz_get_ui(currValues.nums[var1]) == 1)
    {
      mpz_out_str(stdout, 10, tValues.nums[var1]);
      printf("^2 - N = ");
      printFactorization(factorizationTable, factorBase, var1);
      printf("\n");
    }
  }

  freeFactorizationTable(factorizationTable);
  mpz_clear(primeMP);
  mpz_clear(primePowerMP);
  mpz_clear(fTmodpe);
  freeBigNumList(currValues);
  freeBigNumList(tValues);
  freeBigNumList(origValues);
  freePrimeList(factorBase);
  mpz_clear(N);
  mpz_clear(lower);
  mpz_clear(upper);
  mpz_clear(numOriginalValues);
}
