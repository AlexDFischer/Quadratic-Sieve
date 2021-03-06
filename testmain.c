#include "sieve.h"

/**
 * Attempts to factor N, checking L(N)^numRelationsPower random numbers and L(N)^smoothnessPower as the smoothness bound.
 */
void factor(mpz_t N, double numRelationsPower, double smoothnessPower)
{
  mpz_t lower, upper, dummy, remainder;
  mpz_init(lower);
  mpz_init(upper);
  mpz_init(dummy); // this variable is used as a cartch-all bignum to avoid having to repeatedly re-initialize new ones
  mpz_init(remainder);

  prime_t B = LofNpow(N, smoothnessPower);
  printf("We will generate primes up to %lu\n", B);
  PrimeList factorBase = primesLEq(B);
  size_t numPrimes = factorBase.len;
  printPrimes(factorBase);
  printf("We generated %lu primes\n", numPrimes);
  lowerBoundT(lower, N);
  lowerBoundT(upper, N);
  mpz_add_ui(upper, upper, LofNpow(N, numRelationsPower));

  // get number of value's we'll try as a size_t
  mpz_t numOriginalValues;
  mpz_init_set(numOriginalValues, upper);
  mpz_sub(numOriginalValues, numOriginalValues, lower);
  // len is the number of values of f(T) that we'll test
  size_t len = mpz_get_ui(numOriginalValues);
  len++;
  gmp_printf("We will try values of T from %Zd to %Zd (%lu different values).\n", lower, upper, len);
  size_t primeIndex = 0;
  size_t var1 = 0, var2 = 0;
  // first do a little trial division
  for (primeIndex = 0; primeIndex < factorBase.len; primeIndex++)
  {
    //printf("prime %lu: %lu\n", primeIndex, factorBase.primes[primeIndex]);
    if (mpz_fdiv_r_ui(dummy, N, factorBase.primes[primeIndex]) == 0)
    {
      printf("N is divisible by %lu\n", factorBase.primes[primeIndex]);
      exit(0);
    }
  }

  BigNumList tValues = initList(len), origValues = initList(len), currValues = initList(len);
  // set T values
  mpz_set(tValues.nums[0], lower);
  for (var1 = 1; var1 < len; var1++)
  {
    mpz_add_ui(tValues.nums[var1], tValues.nums[var1-1], 1);
  }
  // set f(T) values
  for (var1 = 0; var1 < len; var1++)
  {
    mpz_mul(origValues.nums[var1], tValues.nums[var1], tValues.nums[var1]);
    mpz_sub(origValues.nums[var1], origValues.nums[var1], N);
    mpz_set(currValues.nums[var1], origValues.nums[var1]);
  }

  // first try with prime p = 2 because it's a special case
  // outer loop: loop over powers of our prime 2
  size_t prime = 2;
  size_t primePower = prime;
  mpz_t primeMP, primePowerMP;
  mpz_init_set_ui(primeMP, prime);
  mpz_init_set_ui(primePowerMP, primePower);
  int foundSolutions = 1;
  do
  {
    foundSolutions = 0;
    // FIRST: find solutions to f(T) = 0 mod 2^e
    for (var1 = 0; var1 < primePower && var1 < currValues.len; var1++)
    {
      if (mpz_divisible_p(origValues.nums[var1], primePowerMP))
      {
        foundSolutions++;
        divideAtInterval(currValues, var1, primePower, primeMP);
      }
    }
    // now try the next prime power
    primePower *= prime;
    mpz_mul(primePowerMP, primePowerMP, primeMP);
  } while (foundSolutions > 0);

  // now for the more general case: sieve over all of the primes in our factor base
  size_t sol1index, sol2index, prevSol1index, prevSol2index;
  for (primeIndex = 1; primeIndex < factorBase.len; primeIndex++)
  {
    prime = factorBase.primes[primeIndex];
    primePower = prime;
    mpz_set_ui(primeMP, prime);
    mpz_set_ui(primePowerMP, primePower);

    // middle loop: loop over powers of our prime p
    foundSolutions = 1;
    prevSol1index = 0;
    prevSol2index = 0;
    do
    {
      foundSolutions = findSolutions(prevSol1index, &sol1index, &sol2index, prime, primeMP, primePower, primePowerMP, origValues, tValues, dummy);
      if (!foundSolutions)
      {
        foundSolutions = findSolutions(prevSol2index, &sol1index, &sol2index, prime, primeMP, primePower, primePowerMP, origValues, tValues, dummy);
      }
      if (foundSolutions)
      {
        divideAtInterval(currValues, sol1index, primePower, primeMP);
        divideAtInterval(currValues, sol2index, primePower, primeMP);
        prevSol1index = sol1index;
        prevSol2index = sol2index;
      }

      // lift our previous solutions mod p^(e-1) to a solution mod p^e
      /*for (var1 = prevSol1index; var1 <= (prevSol1index + 1) / 2 + primePower && var1 < currValues.len  && !foundSolutions; var1 += primePower / prime)
      {
        if (mpz_divisible_p(origValues.nums[var1], primePowerMP))
        {
          foundSolutions = 1;
          // assuming N is not a multiple of p, T^2-N=0 mod p will have 2 unique solutions
          sol1index = var1;
          while (sol1index >= primePower)
          {
            sol1index -= primePower;
          }
          mpz_set_ui(dummy, primePower);
          mpz_sub(dummy, dummy, tValues.nums[sol1index]);
          mpz_sub(dummy, dummy, tValues.nums[sol1index]);
          mpz_mod(dummy, dummy, primePowerMP); // if T is a solution to f(T=0) mod p then dummy = p - 2T
          sol2index = sol1index + mpz_get_ui(dummy); // if f(T) = 0 mod p then f(p-T)=(p-T)^2-N=T^2-N=0 mod p
          while (sol2index >= primePower)
          {
            sol2index -= primePower;
          }
          divideAtInterval(currValues, sol1index, primePower, primeMP);
          divideAtInterval(currValues, sol2index, primePower, primeMP);
          //printf("for prime power=%lu", primePower);
          //printBigNumList(currValues);
          prevSol1index = sol1index;
          prevSol2index = sol2index;
          //printf("did division: p^e is %lu, sol1index = %lu, sol2index = %lu\n", primePower, sol1index, sol2index);
          //printBigNumList(currValues);
        }
      }

      if (!foundSolutions)
      {
        // try again starting at prevSol2index. TODO this could be much more modular
        for (var1 = prevSol2index; var1 <= (prevSol2index + 1) / 2 + primePower && var1 < currValues.len  && !foundSolutions; var1 += primePower / prime)
        {
          if (mpz_divisible_p(origValues.nums[var1], primePowerMP))
          {
            foundSolutions = 1;
            // assuming N is not a multiple of p, T^2-N=0 mod p will have 2 unique solutions
            sol1index = var1;
            while (sol1index >= primePower)
            {
              sol1index -= primePower;
            }
            mpz_set_ui(dummy, primePower);
            mpz_sub(dummy, dummy, tValues.nums[sol1index]);
            mpz_sub(dummy, dummy, tValues.nums[sol1index]);
            mpz_mod(dummy, dummy, primePowerMP); // if T is a solution to f(T=0) mod p then dummy = p - 2T
            sol2index = sol1index + mpz_get_ui(dummy); // if f(T) = 0 mod p then f(p-T)=(p-T)^2-N=T^2-N=0 mod p
            while (sol2index >= primePower)
            {
              sol2index -= primePower;
            }
            divideAtInterval(currValues, sol1index, primePower, primeMP);
            divideAtInterval(currValues, sol2index, primePower, primeMP);
            prevSol1index = sol1index;
            prevSol2index = sol2index;
          }
        }
      }*/
      // now try the next prime power
      primePower *= prime;
      mpz_mul(primePowerMP, primePowerMP, primeMP);
    } while (foundSolutions);
  }
  // find all values in currValues that sieved down to 1
  size_t numRelations = 0;
  for (var1 = 0; var1 < currValues.len; var1++)
  {
    if (mpz_get_ui(currValues.nums[var1]) == 1)
    {
      numRelations++;
    }
  }
  size_t relationsList[numRelations];
  var2 = 0;
  for (var1 = 0; var1 < currValues.len; var1++)
  {
    if (mpz_get_ui(currValues.nums[var1]) == 1)
    {
      relationsList[var2] = var1;
      var2++;
    }
  }
  Matrix matrix = initMatrix(numPrimes, numRelations);
  size_t relationIndex, power;
  for (relationIndex = 0; relationIndex < numRelations; relationIndex++)
  {
    gmp_printf("T=%Zd", tValues.nums[relationsList[relationIndex]]);
    for (primeIndex = 0; primeIndex < numPrimes; primeIndex++)
    {
      mpz_set(dummy, origValues.nums[relationsList[relationIndex]]);
      prime = factorBase.primes[primeIndex];
      mpz_set_ui(primeMP, prime);
      power = 0;
      // determine power of prime in prime factorization of origValues[relationsList[relationIndex]]
      mpz_fdiv_qr(dummy, remainder, dummy, primeMP);
      while (mpz_get_ui(remainder) == 0)
      {
        power++;
        mpz_fdiv_qr(dummy, remainder, dummy, primeMP);
      }
      if (power % 2 == 1)
      {
        set1(matrix, primeIndex, relationIndex);
      }
    }
    printf("\n");
  }
  printf("matrix representing the relations we found:\n");
  printMatrix(matrix);

  mpz_clear(remainder);
  mpz_clear(primeMP);
  mpz_clear(primePowerMP);
  freeBigNumList(currValues);
  freeBigNumList(tValues);
  freeBigNumList(origValues);
  freePrimeList(factorBase);
  mpz_clear(lower);
  mpz_clear(upper);
  mpz_clear(numOriginalValues);
}
