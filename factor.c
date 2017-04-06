#include "sieve.h"

/**
 * Attempts to factor N, checking L(N)^numRelationsPower random numbers and L(N)^smoothnessPower as the smoothness bound.
 */
void factor(mpz_t N, double numRelationsPower, double smoothnessPower)
{
  struct timeval t1, t2, t3, t4;
  gettimeofday(&t1, NULL);
  mpz_t lower, upper, dummy, remainder;
  mpz_init(lower);
  mpz_init(upper);
  mpz_init(dummy); // this variable is used as a catch-all bignum to avoid having to repeatedly re-initialize new ones
  mpz_init(remainder);

  prime_t B = LofNpow(N, smoothnessPower);
  PrimeList factorBase = primesLEq(B);
  size_t numPrimes = factorBase.len;
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
  //gmp_printf("We will try values of T from %Zd to %Zd (%lu different values).\n", lower, upper, len);
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
      // searches for solutions starting at prevSol1index, using hensel's lemma
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
      // now try the next prime power
      primePower *= prime;
      mpz_mul(primePowerMP, primePowerMP, primeMP);
    } while (foundSolutions);
  }
  // find all values in currValues that sieved down to 1 and store their indicies
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
  printf("done with relation building\n");
  gettimeofday(&t2, NULL);
  Matrix matrix = initMatrix(numPrimes, numRelations);
  size_t relationIndex, power;
  for (relationIndex = 0; relationIndex < numRelations; relationIndex++)
  {
    //gmp_printf("T=%Zd", tValues.nums[relationsList[relationIndex]]);
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
  }
  Matrix basis = kernelBasis(matrix);
  printf("done with linear algebra\n");
  gettimeofday(&t3, NULL);
  mpz_t a, b, aminusb, gcd;
  mpz_init(a);
  mpz_init(b);
  mpz_init(aminusb);
  mpz_init(gcd);
  // outer loop: iterate through each basis vector of the null space we found;
  for (var1 = 0; var1 < basis.numRows; var1++)
  {
    mpz_set_ui(a, 1);
    mpz_set_ui(b, 1);
    // inner loop: iterate through each element of the null space vector, using it to construct 2 squares congruent mod N
    for (var2 = 0; var2 < numRelations; var2++)
    {
      if (get(basis, var1, var2))
      {
        mpz_mul(a, a, origValues.nums[relationsList[var2]]);
        mpz_mul(b, b, tValues.nums[relationsList[var2]]);
      }
    }
    if (mpz_root(a, a, 2))
    mpz_sub(aminusb, a, b);
    mpz_gcd(gcd, aminusb, N);
    if (mpz_cmp(gcd, N) && mpz_cmp_ui(gcd, 1))
    {
      gmp_printf("nontrivial factor found: %Zd\n", gcd);
      gettimeofday(&t4, NULL);
      printf("relation building took: %ld ms.\n", timeDifferenceMS(&t2, &t1));
      printf("linear algebra took:    %ld ms.\n", timeDifferenceMS(&t3, &t2));
      printf("gcd computation took:   %ld ms.\n", timeDifferenceMS(&t4, &t3));
      printf("total time was:         %ld ms.\n", timeDifferenceMS(&t4, &t1));
      exit(0);
    }
  }
  printf("no nontrivial factors found :( try more relations\n");

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(aminusb);
  mpz_clear(remainder);
  mpz_clear(gcd);
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

void freeStuff()
{

}

void divideAtInterval(BigNumList values, size_t initialIndex, size_t offset, mpz_t divisor)
{
  size_t i;
  for (i = initialIndex; i < values.len; i += offset)
  {
    mpz_fdiv_q(values.nums[i], values.nums[i], divisor);
  }
}

int findSolutions(size_t startIndex, size_t *sol1index,
  size_t *sol2index, size_t prime, mpz_t primeMP, size_t primePower, mpz_t primePowerMP,
  BigNumList origValues, BigNumList tValues, mpz_t dummy)
{
  size_t var1;
  // start at startIndex and advance by p^(e-1), which is valid by hensel's lemma
  for (var1 = startIndex; var1 <= (startIndex + 1) / 2 + primePower && var1 < tValues.len; var1 += primePower / prime)
  {
    if (mpz_divisible_p(origValues.nums[var1], primePowerMP))
    {
      // assuming N is not a multiple of p, T^2-N=0 mod p will have 2 unique solutions
      *sol1index = var1;
      while (*sol1index >= primePower)
      {
        *sol1index -= primePower;
      }
      // if T1 is first solution, second solution is -T1 (mod p^e)
      mpz_set_ui(dummy, primePower);
      mpz_sub(dummy, dummy, tValues.nums[*sol1index]);
      mpz_sub(dummy, dummy, tValues.nums[*sol1index]);
      mpz_mod(dummy, dummy, primePowerMP);
      *sol2index = *sol1index + mpz_get_ui(dummy);
      while (*sol2index >= primePower)
      {
        *sol2index -= primePower;
      }
      return 1;
    }
  }
  return 0;
}

long int timeDifferenceMS(struct timeval *t1, struct timeval *t0)
{
  return (t1->tv_sec - t0->tv_sec) * 1000 + (t1->tv_usec - t0->tv_usec) / 1000;
}
