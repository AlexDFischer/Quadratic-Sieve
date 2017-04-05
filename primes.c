#include "sieve.h"

/**
 * Returns a list of primes less than B, using the sieve of Erastosthenes
 */
PrimeList primesLEq(prime_t B)
{
  // isPrime[n] is 1 iff n is prime
  char isPrime[B + 1];
  prime_t n;
  for (n = 0; n <= B; n++)
  {
    isPrime[n] = (char) (n % 2); // initialize even numbers to 0
  }
  isPrime[2] = 1;
  prime_t p = 2;
  prime_t sieveBound = (prime_t) sqrt(B) + 1;
  for (p = 3; p <= sieveBound; p += 2)
  {
    if (isPrime[p])
    {
      for (n = 2 * p; n <= B; n += p)
      {
        isPrime[n] = 0;
      }
    }
  }
  PrimeList result;
  prime_t count = 1;
  for (n = 3; n <= B; n += 2)
  {
    count += isPrime[n];
  }
  result.len = count;
  result.primes = malloc(count * sizeof(prime_t));
  result.primes[0] = 2;
  n = 1;
  p = 1;
  while (n < count)
  {
    do
    {
      p += 2;
    } while (!isPrime[p]);
    result.primes[n] = p;
    n++;
  }
  return result;
}

void freePrimeList(PrimeList primes)
{
  free(primes.primes);
}

void printPrimes(PrimeList primes)
{
  int i;
  for (i = 0; i < primes.len; i++)
  {
    printf("%d: %lu\n", i + 1, primes.primes[i]);
  }
}
