#include "sieve.h"

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
