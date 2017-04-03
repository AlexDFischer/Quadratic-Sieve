#include "sieve.h"

int main(int argc, char **argv)
{
  mpz_t N, lower, upper, dummy;
  mpz_init(N);
  gmp_sscanf(argc > 1 ? argv[1] : "12345", "%Zd", N);
  mpz_init(lower);
  mpz_init(upper);
  mpz_init(dummy);

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

  FactorizationTable factorizationTable = initFactorizationTable(len, numPrimes);
  size_t primeIndex = 0;
  size_t prime = factorBase.primes[primeIndex];
  size_t primePower = prime;
  mpz_t primeMP, primePowerMP, fTmodpe;
  mpz_init_set_ui(primeMP, prime);
  mpz_init_set_ui(primePowerMP, primePower);
  mpz_init(fTmodpe);
  size_t var1 = 0, var2 = 0;
  // first do a little trial division
  for (primeIndex = 0; primeIndex < factorBase.len; primeIndex++)
  {
    //printf("prime %lu: %lu\n", primeIndex, factorBase.primes[primeIndex]);
    if (mpz_fdiv_r_ui(dummy, N, factorBase.primes[primeIndex]) == 0)
    {
      printf("N is divisible by %lu\n", factorBase.primes[primeIndex]);
      exit(0);
    } else
    {
      //printf("N is not divisible by %lu: remainder is %lu\n", factorBase.primes[primeIndex], mpz_fdiv_r_ui(dummy, N, factorBase.primes[primeIndex]));
    }
  }

  printf("We will try from T=");
  mpz_out_str(stdout, 10, lower);
  printf(" to T=");
  mpz_out_str(stdout, 10, upper);
  printf(".  We will try %lu different values of T.\n", len);
  BigNumList tValues = initList(len), origValues = initList(len), currValues = initList(len);
  //size_t i;
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

  //printBigNumList(tValues);
  //printBigNumList(origValues);

  // outer loop: loop over powers of our prime 2

  // number of solutions (unique modulo 2^e) to the congruence f(T) = 0 mode 2^e
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
      foundSolutions = 0;
      // lift our previous solutions mod p^(e-1) to a solution mod p^e
      for (var1 = prevSol1index; var1 <= (prevSol1index + 1) / 2 + primePower && var1 < currValues.len  && !foundSolutions; var1 += primePower / prime)
      {
        mpz_mod(fTmodpe, origValues.nums[var1], primePowerMP);
        if (mpz_get_ui(fTmodpe) == 0)
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
          /*
          gmp_printf("for prime power = %Zd, solutions are T = %Zd", primePowerMP, tValues.nums[sol1index]);
          if (sol2index < tValues.len)
          {
            gmp_printf(" and %Zd", tValues.nums[sol2index]);
          }
          printf("\n");
          */
          divideAtInterval(currValues, factorizationTable, sol1index, primePower, primeMP, primeIndex, tValues);
          divideAtInterval(currValues, factorizationTable, sol2index, primePower, primeMP, primeIndex, tValues);
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
          mpz_mod(fTmodpe, origValues.nums[var1], primePowerMP);
          if (mpz_get_ui(fTmodpe) == 0)
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
            /*
            gmp_printf("for prime power = %Zd, solutions are T = %Zd", primePowerMP, tValues.nums[sol1index]);
            if (sol2index < tValues.len)
            {
              gmp_printf(" and %Zd", tValues.nums[sol2index]);
            }
            printf("\n");
            */
            divideAtInterval(currValues, factorizationTable, sol1index, primePower, primeMP, primeIndex, tValues);
            divideAtInterval(currValues, factorizationTable, sol2index, primePower, primeMP, primeIndex, tValues);
            //printf("for prime power=%lu", primePower);
            //printBigNumList(currValues);
            prevSol1index = sol1index;
            prevSol2index = sol2index;
            //printf("did division: p^e is %lu, sol1index = %lu, sol2index = %lu\n", primePower, sol1index, sol2index);
            //printBigNumList(currValues);
          }
        }
        /*
        if (!foundSolutions)
        {
          printf("did not find any solutions for primepower=%lu\n", primePower);
        }
        */
      }
      /*
      // FIRST: find solutions to f(T) = 0 mod p^e
      for (var1 = 0; var1 < primePower && var1 < currValues.len; var1++)
      {
        mpz_mod(fTmodpe, origValues.nums[var1], primePowerMP);
        if (mpz_get_ui(fTmodpe) == 0)
        {
          foundSolutions++;
          divideAtInterval(currValues, factorizationTable, var1, primePower, primeMP, primeIndex);
          // INNER LOOP: now advance by p^e and divide each entry by p, incrementing
          // the prime power in the factorization table as well

        }
      }
      */
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
      mpz_out_str(stdout, 10, tValues.nums[relationsList[var2]]);
      printf("^2 - N = ");
      printFactorization(factorizationTable, factorBase, relationsList[var2]);
      printf("\n");
      var2++;
    }
  }
  /*
  for (var1 = 0; var1 < tValues.len; var1++)
  {
    printf("%lu:", var1);
    mpz_out_str(stdout, 10, tValues.nums[var1]);
    printf(":");
    mpz_out_str(stdout, 10, currValues.nums[var1]);
    printf("\n");
  }
  */
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
