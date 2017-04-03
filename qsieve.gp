{
  L(N)=exp(sqrt(log(N)*log(log(N))));
}

{
  B(N)=ceil(L(N)^(1/sqrt(2)));
}

{
  f(T,N)=T^2-N;
}

{
qsieve(N)=
  local(numValues, tValues, origValues, currValues, factorBase, numPrimes, curPrime, power, solutions);
  numValues = ceil(L(N)^sqrt(2));
  tValues = vector(numValues, x, floor(sqrt(N)) + x);
  origValues = vector(numValues, x, f(tValues[x], N));
  currValues = vector(numValues, x, f(tValues[x], N));
  factorBase = primes([2,B(N)]);
  numPrimes = #factorBase;
  factorizationTable = matrix(numValues, numPrimes, x, y, 0);
  for(primeIndex = 2, numPrimes,
    curPrime = factorBase[primeIndex];
    if (N % curPrime == 0,
      print(curPrime, " is a factor of ", N);
      return();
    );
    power = 1;
    solutions = lift(polrootspadic(f(T,N), curPrime, power));
    print(currValues);
    while(#solutions > 0 && (solutions[1] < tValues[#tValues] || solutions[2] < tValues[#tValues]),
      currValues = divide(currValues, factorizationTable, (solutions[1] - (tValues[1] % curPrime^power - 1)-1) % curPrime^power + 1, curPrime, curPrime^power, primeIndex);
      currValues = divide(currValues, factorizationTable, (solutions[2] - (tValues[1] % curPrime^power - 1)-1) % curPrime^power + 1, curPrime, curPrime^power, primeIndex);
      power += 1;
      solutions = lift(polrootspadic(f(T,N), curPrime, power));
      print(currValues);
    );
  );
  print(factorizationTable);
}

{
divide(currValues, factorizationTable, solution, curPrime, primePower, primeIndex) =
  while (solution <= #currValues,
    currValues[solution] /= curPrime;
    factorizationTable[solution, primeIndex] += 1;
    solution += primePower;
  );
  return(currValues);
}
