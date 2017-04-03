qsieve(N) =
{
  L = N -> exp(sqrt(log(N)*log(log(N))));
  B = N -> ceil(L(N)^(1/sqrt(2)));
  f = T -> T^2 - N;
  tValues = vector(T, floor(sqrt(N)) + 1, floot(sqrt(N) + 1 + ceil(L(N)^(1/sqrt(2)))));
  print (tValues);
}
