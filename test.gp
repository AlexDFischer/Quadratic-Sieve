isBsmooth(B, N) =
{
  local(factors);
  factors = factor(N)[,1]~;
  return(factors[#factors] <= B);
}

B=(N)->ceil(L(N)^(1/sqrt(2)))

log2 = x->log(x)/log(2)

L(N)=exp(sqrt(log(N)*log(log(N))))
