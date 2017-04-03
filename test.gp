{
isBsmooth(B, N) =
  local(factors);
  factors = factor(N)[,1]~;
  return(factors[#factors] <= B);
}

B=(N)->ceil(L(N)^(1/sqrt(2)))

log2 = x->log(x)/log(2)

L(N)=exp(sqrt(log(N)*log(log(N))))

{
printRelations(N) =
  for (T=floor(sqrt(N))+1, floor(sqrt(N))+1+ceil(L(N)^sqrt(2)),
    if (isBsmooth(B(N), T^2-N), print("T=",T));
  );
}
