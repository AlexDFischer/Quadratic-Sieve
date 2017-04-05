#include "sieve.h"

int main(int argc, char **argv)
{
  mpz_t N;
  mpz_init(N);
  float numRelationsPower = sqrt(2), smoothnessPower = 1 / sqrt(2);
  switch (argc)
  {
    case 1:
      mpz_set_ui(N, 310069);
      break;
    case 2:
      gmp_sscanf(argv[1], "%Zd", N);
      break;
    case 4:
      gmp_sscanf(argv[1], "%Zd", N);
      sscanf(argv[2], "%f", &numRelationsPower);
      sscanf(argv[3], "%f", &smoothnessPower);
  }
  factor(N, numRelationsPower, smoothnessPower);
  mpz_clear(N);
}
