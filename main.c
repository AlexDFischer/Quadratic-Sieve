#include "sieve.h"
#include <string.h>

int main(int argc, char **argv)
{
  if (argc > 1 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
  {
    printf("Tries to factor a number, either the first command line argument or 310069 if none is given. Optional second and third arguments: if the second and third arguments are a and b respectively, it will try L(N)^a random numbers in searching for relations and a smoothness bound of L(N)^b. Default values are sqrt(2) and 1/sqrt(2), respectively.\n");
    exit(0);
  }
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
