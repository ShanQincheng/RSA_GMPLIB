#include "random_number.h"


const unsigned long int COMPLEX = 1024; 

void RandomMpz_t(mpz_t *randomNum, const mpz_t limitNum, unsigned long int seedParameter)
{
  unsigned long int seed = 0;
  srand((unsigned)time(NULL));
  seed = rand() * seedParameter;
  gmp_randstate_t state;
  gmp_randinit_mt(state);  // random state initialization
  gmp_randseed_ui(state, seed); // random state seeding

  if( mpz_cmp_ui(limitNum, 0) == 0)  // Generate a random number from 2^(n-1) to 2^n-1  
  {
    mp_bitcnt_t n = COMPLEX;
    mpz_rrandomb(*randomNum, state, n);  // Generate a random integer with long strings of zeros and ones in the binary representation. 2^(n-1) ~ 2^n-1
  }
  else if( mpz_cmp_ui(limitNum, 0) > 0 )  //  Generate a random integer in the range 0 to (limitNum-1) 
  {
    mpz_urandomm(*randomNum, state, limitNum);  // Generate a uniform random integer in the range 0 to limitNum - 1, inclusive   
  }
  else assert(0);
  
  gmp_randclear(state);
  return;
}

void RandomInteger(int *randomNum)  // random an integer in the range  0 to 19 
{
  unsigned long int seed = 0;
  srand((unsigned)time(NULL));
  *randomNum =  rand() % 20;

  return;
}
 
