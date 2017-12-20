#include "miller_rabin_primality_test.h"


bool WITNESS(const mpz_t a, const mpz_t n)
{
  mpz_t u, n_1, x[2];  //   
  unsigned long int t = 0, base = 2; 
   
  mpz_inits(u, n_1, x[0], x[1], NULL);
  /*
    let n - 1 = 2^t * u, where t >= 1 and u is odd
  */
  mpz_sub_ui(n_1, n, 1);
  do {
    mpz_t d;  // divisor
    mpz_init(d);
 
    t++; 
    mpz_ui_pow_ui(d, base, t);
    mpz_divexact(u, n_1, d);
  }while(!mpz_tstbit(u, 0)); 
  /*
    X0 <- MODULAR-EXPONENTIATION(a, u, n)
  */ 
  mpz_powm(x[0], a, u, n);  // X0 = a^u mod n
  for(int i = 0; i < t; i++)
  {
    mpz_powm_ui(x[1], x[0], base, n);  // do Xi <- Xi-1^2 mod n
    if(!mpz_cmp_ui(x[1], 1) && mpz_cmp_ui(x[0], 1) != 0 && mpz_cmp_ui(x[0], -1) != 0 && mpz_cmp(x[0], n_1) != 0)
      return true;
    else
      mpz_set(x[0], x[1]);
  } 
  if(mpz_cmp_ui(x[1], 1) != 0)
    return true;
  else return false; 
    
}

bool MILLER_RABIN(const mpz_t n)  // prime test
{
   int s = 0;
   RandomInteger(&s);
   mpz_t a;
   mpz_init(a);

   for(int i = 0; i < s; i++)
   {
     //RandomMpz_t(&a, n);
     //RandomInteger(&a);
     mpz_set_ui(a, i+2); 
     if( WITNESS(a, n) )  return true;
   }
   return false; 
}


