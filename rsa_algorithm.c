#include "rsa_algorithm.h"


void RSAEncryption(const mpz_t M, const unsigned long int e, const mpz_t n, mpz_t* C)
{
 //mpz_t M_mpz;
 //mpz_init(M_mpz);
 //mpz_set_ui(M_mpz, M); 
 mpz_powm_ui(*C, M, e, n); 
}

void RSADecrypt(const mpz_t C, const mpz_t d, const mpz_t n, mpz_t* M)
{
  //mpz_t M_mpz;
  //mpz_init(M_mpz);

  mpz_powm(*M, C, d, n);
  //*M = mpz_get_ui(M_mpz);
} 


