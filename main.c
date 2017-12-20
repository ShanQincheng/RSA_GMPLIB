#include <stdio.h>
#include "rsa_algorithm.h"
#include "select_parameter.h"
#include "random_number.h"
#include "miller_rabin_primality_test.h"


int main(void)
{
  mpz_t p; // random COMPLEX bit prime p
  mpz_t q; // random COMPLEX bit prime q
  mpz_t n; // by the equation n = pq 
  mpz_t _n; // by equation _n = (p - 1)(q - 1)
  //unsigned long int d;
  mpz_t d;
  unsigned long int e; // a small odd integer e
  mpz_t limitNum; // set as random type or random number max limit
  mp_bitcnt_t bit_index = 0;  // the least significant bit is number 0 for function int mpz_tstbit()
  struct PublicKeyPair
  { 
    unsigned long int e;
    mpz_t n;
  };
  struct PublicKeyPair public_key_pair;
  
  struct PrivateKeyPair
  {
    //unsigned long int d;
    mpz_t d;
    mpz_t n;
  };
  struct PrivateKeyPair private_key_pair;
 
  mpz_inits(limitNum, d, NULL);
 // mpz_init2(d, COMPLEX);
  mpz_init2(p, COMPLEX);
  mpz_init2(q, COMPLEX);
  bit_index = 2 * COMPLEX;
  mpz_init2(n, bit_index);
  mpz_init2(_n, bit_index);
  bit_index = 0;
  /*
     Random COMPLEX bit p
  */
 // mpz_set_ui(limitNum, 0);
 // RandomMpz_t(&p, limitNum);
  Select_p(&p);
  if(mpz_tstbit(p, bit_index))
  {
    while( MILLER_RABIN(p) ) mpz_add_ui(p, p, 2);
  }else {
    mpz_add_ui(p, p, 1);
    while( MILLER_RABIN(p) ) mpz_add_ui(p, p, 2);
  }
  /*
    Random COMPLEX bit q
  */
  //RandomMpz_t(&q, limitNum);
  Select_q(&q);
  if(mpz_tstbit(q, bit_index))
  {
    while( MILLER_RABIN(q) ) mpz_add_ui(q, q, 2);
  }else {
    mpz_add_ui(q, q, 1);
    while( MILLER_RABIN(q) ) mpz_add_ui(q, q, 2);
  }
  gmp_printf("p == %Zd\n", p);
  gmp_printf("q == %Zd\n", q);
  Select_n(p, q, &n);  // n = pq
  gmp_printf("n (p * q) == %Zd\n", n);
  Select_Euler_n(p, q, &_n); // Euler n = (p-1)(q-1)
  gmp_printf("_n (p-1)*(q-1) == %Zd\n", _n);
  Select_e(&e, _n); //  gcd(_n, e)
  mpz_init(public_key_pair.n);
  mpz_set(public_key_pair.n, n); // P = (e. n) as RSA public key
  public_key_pair.e = e;

  Select_d(e, _n, &d);
  //Select_d(_n, e, &d);
  //unsigned long int ddd = 0; 
  //Select_d_simple(e, _n, &ddd);
  //mpz_set_ui(d, ddd);
   
 // private_key_pair.d = d;
  mpz_init(private_key_pair.d);
  mpz_set(private_key_pair.d, d);
  mpz_init(private_key_pair.n);
  mpz_set(private_key_pair.n, n);  
  gmp_printf("public_key_n is %Zd\n", public_key_pair.n);
  printf("public_key_e is %lu\n", public_key_pair.e);
  gmp_printf("private_key_n is %Zd\n", public_key_pair.n);
  //printf("private_key_d is %lu\n", private_key_pair.d);
  gmp_printf("private_key_d is %Zd\n", private_key_pair.d);

  mpz_t C; // Ciphertext
  mpz_t M, N; // Cleartext
  mpz_inits(C, M, N, NULL);

  //M = COMPLEX;
  mpz_set_ui(M, 5201314);
  RSAEncryption(M, public_key_pair.e, public_key_pair.n, &C);
  gmp_printf("The Ciphertext is %Zd\n", C);
  RSADecrypt(C, private_key_pair.d, private_key_pair.n, &N);
  gmp_printf("The Clear text is %Zd\n", N); 

  mpz_clear(p);
  return 0; 
}


