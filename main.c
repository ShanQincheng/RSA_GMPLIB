#include <stdio.h>
#include "rsa_algorithm.h"
#include "select_parameter.h"
#include "random_number.h"
#include "miller_rabin_primality_test.h"

int main(void)
{
  mpz_t p; // random COMPLEX bit prime p
  mpz_t q; // random COMPLEX bit prime q, q != p
  mpz_t n; // by the equation n = pq 
  mpz_t _n; // by equation _n = (p - 1)(q - 1)
  mpz_t d;  // d is the multiplicative inverse of e, modulo _n; we can use Extend_EUCLID algorithm to solve this problem
  unsigned long int e; // a small odd integer e
  mpz_t limitNum; // set as random type or random number max limit
  mp_bitcnt_t bit_index = 0;  // the least significant bit is number 0 for function int mpz_tstbit()
  struct PublicKeyPair  // rsa public key pair struct
  { 
    unsigned long int e;
    mpz_t n;
  };
  struct PublicKeyPair public_key_pair;
  
  struct PrivateKeyPair  // rsa private key pair struct
  {
    mpz_t d;
    mpz_t n;
  };
  struct PrivateKeyPair private_key_pair;
 
  mpz_inits(limitNum, d, NULL);  // use mpz_inits or mpz_init2 function to initialize mpz_t type parameters before use them
  mpz_init2(p, COMPLEX);
  mpz_init2(q, COMPLEX);
  bit_index = 2 * COMPLEX; // x bits number multiple y bits number , the result will be a (x+y) bits number
  mpz_init2(n, bit_index); 
  mpz_init2(_n, bit_index); // the same reason to n
  bit_index = 0; // recovery bit_index to 0
  
  Select_p(&p);  // use random mpz_z type number algorithm to choose p
  Select_q(&q);  // the sam reason to p
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
  gmp_printf("public_key_n is %Zd\n", public_key_pair.n);
  printf("public_key_e is %lu\n", public_key_pair.e);

  Select_d(e, _n, &d);  // use EXTEND EUCLID algorithm choose the unique d
  mpz_init(private_key_pair.d);
  mpz_set(private_key_pair.d, d);
  mpz_init(private_key_pair.n);
  mpz_set(private_key_pair.n, n);  
  gmp_printf("private_key_n is %Zd\n", private_key_pair.n);
  gmp_printf("private_key_d is %Zd\n", private_key_pair.d);

  mpz_t C; // Ciphertext
  mpz_t M, N; // Cleartext
  mpz_inits(C, M, N, NULL);

  mpz_set_ui(M, 5201314);
  RSAEncryption(M, public_key_pair.e, public_key_pair.n, &C); // use rsa algorithm encrypt cleartext M, and then store ciphertext in variable C
  gmp_printf("The Ciphertext is %Zd\n", C);
  RSADecrypt(C, private_key_pair.d, private_key_pair.n, &N);  // use rsa algorithm decrypt ciphertext C, and store cleartext in variable N 
  gmp_printf("The Clear text is %Zd\n", N); 

  mpz_clear(p);
  return 0; 
}


