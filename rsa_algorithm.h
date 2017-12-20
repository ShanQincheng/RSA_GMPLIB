#ifndef _RSA_ALGORITHM_HEAD_
#define _RSA_ALGORITHM_HEAD_

#include <gmp.h>


void RSAEncryption(const mpz_t, const unsigned long int, const mpz_t, mpz_t*);

void RSADecrypt(const mpz_t, const mpz_t, const mpz_t, mpz_t*);

#endif
