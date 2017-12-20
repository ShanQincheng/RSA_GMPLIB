#ifndef _MILLER_RABIN_PRIMALITY_TEST_HEAD_
#define _MILLER_RABIN_PRIMALITY_TEST_HEAD_

#include <stdbool.h>
#include <gmp.h>
#include "random_number.h"

bool WITNESS(const mpz_t, const mpz_t);  // MILLER RABIN primality test core code
bool MILLER_RABIN(const mpz_t);  // use random number times to primality test the tested odd 


#endif



