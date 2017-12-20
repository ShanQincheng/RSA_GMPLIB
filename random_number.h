#ifndef _RANDOM_NUMBER_HEAD_
#define _RANDOM_NUMBER_HEAD_


#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <gmp.h>

extern const unsigned long int COMPLEX;
void RandomMpz_t(mpz_t *randomNum, const mpz_t limitNum, unsigned long int seedParameter);
void RandomInteger(int *randomNum);  // random an integer in the range  0 to 19 


#endif
