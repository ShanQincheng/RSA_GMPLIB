#ifndef _SELECT_PARAMETER_HEAD_
#define _SELECT_PARAMETER_HEAD_

#include <gmp.h>
#include "random_number.h" 


void Select_p(mpz_t*);
void Select_q(mpz_t*);
void Select_n(const mpz_t, const mpz_t, mpz_t*);
void Select_Euler_n(const mpz_t, const mpz_t, mpz_t*);
void Select_e(unsigned long int*, const mpz_t);
void Select_d(const unsigned long int , mpz_t, mpz_t*);

#endif

