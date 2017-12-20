#include "select_parameter.h"


void Select_p(mpz_t* p)
{
  unsigned long int seedParameter = 3.1415926553660;
  mpz_t limitNum;
  mpz_init(limitNum);
  mpz_set_ui(limitNum, 0);
  RandomMpz_t(p, limitNum, seedParameter);
}

void Select_q(mpz_t* q)
{
  unsigned long int seedParameter = 15659269653;
  mpz_t limitNum;
  mpz_init(limitNum);
  mpz_set_ui(limitNum, 0);
  RandomMpz_t(q, limitNum, seedParameter);
}


void Select_n(const mpz_t p, const mpz_t q, mpz_t *n)
{
  mpz_mul(*n, p, q);  // set n to p times q  
}

void Select_Euler_n(const mpz_t p, const mpz_t q, mpz_t *_n)
{
  mpz_t _p, _q;
  mpz_inits(_p, _q, NULL);

  mpz_sub_ui(_p, p, 1);
  mpz_sub_ui(_q, q, 1);
  mpz_mul(*_n, _p, _q);
}

void Select_e(unsigned long int *e, const mpz_t _n)
{
  const unsigned long int COMPLEX_e = 65536;
  mpz_t rop;
 
  *e = COMPLEX_e;  // initial  e, test start from COMPLEX_e 
  mpz_init(rop);
  mpz_set_ui(rop, 0);
  while( mpz_cmp_ui(rop, 1 ) != 0 )
  {
    (*e) += 1;
    mpz_gcd_ui(rop, _n, *e);  
  } 
}

void Select_d(const unsigned long int e, mpz_t _n, mpz_t *d)
{
  mpz_t mpz_t_e, _d, y;
  mpz_inits(mpz_t_e, _d, y, NULL);
  mpz_set_ui(mpz_t_e, e);

  //EXTENDED_ECUCLID(mpz_t_e, _n, &_d, d, &y);
//  EXTENDED_ECUCLID(_n, mpz_t_e, &_d, d, &y);
//  mpz_set_ui(*d, EXTENDED_ECUCLID(mpz_t_e, _n)[1]); 
//  xgcd(mpz_t_e, _n, &_d, d, &y);
  //xgcd(_n, mpz_t_e, &_d, d, &y);
  //EXTENDED_ECUCLID(_n, mpz_t_e, &_d, d, &y);
  mpz_gcdext(_d, *d, y, mpz_t_e, _n); 
}


