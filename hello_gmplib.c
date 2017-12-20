#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <gmp.h>

const unsigned long int COMPLEX = 1024; 
const unsigned long int COMPLEX_e = 65536;
bool WITNESS(const mpz_t, const mpz_t);

void mpzDiv(const mpz_t a, const mpz_t b, mpz_t* result)
{
  mpz_cdiv_q(*result, a, b);
  return;
}

void mpzMul(const mpz_t a, const mpz_t b, mpz_t* result)
{
  mpz_mul(*result, a, b);
  return;
}

void mpzSub(const mpz_t a, const mpz_t b, mpz_t* result)
{
  mpz_sub(*result, a, b);
  return;
}

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

void xgcd(mpz_t a, mpz_t b, mpz_t *d, mpz_t* x, mpz_t* y)
{
  mpz_t aa[2], bb[2], q, r;
  mpz_t temp_mpz;
  mpz_inits(aa[0], aa[1], bb[0], bb[1], q, r, NULL);
  mpz_set_ui(aa[0], 1);
  mpz_set_ui(aa[1], 0);
  mpz_set_ui(bb[0], 0);
  mpz_set_ui(bb[1], 1);
  mpz_init(temp_mpz);
  while(1) {
    mpz_cdiv_q(q, a, b); //  q = a / b
    //mpz_set(temp_mpz, a);  
    mpz_powm_ui(a, a, 1, b);  // a = a % b 
    mpz_mul(temp_mpz, q, aa[1]); //  
    mpz_sub(aa[0], aa[0], temp_mpz);  // aa[0] = aa[0] - q*aa[1]
    mpz_mul(temp_mpz, q, bb[1]);
    mpz_sub(bb[0], bb[0], temp_mpz);  // bb[0] = bb[0] - q*bb[1]
    if(mpz_cmp_ui(a, 0) == 0)  // if a==0
    {
      mpz_set(*d, b);  // d = b;
      mpz_set(*x, aa[1]); // x = aa[1]
      mpz_set(*y, bb[1]);  // y = bb[1]
      return;
    } 
    mpz_cdiv_q(q, b, a);  // q = b / a 
    mpz_powm_ui(b, b, 1, a);  // b = b % a
    mpz_mul(temp_mpz, q, aa[0]); //  
    mpz_sub(aa[1], aa[1], temp_mpz);  // aa[1] = aa[1] - q*aa[0]
    mpz_mul(temp_mpz, q, bb[0]);
    mpz_sub(bb[1], bb[1], temp_mpz);  // bb[1] = bb[1] - q*bb[0]
    if(mpz_cmp_ui(b, 0) == 0)
    {
      mpz_set(*d, a);
      mpz_set(*x, aa[0]);
      mpz_set(*y, bb[0]);
      return;
    }	
  }
}

void EXTENDED_ECUCLID(mpz_t a, mpz_t b, mpz_t *d, mpz_t *x, mpz_t *y)
{
  mpz_t q, r, x1, x2, y1, y2;
  mpz_inits(q, r, x1, x2, y1, y2, NULL);

  mpz_t result_mul;
  mpz_init(result_mul);
  if( mpz_cmp_ui(b, 0) == 0)
  {
    mpz_set(*d, a);
    mpz_set_ui(*x, 1);
    mpz_set_ui(*y, 0);
    return; 
  }
  mpz_set_ui(x2, 1); mpz_set_ui(x1, 0); mpz_set_ui(y2, 0); mpz_set_ui(y1, 1);
  while( mpz_cmp_ui(b, 0) > 0)
  {
    mpz_cdiv_q(q, a, b);  // q = a / b
    mpz_mul(result_mul, q, b);  
    mpz_sub(r, a, result_mul); // r = a - q * b
    mpz_mul(result_mul, q, x1);
    mpz_sub(*x, x2, result_mul); // *x = x2 - q * x1
    mpz_mul(result_mul, q, y1);
    mpz_sub(*y, y2, result_mul);  // *y = y2 - q * y1
    mpz_set(a, b); mpz_set(b, r);  // a = b, b = r
    mpz_set(x2, x1);  // x2 = x1
    mpz_set(x1, *x);  // x1 = *x
    mpz_set(y2, y1);  // y2 = y1
    mpz_set(y1, *y);  // y1 = *y
  }
  mpz_set(*d, a);  // *d = a
  mpz_set(*x, x2);  // *x = x2
  mpz_set(*y, y2);  // *y = y2
  return;
}

//void Select_d(const unsigned long int e, const mpz_t _n, unsigned long int *d)
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

void Select_d_simple(const unsigned long int e, const mpz_t _n, unsigned long int *d)
{
  mpz_t mpz_t_ed;
  mpz_init(mpz_t_ed);
  *d = 0;
  do {
    (*d)++;
    mpz_set_ui(mpz_t_ed, e*(*d));
    mpz_powm_ui(mpz_t_ed, mpz_t_ed, 1, _n);
  }while(mpz_cmp_si(mpz_t_ed, 1) != 0);
 
} 

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


int main(void)
{
  mpz_t p; // random COMPLEX bit prime p
  mpz_t q; // random COMPLEX bit prime q
  mpz_t n; // by the equation n = pq 
  mpz_t _n; // by equation _n = (p - 1)(q - 1)
  //unsigned long int d;
  static  mpz_t d;
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
  mpz_set_ui(M, 500);
  RSAEncryption(M, public_key_pair.e, public_key_pair.n, &C);
  gmp_printf("The Ciphertext is %Zd\n", C);
  RSADecrypt(C, private_key_pair.d, private_key_pair.n, &N);
  gmp_printf("The Clear text is %Zd\n", N); 

//  mpz_t result_test, a_test, b_test;
 // mpz_inits(result_test, a_test, b_test, NULL);
 // mpz_set_ui(a_test, 10);
  //mpz_set_ui(b_test, 20);  
  //mpzMul(a_test, b_test, &result_test);
  //gmp_printf("%Zd\n", result_test);
 //mpz_inits(n, e, NULL);
 // unsigned long int op[] = {561, 617, 653, 879, 917, 973};
 // mpz_set_ui(e, 2);
//  for(int i = 0; i < 6; i++)
//  {
//    mpz_set_ui(n, op[i]);
//    if( MILLER_RABIN(n) )
//        printf("%ld is a COMPOSITE\n", op[i]);
//    else printf("%ld is a PRIME\n", op[i]);
//  }  

  mpz_clear(p);
  return 0; 
}


