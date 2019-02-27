/*
    Copyright (c) 2011 Luca De Feo.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
	Modifed by Dieter Fishbein
*/
  
#include <pthread.h>
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <openssl/rand.h>
#include <sys/time.h>
#define DEBUG
#define GF_TMP_REGS 9
#define MAX(a,b) (((a)>(b))? (a):(b))


/****************** TYPES *****************/

struct GF_params;
typedef struct GF_params GF_params;
int count=1;
int countAdd =0;
int dieter=0;
int version=2;
int count1=0;
int count2=0;
int count3=0;
int count2iso=0;
int count2mul=0;
int count3iso=0;
int count3mul=0;
int errorCount = 0;
int countMul = 0;
int countBarrett = 0;
int countIF = 0;
int size=0;
int zero=0;
int addCount=0;
int barrettCount;

//Barreet precomputations
mpz_t mB;
mpz_t kF;
int * kB;

int negFlag;


mpz_t prime;
mpz_t tmp;
int MPZ_MEMORY;

// Elements of GF(p^2)
typedef struct {
  GF_params* parent;
  mpz_t a, b;
} GF;


// The field GF(p^2)
// basically its characteristic and some work registers
struct GF_params {
    mpz_t p, tmp1, tmp2, tmp3, tmp4, tmp5;
  GF GFtmp[GF_TMP_REGS];
  gmp_randstate_t state;
  int initialized;
};


//threading stuff
typedef struct {
    GF *a, *x, *y;
    long int *u;
} arith;

void *add_t(void *ref){
    arith *arithInfo;
    arithInfo = (arith *) ref;
    add_GF(arithInfo->a, *(arithInfo->x), *(arithInfo->y));
}

/****some global variables for threading*****/

typedef struct queue_point queue_point;

typedef struct mt_info {
    
    queue_point *qp2, *qp3;
    
}mt_info;

mt_info *t;
pthread_t threads[2];



/******** IMPLEMENTATION OF GF(p^2) ********/

// Memory management
void init_GF(GF* x, GF_params* parent) {
  mpz_init2(x->a, MPZ_MEMORY*64);
  mpz_init2(x->b, MPZ_MEMORY*64);
    
 //    mpz_init(x->a);
  //   mpz_init(x->b);
    
    
  x->parent = parent;
}

void clear_GF(GF *x) {
  mpz_clear(x->a);
  mpz_clear(x->b);
}

// Initialization of GF(p,2)
int setup_GF(GF_params* field, const char* characteristic) {
  if (!characteristic) {
    // If p is NULL, use 2^387 * 3^242 - 1 
    // as default value
    mpz_init_set_ui(field->p, 3);
    mpz_pow_ui(field->p, field->p, 242);
    mpz_mul_2exp(field->p, field->p, 387);
    mpz_sub_ui(field->p, field->p, 1);
  } else {
	mpz_init( field->p );
    mpz_set_str(field->p, characteristic, 0);
  }
  // Check that the Legendre symbol of -1 is -1 (p = 3 mod 4)
  if (mpz_fdiv_ui(field->p, 4) != 3) {
    mpz_clear(field->p);
    return 0;
  }

  gmp_randinit_default(field->state);
    mpz_init2(field->tmp1, MPZ_MEMORY*64); mpz_init2(field->tmp2, MPZ_MEMORY*64); mpz_init2(field->tmp3, MPZ_MEMORY*64); mpz_init2(field->tmp4, MPZ_MEMORY*64); mpz_init2(field->tmp5, MPZ_MEMORY*64);

  int i;
  for (i = 0 ; i < GF_TMP_REGS ; i++)
    init_GF(&field->GFtmp[i], field);
  field->initialized = 1;
  
  return 1;
}

void free_GF(GF_params* field) {
  if (field->initialized) {
    int i;
    for (i = 0 ; i < GF_TMP_REGS ; i++)
      clear_GF(&field->GFtmp[i]);
    mpz_clear(field->p); mpz_clear(field->tmp1); mpz_clear(field->tmp2); mpz_clear(field->tmp3); //mpz_clear(field->tmp4); mpz_clear(field->tmp5);
    field->initialized = 0;
  }
}

// IO of elements
// outputs are strings in base 16
// (to save some bytes)
void set_GF(GF* x, const char* a, const char* b) {
  mpz_set_str(x->a, a, 0);
  mpz_set_str(x->b, b, 0);
}

void set_GFc(GF* x, const char* p) {
	mpz_set_str(x->parent->p, p, 0);	
}


void get_GF(char *a, char *b, const GF x) {
  gmp_sprintf(a, "%#Zx", x.a);
  gmp_sprintf(b, "%#Zx", x.b);
}


// Arithmetic modulo X^2 + 1

/*
  There seems to be a bug in GMP 4.2.1 that makes mpz_mod give
  unpredictable results when the mpz_t holding the result is the same
  as one of the operands.
*/

void copy_GF(GF* res, const GF x){
  mpz_set(res->a, x.a);
  mpz_set(res->b, x.b);
  res->parent = x.parent;
}

void add_GF(GF *res, const GF x, const GF y) {
  
    
    if( (x.a)->_mp_size==0 && (y.a)->_mp_size==0){
        (res->a)->_mp_size = 0;
    }else{    
    
    mpz_add(res->a,x.a,y.a);
    
    //  if( res->_mp_alloc < 13 )
    //      _mpz_realloc(res,13);
    
    // printf("\nIN FUNC: res size: %d\n", res->_mp_size);
    //  gmp_printf("\nIN FUNC: res : %Zd\n", res);
        
    
    if( mpz_cmp(res->a,prime) >= 0){
        
        asm   (  
               "movq	$0xffffffffffffffff, %%r8;"
               "subq	%%r8, %0;"
               
             //  "movq	$0xffffffffffffffff, %%r8;"     //problem line
               "sbbq	%%r8, %1;"   
               
              // "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %2;"
               
              // "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %3;"
               
              // "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %4;"
               
               //"movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %5;"
               
               "movq   $0xf007669a5ce89647, %%r8;"
               "sbbq	%%r8, %6;"
               
               "movq	$0xade00d91484504f9, %%r8;"
               "sbbq	%%r8, %7;"
               
               "movq	$0x0979d570c24486e3, %%r8;"
               "sbbq	%%r8, %8;"
               
               "movq	$0x8bbae3679a4c7025, %%r8;"
               "sbbq	%%r8, %9;"
               
               "movq	$0xa06a805a9f6808b4, %%r8;"
               "sbbq	%%r8, %10;"
               
               "movq	$0xe69ebefa87fabdfa, %%r8;"
               "sbbq	%%r8, %11;"
               
               "movq	$0x5, %%r8;"
               "sbbq   %%r8, %12;"
               
               
               : "+X" (((res->a)->_mp_d)[0]), "+X" (((res->a)->_mp_d)[1]), "+X" (((res->a)->_mp_d)[2]), "+X" (((res->a)->_mp_d)[3]), "+X" (((res->a)->_mp_d)[4]), "+X" (((res->a)->_mp_d)[5]), "+X" (((res->a)->_mp_d)[6]), "+X" (((res->a)->_mp_d)[7]), "+X" (((res->a)->_mp_d)[8]), "+X" (((res->a)->_mp_d)[9]), "+X" (((res->a)->_mp_d)[10]), "+X" (((res->a)->_mp_d)[11]), "+X" (((res->a)->_mp_d)[12])
               :
               :"%r8"
               );
        
        //  printf("IN FUNC RESULT:\n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
    }
    
    if( ((res->a)->_mp_d)[12]==0 )
        (res->a)->_mp_size = 12;
    else
        (res->a)->_mp_size = 13;
    
    }
    
    
    
    
    if( (x.b)->_mp_size==0 && (y.b)->_mp_size==0){
        (res->b)->_mp_size = 0;
    }else  if( (x.b)->_mp_size==1 && (y.b)->_mp_size==0){
        
        mpz_set_ui(res->b,1);
    } else{
    
    
    mpz_add(res->b, x.b, y.b);
    
    //  if( res->_mp_alloc < 13 )
    //      _mpz_realloc(res,13);
    
    // printf("\nIN FUNC: res size: %d\n", res->_mp_size);
    //  gmp_printf("\nIN FUNC: res : %Zd\n", res);
    
    if( mpz_cmp(res->b,prime) >= 0){
        
        asm   (  
               "movq	$0xffffffffffffffff, %%r8;"
               "subq	%%r8, %0;"
               
           //    "movq	$0xffffffffffffffff, %%r8;"     //problem line
               "sbbq	%%r8, %1;"   
               
            //   "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %2;"
               
             //  "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %3;"
               
             //  "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %4;"
               
             //  "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %5;"
               
               "movq   $0xf007669a5ce89647, %%r8;"
               "sbbq	%%r8, %6;"
               
               "movq	$0xade00d91484504f9, %%r8;"
               "sbbq	%%r8, %7;"
               
               "movq	$0x0979d570c24486e3, %%r8;"
               "sbbq	%%r8, %8;"
               
               "movq	$0x8bbae3679a4c7025, %%r8;"
               "sbbq	%%r8, %9;"
               
               "movq	$0xa06a805a9f6808b4, %%r8;"
               "sbbq	%%r8, %10;"
               
               "movq	$0xe69ebefa87fabdfa, %%r8;"
               "sbbq	%%r8, %11;"
               
               "movq	$0x5, %%r8;"
               "sbbq   %%r8, %12;"
               
               
               : "+X" (((res->b)->_mp_d)[0]), "+X" (((res->b)->_mp_d)[1]), "+X" (((res->b)->_mp_d)[2]), "+X" (((res->b)->_mp_d)[3]), "+X" (((res->b)->_mp_d)[4]), "+X" (((res->b)->_mp_d)[5]), "+X" (((res->b)->_mp_d)[6]), "+X" (((res->b)->_mp_d)[7]), "+X" (((res->b)->_mp_d)[8]), "+X" (((res->b)->_mp_d)[9]), "+X" (((res->b)->_mp_d)[10]), "+X" (((res->b)->_mp_d)[11]), "+X" (((res->b)->_mp_d)[12])
               :
               :"%r8"
               );
        
        //    printf("IN FUNC RESULT:\n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
    }
    
    if( ((res->b)->_mp_d)[12]==0 )
        (res->b)->_mp_size = 12;
    else
        (res->b)->_mp_size = 13;
    
    }
    
 /*   mpz_add(x.parent->tmp1, x.a, y.a);
    mpz_mod(res->a, x.parent->tmp1, x.parent->p);
    mpz_add(x.parent->tmp1, x.b, y.b);
    mpz_mod(res->b, x.parent->tmp1, x.parent->p);
*/
    
    //  mpz_addm_x86(res->a, x.a, y.a);
    //  mpz_addm_x86_2(res->b, x.b, y.b);
    
  /*  int flag=0;
    if( (x.a)->_mp_size==0 ^ (y.a)->_mp_size==0 ){
        flag=1;
        printf("\n********** %d ***********\n", addCount);
        gmp_printf("\nx.a: %Zd\n", x.a);
        gmp_printf("\ny.a: %Zd\n", y.a);
        printf("\nx.a size: %d\n", x.a->_mp_size);
        printf("\ny.a size: %d\n", y.a->_mp_size);  
    } 
*/
 //   mpz_add(x.parent->tmp1, x.a, y.a);
 //   mpz_mod(x.parent->tmp2, x.parent->tmp1, x.parent->p);
  
  
 //   mpz_addm_x86(res->a, x.a, y.a);
    
 /*   if( mpz_cmp(res->a,x.parent->tmp2) !=0 )
        printf("\nFAIL\n");
  */ 
  /*  if(flag){
        gmp_printf("\nres: %Zd\n", res->a);
        printf("\nres size: %d\n", res->a->_mp_size);
        printf("\n************************\n", addCount);
    }
    */
   // mpz_addm_x86_2(res->b, x.b, y.b);
    
    res->parent = x.parent;
    
  //  addCount++;
}    

void add_GF_r(GF *res, const GF x, const GF y) {
     mpz_addm_x86(res->a, x.a, y.a);
      mpz_addm_x86(res->b, x.b, y.b);
    
    
    res->parent = x.parent;
    
   // addCount++;
   
}

void add_GF_t(GF *res, const GF x, const GF y, mpz_t *tmp) {
    
/*    mpz_t tmp1;
    mpz_init(tmp1); */ 
    
    mpz_add(tmp[0], x.a, y.a);
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_add(tmp[0], x.b, y.b);
    mpz_mod(res->b, tmp[0], x.parent->p);
    res->parent = x.parent;
    
  //  mpz_clear(tmp1);
}

void add_GF_ui(GF *res, const GF x, unsigned long int u) {
  mpz_add_ui(x.parent->tmp1, x.b, u);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  mpz_set(res->a, x.a);
  res->parent = x.parent;
}

void sub_GF(GF *res, const GF x, const GF y) {
  mpz_sub(x.parent->tmp1, x.a, y.a);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_sub(x.parent->tmp1, x.b, y.b);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

//specify mpz tmp variables in tmp
void sub_GF_t(GF *res, const GF x, const GF y, mpz_t *tmp) {
    
 //   mpz_t tmp1;
 //   mpz_init(tmp1);    
    
    mpz_sub(tmp[0], x.a, y.a);
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_sub(tmp[0], x.b, y.b);
    mpz_mod(res->b, tmp[0], x.parent->p);
    res->parent = x.parent;
    
  //  mpz_clear(tmp1);
    
}

void sub_GF_ui(GF *res, const GF x, unsigned long int u) {
  mpz_sub_ui(x.parent->tmp1, x.b, u);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  mpz_set(res->a, x.a);
  res->parent = x.parent;
}
//specify tmps
void sub_GF_ui_t(GF *res, const GF x, unsigned long int u, mpz_t *tmp) {
    
   // mpz_t tmp1;
   // mpz_init(tmp1);
    
    mpz_sub_ui(tmp[0], x.b, u);
    mpz_mod(res->b, tmp[0], x.parent->p);
    mpz_set(res->a, x.a);
    res->parent = x.parent;
    
 //   mpz_clear(tmp1);
}

void neg_GF(GF *res, const GF x) {
  if (mpz_sgn(x.a) == 0)
    mpz_set(res->a, x.a);
  else
     mpz_sub(res->a, x.parent->p, x.a); 
                                
  if (mpz_sgn(x.b) == 0)
    mpz_set(res->b, x.b);
  else
    mpz_sub(res->b, x.parent->p, x.b);
  res->parent = x.parent;
}

void scalar_GF(GF *res, const GF x, mpz_t s) {
  mpz_mul(x.parent->tmp1, x.a, s);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mul(x.parent->tmp1, x.b, s);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void scalar_GF_si(GF *res, const GF x, long int s) {
  mpz_mul_si(x.parent->tmp1, x.a, s);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mul_si(x.parent->tmp1, x.b, s);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void scalar_GF_si_t(GF *res, const GF x, long int s, mpz_t *tmp) {
    
  //  mpz_t tmp1, tmp2, tmp3;
  //  mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);    
    
    mpz_mul_si(tmp[0], x.a, s);
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_mul_si(tmp[0], x.b, s);
    mpz_mod(res->b, tmp[0], x.parent->p);
    res->parent = x.parent;
    
 //   mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
}


void mul_GF(GF *res, const GF x, const GF y) {
    
/*    printf("\n****** %d ********", countMul) ;
    printf("\nInputs to mul_GF ") ; 
    gmp_printf("\n x.a: %Zd", x.a);
    gmp_printf("\n x.b: %Zd", x.b);
    gmp_printf("\n y.a: %Zd", y.a);  
    gmp_printf("\n y.b: %Zd\n", y.b);  
    printf("\nInputs to mpz_mul ", countMul) ;  
  */ 
   // mpz_add(y.parent->tmp1, x.a, x.b);
   // mpz_sub(y.parent->tmp2, y.b, y.a);
    
    
  //  mpz_mul_x86_1(y.parent->tmp3, y.parent->tmp1, y.parent->tmp2);
//    mpz_mul(y.parent->tmp3, y.parent->tmp1, y.parent->tmp2);
    
  //  int flag=0;    
 // if( ((x.a)->_mp_size!=0 && (x.a)->_mp_size!=12 && (x.a)->_mp_size!=13) || ( (y.b)->_mp_size!=0  && (y.b)->_mp_size != 1 && (y.b)->_mp_size != 12 && (y.b)->_mp_size != 13)){
   //   printf("\n*************** %d **********\n", countMul);
/*         gmp_printf("\n x.a: %Zd\n", x.a);
         gmp_printf("\n y.b: %Zd\n", y.b);
    
       gmp_printf("\n x.a (h): %Zx\n", x.a);
      gmp_printf("\n y.b (h): %Zx\n", y.b);
    
    printf("\n x.a size: %d\n", (x.a)->_mp_size);
    printf("\n y.b size: %d\n", (y.b)->_mp_size);
      flag=1;
*///  }
   
    //       printf("\ntmp1 (h) \n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", ((y.parent->tmp1)->_mp_d)[12], ((y.parent->tmp1)->_mp_d)[11], ((y.parent->tmp1)->_mp_d)[10], ((y.parent->tmp1)->_mp_d)[9], ((y.parent->tmp1)->_mp_d)[8], ((y.parent->tmp1)->_mp_d)[7], ((y.parent->tmp1)->_mp_d)[6], ((y.parent->tmp1)->_mp_d)[5], ((y.parent->tmp1)->_mp_d)[4], ((y.parent->tmp1)->_mp_d)[3], ((y.parent->tmp1)->_mp_d)[2], ((y.parent->tmp1)->_mp_d)[1], ((y.parent->tmp1)->_mp_d)[0]);
    //      printf("\ntmp2 (h) \n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", ((y.parent->tmp2)->_mp_d)[12], ((y.parent->tmp2)->_mp_d)[11], ((y.parent->tmp2)->_mp_d)[10], ((y.parent->tmp2)->_mp_d)[9], ((y.parent->tmp2)->_mp_d)[8], ((y.parent->tmp2)->_mp_d)[7], ((y.parent->tmp2)->_mp_d)[6], ((y.parent->tmp2)->_mp_d)[5], ((y.parent->tmp2)->_mp_d)[4], ((y.parent->tmp2)->_mp_d)[3], ((y.parent->tmp2)->_mp_d)[2], ((y.parent->tmp2)->_mp_d)[1], ((y.parent->tmp2)->_mp_d)[0]);
    
    //   printf("\nFLAG\n");
    // }
    
    
   // mpz_set(y.parent->tmp4, x.a);
   //  mpz_set(y.parent->tmp5, y.b);
    
//mpz_mul_x86_2(y.parent->tmp1, x.a, y.b, y.parent->tmp4, y.parent->tmp5);
//   mpz_mul(y.parent->tmp1, x.a, y.b);
    
  /*  if( mpz_cmp(y.parent->tmp4,x.a)!=0 || mpz_cmp(y.parent->tmp5,y.b)!=0 )
        printf("\nFAILED\n");
    
    //  if(flag==1 || ( (y.parent->tmp1)->_mp_size !=12 && (y.parent->tmp1)->_mp_size !=13 && (y.parent->tmp1)->_mp_size !=24 && (y.parent->tmp1)->_mp_size !=25 && (y.parent->tmp1)->_mp_size !=0)){
        gmp_printf("\n res: %Zd\n", y.parent->tmp1);
    // printf("\n res (hex): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", ((y.parent->tmp3)->_mp_d)[25], ((y.parent->tmp3)->_mp_d)[24], ((y.parent->tmp3)->_mp_d)[23], ((y.parent->tmp3)->_mp_d)[22], ((y.parent->tmp3)->_mp_d)[21], ((y.parent->tmp3)->_mp_d)[20], ((y.parent->tmp3)->_mp_d)[19], ((y.parent->tmp3)->_mp_d)[18], ((y.parent->tmp3)->_mp_d)[17], ((y.parent->tmp3)->_mp_d)[16], ((y.parent->tmp3)->_mp_d)[15], ((y.parent->tmp3)->_mp_d)[14], ((y.parent->tmp3)->_mp_d)[13], ((y.parent->tmp3)->_mp_d)[12], ((y.parent->tmp3)->_mp_d)[11], ((y.parent->tmp3)->_mp_d)[10], ((y.parent->tmp3)->_mp_d)[9], ((y.parent->tmp3)->_mp_d)[8], ((y.parent->tmp3)->_mp_d)[7], ((y.parent->tmp3)->_mp_d)[6], ((y.parent->tmp3)->_mp_d)[5], ((y.parent->tmp3)->_mp_d)[4], ((y.parent->tmp3)->_mp_d)[3], ((y.parent->tmp3)->_mp_d)[2], ((y.parent->tmp3)->_mp_d)[1], ((y.parent->tmp3)->_mp_d)[0]);
     printf("\n res size: %d\n", (y.parent->tmp1)->_mp_size);
  //   printf("\n*****************************\n");
      // printf("\nFLAG\n");
   //  }
    
    */
    
  /*  printf("\n*************** %d **********\n", countMul);
    gmp_printf("\n y.a: %Zd\n", y.a);
    gmp_printf("\n x.b: %Zd\n", x.b);
    printf("\n y.a size: %d\n", (y.a)->_mp_size);
    printf("\n x.b size: %d\n", (x.b)->_mp_size);
    */
    
 //   mpz_mul_x86_2(y.parent->tmp2, y.a, x.b);
    //mpz_mul_x86(y.parent->tmp2, y.a, x.b);
    
  /*  gmp_printf("\n res: %Zd\n", y.parent->tmp2);
    printf("\n res size: %d\n", (y.parent->tmp2)->_mp_size);
    printf("\n*****************************\n");
    */
 /*   mpz_sub(y.parent->tmp3, y.parent->tmp3, y.parent->tmp1);
    mpz_add(y.parent->tmp1, y.parent->tmp1, y.parent->tmp2);
    mpz_add(y.parent->tmp3, y.parent->tmp3, y.parent->tmp2);
   */ 
 /*   printf("\n*************** %d **********\n", countMul);
    gmp_printf("\n tmp1: %Zd\n", y.parent->tmp1);
    gmp_printf("\n tmp3: %Zd\n", y.parent->tmp3);
    printf("\n tmp1 size: %d\n", (y.parent->tmp1)->_mp_size);
    printf("\n tmp3 size: %d\n", (y.parent->tmp3)->_mp_size);
   */ 
   // mpz_mod(res->a, y.parent->tmp1, y.parent->p);
   // barrett(res->a, y.parent->tmp1);
    
   // barrett_neg(res->b, y.parent->tmp3);
 //   mpz_mod(res->b, y.parent->tmp3, y.parent->p);
    
     //printf("\n*****************************\n");
    

    
 /*   printf("\n\n Outputs from mul_GF ") ; 
    gmp_printf("\n res.a: %Zd", res->a);
    gmp_printf("\n res.b: %Zd", res->b);
    printf("\n*****************", countMul) ; 
   */
 
    mpz_add(y.parent->tmp1, x.a, x.b);
    mpz_sub(y.parent->tmp2, y.b, y.a);
    mpz_mul_x86_1(y.parent->tmp3, y.parent->tmp1, y.parent->tmp2);    //could combine these 3 muls together with simd (maybe)
    mpz_mul_x86_2(y.parent->tmp1, x.a, y.b);
    mpz_mul_x86_2(y.parent->tmp2, y.a, x.b);
    mpz_sub(y.parent->tmp3, y.parent->tmp3, y.parent->tmp1);
    mpz_add(y.parent->tmp1, y.parent->tmp1, y.parent->tmp2);
    mpz_add(y.parent->tmp3, y.parent->tmp3, y.parent->tmp2);
    barrett(res->a, y.parent->tmp1);
    barrett_neg(res->b, y.parent->tmp3);

   
/*     mpz_add(y.parent->tmp1, x.a, x.b);
     mpz_sub(y.parent->tmp2, y.b, y.a);
     mpz_mul(y.parent->tmp3, y.parent->tmp1, y.parent->tmp2);
     mpz_mul(y.parent->tmp1, x.a, y.b);
     mpz_mul(y.parent->tmp2, y.a, x.b);
     mpz_sub(y.parent->tmp3, y.parent->tmp3, y.parent->tmp1);
     mpz_add(y.parent->tmp1, y.parent->tmp1, y.parent->tmp2);
     mpz_add(y.parent->tmp3, y.parent->tmp3, y.parent->tmp2);
     mpz_mod(res->a, y.parent->tmp1, y.parent->p);
     mpz_mod(res->b, y.parent->tmp3, y.parent->p);
 */  
    
    res->parent = x.parent;
    countMul++;
    
}

//12106450298025607739342710490884651867672811151897274
//x.b: 12337789618620744929184656802984297727262935111049844



//mul version where you can specify the temp variables. It will take the tmp veriables in the parent of the fourth object
void mul_GF_t(GF *res, const GF x, const GF y, mpz_t *tmp) {
    
   // mpz_t tmp1, tmp2, tmp3;
   // mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
    
    mpz_add(tmp[0], x.a, x.b);	
    
    mpz_sub(tmp[1], y.b, y.a);   
	
  //  mpz_mul_x86(tmp[2], tmp[0], tmp[1], x.parent->tmp2);
    
    mpz_mul(tmp[0], x.a, y.b);
    mpz_mul(tmp[1], y.a, x.b);
    
    mpz_sub(tmp[2], tmp[2], tmp[0]);
    mpz_add(tmp[2], tmp[2], tmp[1]);
    mpz_add(tmp[0], tmp[0], tmp[1]);
    
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_mod(res->b, tmp[2], x.parent->p);
    res->parent = x.parent;
    
 //   mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3); 
    
}

void sqr_GF(GF *res, const GF x) {
  mpz_mul(x.parent->tmp1, x.a, x.b);
    
 //   if(dieter==1)  
 //       gmp_printf("\n\ntmp1_sqr: %Zd\n\n", x.parent->tmp1);    
    
  mpz_add(x.parent->tmp1, x.parent->tmp1, x.parent->tmp1);
   
  

  mpz_add(x.parent->tmp2, x.b, x.a);
  mpz_sub(x.parent->tmp3, x.b, x.a);
  mpz_mul(x.parent->tmp2, x.parent->tmp2, x.parent->tmp3);

  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mod(res->b, x.parent->tmp2, x.parent->p);
  res->parent = x.parent;
}

//specify the GF that gives you your mpz temp variables
void sqr_GF_t(GF *res, const GF x, mpz_t *tmp) {
    
  //  mpz_t tmp1, tmp2, tmp3;
   // mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
    
    mpz_mul(tmp[0], x.a, x.b);
    mpz_add(tmp[0], tmp[0], tmp[0]);
    
    mpz_add(tmp[1], x.b, x.a);
    mpz_sub(tmp[2], x.b, x.a);
    mpz_mul(tmp[1], tmp[1], tmp[2]);
    
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_mod(res->b, tmp[1], x.parent->p);
    res->parent = x.parent;
    
   // mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3); 
}


int inv_GF(GF *res, const GF x) {
  mpz_mul(x.parent->tmp1, x.a, x.a);
  mpz_addmul(x.parent->tmp1, x.b, x.b);
  if (!mpz_invert(x.parent->tmp3, x.parent->tmp1, x.parent->p))
	  return 0;

  mpz_mul(x.parent->tmp1, x.b, x.parent->tmp3);
  mpz_neg(x.parent->tmp3, x.parent->tmp3);
  mpz_mul(x.parent->tmp2, x.a, x.parent->tmp3);

  mpz_mod(res->a, x.parent->tmp2, x.parent->p);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);

  res->parent = x.parent;
  return 1;
}

int div_GF(GF *res, const GF x, const GF y) {
 if (!inv_GF(&x.parent->GFtmp[0], y)) return 0;
  mul_GF(res, x, x.parent->GFtmp[0]);
  return 1;
}

// Miscellaneaous
int cmp_GF(const GF x, const GF y) {
  int c = mpz_cmp(x.a, y.a);
  if (c == 0) c = mpz_cmp(x.b, y.b);
  return c;
}

int is_one_GF(const GF x) {
  return (mpz_sgn(x.a) == 0) && (mpz_cmp_ui
                                 (x.b, 1) == 0);
}

int is_zero_GF(const GF x) {
  return (mpz_sgn(x.a) == 0) && (mpz_sgn(x.b) == 0);
}

void random_GF(GF *res) {
  mpz_urandomm(res->a, res->parent->state, res->parent->p);
  mpz_urandomm(res->b, res->parent->state, res->parent->p);
}

void print_GF(const GF x, char * a) {
  gmp_printf("%s: %Zd*x + %Zd\n",a, x.a, x.b);
}

//returns the characteristic of the field containing the given GF element
char * getChar( GF *res) {
	char * sa1= malloc( mpz_sizeinbase(res->parent->p, 10)+10 );
    mpz_get_str(sa1, 10, res->parent->p);
	return sa1;
}


char * getA(GF *x) {
    char * num;
    num  = malloc(sizeof(char)*1000);
    mpz_get_str(num,10,(*x).a);
    return num;
}

char * getB(  GF *x ) {
    char * num;
    num  = malloc(sizeof(char)*1000);
    mpz_get_str(num,10,(*x).b);
    return num;
}

//Returns 1 if a=b and another int otherwise
int equals(GF *a, GF *b){
    GF tmp1;
    
    GF_params *parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    init_GF( &tmp1, parent );
    sub_GF(&tmp1, *a, *b);
    return is_zero_GF(tmp1);
    
}

/***************X86 ASSEMBLY OPTIMIZATIONS***************/

void barrett_neg(mpz_t res, mpz_t a){
    
    if( a->_mp_size > 0 )
        barrett(res,a);
    else{
        a->_mp_size = (-1)*(a->_mp_size); 
       
        //modified Barrett until end of function
        mpz_mul_x86_3(res, a);
        
        switch(res->_mp_size){
                
            case 36:
                
                
                asm(
                    "shrq $0x6, %0;" 
                    
                    "movq 200(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %0;"
                    "shrq $0x6, %1;" 
                    
                    "movq 208(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %1;"
                    "shrq $0x6, %2;" 
                    
                    "movq 216(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %2;"
                    "shrq $0x6, %3;" 
                    
                    "movq 224(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %3;"
                    "shrq $0x6, %4;" 
                    
                    "movq 232(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %4;"
                    "shrq $0x6, %5;" 
                    
                    "movq 240(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %5;"
                    "shrq $0x6, %6;" 
                    
                    "movq 248(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %6;"
                    "shrq $0x6, %7;" 
                    
                    "movq 256(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %7;"
                    "shrq $0x6, %8;" 
                    
                    "movq 264(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %8;"
                    "shrq $0x6, %9;" 
                    
                    "movq 272(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %9;"
                    "shrq $0x6, %10;" 
                    
                    "movq 280(%12), %%r8;"
                    "shlq $0x3a, %%r8;" 
                    "orq %%r8, %10;"
                    "shrq $0x6, %11;" 
                    
                    :"=g" ((res->_mp_d)[24]), "=g" ((res->_mp_d)[25]), "=g" ((res->_mp_d)[26]), "=g" ((res->_mp_d)[27]), "=g" ((res->_mp_d)[28]), "=g" ((res->_mp_d)[29]), "=g" ((res->_mp_d)[30]), "=g" ((res->_mp_d)[31]), "=g" ((res->_mp_d)[32]), "=g" ((res->_mp_d)[33]), "=g" ((res->_mp_d)[34]), "=g" ((res->_mp_d)[35])
                    :"r" (res->_mp_d)
                    : "%r8"
                    
                    );

                res->_mp_d[0] =  res->_mp_d[24]; 
                res->_mp_d[1] =  res->_mp_d[25]; 
                res->_mp_d[2] =  res->_mp_d[26];  
                res->_mp_d[3] =  res->_mp_d[27]; 
                res->_mp_d[4] =  res->_mp_d[28]; 
                res->_mp_d[5] =  res->_mp_d[29];  
                res->_mp_d[6] =  res->_mp_d[30]; 
                res->_mp_d[7] =  res->_mp_d[31]; 
                res->_mp_d[8] =  res->_mp_d[32]; 
                res->_mp_d[9] =  res->_mp_d[33]; 
                res->_mp_d[10] =  res->_mp_d[34]; 
                res->_mp_d[11] =  res->_mp_d[35];  
                
                
                res->_mp_size = 12;
                
                //   printf("\nFLAG 36\n");
                
                break;
                
            case 37:
                
                //    printf("\nFLAG 37\n");
                
                asm (  "shrq $0x6, %0;" 
                     
                     
                     "movq 200(%1), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %0;"
                     
                     :"=g" ((res->_mp_d)[24])
                     : "r" (res->_mp_d)
                     :"%r8"
                     );
                
                asm (
                     
                     "shrq $0x6, %1;" 
                     
                     "movq 208(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %1;"
                     "shrq $0x6, %2;" 
                     
                     "movq 216(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %2;"
                     "shrq $0x6, %3;" 
                     
                     "movq 224(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %3;"
                     "shrq $0x6, %4;" 
                     
                     "movq 232(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %4;"
                     "shrq $0x6, %5;" 
                     
                     "movq 240(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %5;"
                     "shrq $0x6, %6;" 
                     
                     "movq 248(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %6;"
                     "shrq $0x6, %7;" 
                     
                     "movq 256(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %7;"
                     "shrq $0x6, %8;" 
                     
                     "movq 264(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %8;"
                     "shrq $0x6, %9;" 
                     
                     "movq 272(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %9;"
                     "shrq $0x6, %10;" 
                     
                     "movq 280(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %10;"
                     "shrq $0x6, %11;" 
                     
                     "movq 288(%13), %%r8;"
                     "shlq $0x3a, %%r8;" 
                     "orq %%r8, %11;"
                     "shrq $0x6, %12;" 
                     
                     
                     :"=g" ((res->_mp_d)[24]), "=g" ((res->_mp_d)[25]), "=g" ((res->_mp_d)[26]), "=g" ((res->_mp_d)[27]), "=g" ((res->_mp_d)[28]), "=g" ((res->_mp_d)[29]), "=g" ((res->_mp_d)[30]), "=g" ((res->_mp_d)[31]), "=g" ((res->_mp_d)[32]), "=g" ((res->_mp_d)[33]), "=g" ((res->_mp_d)[34]), "=g" ((res->_mp_d)[35]), "=g" ((res->_mp_d)[36])
                     :"r" (res->_mp_d)
                     : "%r8"
                     );
                
                
                res->_mp_d[0] =  res->_mp_d[24]; 
                res->_mp_d[1] =  res->_mp_d[25]; 
                res->_mp_d[2] =  res->_mp_d[26];  
                res->_mp_d[3] =  res->_mp_d[27]; 
                res->_mp_d[4] =  res->_mp_d[28]; 
                res->_mp_d[5] =  res->_mp_d[29];  
                res->_mp_d[6] =  res->_mp_d[30]; 
                res->_mp_d[7] =  res->_mp_d[31]; 
                res->_mp_d[8] =  res->_mp_d[32]; 
                res->_mp_d[9] =  res->_mp_d[33]; 
                res->_mp_d[10] =  res->_mp_d[34]; 
                res->_mp_d[11] =  res->_mp_d[35];
                res->_mp_d[12] =  res->_mp_d[36];
                
                //   res->_mp_size = 13;
                break;
                
            default:
                res->_mp_size=0;
                
                
        }

        mpz_set(tmp,res);
        mpz_mul_x86_4(res, tmp);
        

        mpz_sub(res, a, res);
        
        if( mpz_cmp(res, prime) > 0 ){
            

            //2p-res
           
            
            asm(
                
                "movq	0(%7), %%r11;"
                "movq	$0xfffffffffffffffe, %%r8;"
                "subq	%%r11, %%r8;"
                "movq %%r8, %0;"
                
                "movq	8(%7), %%r11;"
                     "movq	$0xffffffffffffffff, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %1;"
                
                "movq	16(%7), %%r11;"
                  "movq	$0xffffffffffffffff, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %2;"
                
                "movq	24(%7), %%r11;"
                   "movq	$0xffffffffffffffff, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %3;"
                
                "movq	32(%7), %%r11;"
                    "movq	$0xffffffffffffffff, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %4;"
                
                "movq	40(%7), %%r11;"
                   "movq	$0xffffffffffffffff, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %5;"
                              
                "movq	48(%7), %%r11;"
                "movq	$0xe00ecd34b9d12c8f, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %6;"
                
                             
                : "=&c" ((res->_mp_d)[0]), "=&d" ((res->_mp_d)[1]), "=&a" ((res->_mp_d)[2]), "=&m" ((res->_mp_d)[3]), "=&m" ((res->_mp_d)[4]), "=&m" ((res->_mp_d)[5]), "=&m" ((res->_mp_d)[6])
                : "r" ((res->_mp_d))
                : "%r8", "%r11"//, "memory", "cc"
                );
            
            asm(
                
                "movq	56(%6), %%r11;"
                "movq	$0x5bc01b22908a09f3, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %0;"
                //            
                "movq	64(%6), %%r11;"
                "movq	$0x12f3aae184890dc7, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %1;"
                
                "movq	72(%6), %%r11;"
                "movq	$0x1775c6cf3498e04a, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %2;"
                //          
                "movq	80(%6), %%r11;"
                "movq	$0x40d500b53ed01169, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %3;"
                
                "movq	88(%6), %%r11;"
                "movq	$0xcd3d7df50ff57bf5, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %4;"
                
                //  "xor %%r11, %%r11;"
                //  "xor %%r8, %%r8;"
                
                "movq	96(%6), %%r11;"
                "movq	$0xb, %%r8;"
                "sbbq	%%r11, %%r8;"
                "movq %%r8, %5;"
                
                
                : "=&m" ((res->_mp_d)[7]), "=&m" ((res->_mp_d)[8]), "=&m" ((res->_mp_d)[9]), "=&m" ((res->_mp_d)[10]), "=&m" ((res->_mp_d)[11]), "=b" ((res->_mp_d)[12])
                : "r" ((res->_mp_d))
                : "%r8", "%r11"//, "memory", "cc"
                );
            
            
       //     if( (res->_mp_d)[12]==0 )
        //        res->_mp_size = 12;
            
              return;
            
            
        }
          

        if( res->_mp_size==0 ){
           mpz_set(res,prime);
        }
        
        //next two assembly blocks add the prime to res (which is negative so thats why were doing a subtraction)
        asm(
            
            "movq	$0xffffffffffffffff, %%r8;"
            "movq	0(%7), %%r11;"
            "subq	%%r11, %%r8;"
            "movq %%r8, %0;"
            
            "movq	8(%7), %%r11;"
                 "movq	$0xffffffffffffffff, %%r8;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %1;"
            
            "movq	16(%7), %%r11;"
              "movq	$0xffffffffffffffff, %%r8;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %2;"
            
            "movq	24(%7), %%r11;"
               "movq	$0xffffffffffffffff, %%r8;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %3;"
            
            "movq	32(%7), %%r11;"
                "movq	$0xffffffffffffffff, %%r8;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %4;"
            
            "movq	40(%7), %%r11;"
               "movq	$0xffffffffffffffff, %%r8;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %5;"
            //              
            "movq	$0xf007669a5ce89647, %%r8;"
            "movq	48(%7), %%r11;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %6;"
            
            
            : "=&c" ((res->_mp_d)[0]), "=&d" ((res->_mp_d)[1]), "=&a" ((res->_mp_d)[2]), "=&S" ((res->_mp_d)[3]), "=&o" ((res->_mp_d)[4]), "=&g" ((res->_mp_d)[5]), "=&g" ((res->_mp_d)[6])
            : "r" ((res->_mp_d))
            : "%r8", "%r11"//, "memory", "cc"
            );
        
        asm(
        
            
            "movq	$0xade00d91484504f9, %%r8;"
            "movq	56(%6), %%r11;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %0;"
            //            
            "movq	$0x0979d570c24486e3, %%r8;"
            "movq	64(%6), %%r11;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %1;"
            
            "movq	$0x8bbae3679a4c7025, %%r8;"
            "movq	72(%6), %%r11;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %2;"
            //          
            "movq	$0xa06a805a9f6808b4, %%r8;"
            "movq	80(%6), %%r11;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %3;"
            
            "movq	$0xe69ebefa87fabdfa, %%r8;"
            "movq	88(%6), %%r11;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %4;"
            
            //  "xor %%r11, %%r11;"
            //  "xor %%r8, %%r8;"
            
            "movq	$0x5, %%r8;"
            "movq	96(%6), %%r11;"
            "sbbq	%%r11, %%r8;"
            "movq %%r8, %5;"
            
            
            :  "=&g" ((res->_mp_d)[7]), "=&g" ((res->_mp_d)[8]), "=&m" ((res->_mp_d)[9]), "=&g" ((res->_mp_d)[10]), "=&g" ((res->_mp_d)[11]), "=b" ((res->_mp_d)[12])
            : "r" ((res->_mp_d))
            : "%r8", "%r11"//, "memory", "cc"
            );
        

        if( (res->_mp_d)[12]!=0 )
            res->_mp_size = 13; 
        else
            res->_mp_size = 12; 
    }
    
    
}

void barrett(mpz_t res, mpz_t a){
/*    
    mpz_t tmp1, tmp2, tmp3;
    
    mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
  */
   
   
/*    printf("\n*************** %d **********\n", countBarrett);
    gmp_printf("\n a: %Zd\n", a);
    gmp_printf("\n mB: %Zd\n", mB);
    printf("\n a size: %d\n", a->_mp_size);
    printf("\n mB size: %d\n", mB->_mp_size);
   */ 
    
   // mpz_mul(res, a, mB);
      mpz_mul_x86_3(res, a);
    
  //  if(mpz_cmp(res,tmp3) != 0)
   //     printf("\nFAILED\n");
    
 /*   gmp_printf("\n res: %Zd\n", res);
    printf("\n res (limbs): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx \n\n", (res->_mp_d)[36], (res->_mp_d)[35], (res->_mp_d)[24], (res->_mp_d)[34], (res->_mp_d)[33], (res->_mp_d)[32], (res->_mp_d)[31], (res->_mp_d)[30], (res->_mp_d)[29], (res->_mp_d)[28], (res->_mp_d)[27], (res->_mp_d)[26], (res->_mp_d)[25], (res->_mp_d)[24], (res->_mp_d)[23], (res->_mp_d)[22], (res->_mp_d)[21], (res->_mp_d)[20], (res->_mp_d)[19], (res->_mp_d)[18], (res->_mp_d)[17], (res->_mp_d)[16], (res->_mp_d)[15], (res->_mp_d)[14], (res->_mp_d)[13], (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
    gmp_printf("\n res (h): %Zx\n", res);
    printf("\n res size: %d\n", res->_mp_size);
    printf("\n*****************************\n", countBarrett);
   */ 
      // the next ~100 lines implement a special purpose division algo based on multiple precision bit shifting (were dividing by a power of 2)
    // the commented out assembly code does not seem to be as fast as the c instructions (not commented out) it is desgined to replace
    
   
   // if( res->_mp_size != 37 && res->_mp_size != 0 && res->_mp_size != 36 && res->_mp_size != 25 && res->_mp_size != 24 ){
 /*      printf("\n*************** %d **************\n", countBarrett);
        gmp_printf("\n res: %Zd\n", res);
    gmp_printf("\n res: %Zx\n", res);
        printf("\n res (limbs): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx \n\n", (res->_mp_d)[36], (res->_mp_d)[35], (res->_mp_d)[24], (res->_mp_d)[34], (res->_mp_d)[33], (res->_mp_d)[32], (res->_mp_d)[31], (res->_mp_d)[30], (res->_mp_d)[29], (res->_mp_d)[28], (res->_mp_d)[27], (res->_mp_d)[26], (res->_mp_d)[25], (res->_mp_d)[24], (res->_mp_d)[23], (res->_mp_d)[22], (res->_mp_d)[21], (res->_mp_d)[20], (res->_mp_d)[19], (res->_mp_d)[18], (res->_mp_d)[17], (res->_mp_d)[16], (res->_mp_d)[15], (res->_mp_d)[14], (res->_mp_d)[13], (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
      //  gmp_printf("\n res (h): %Zx\n", res);
        printf("\n res size: %d\n", res->_mp_size); 
    
*/
 //   }
    
    //size is available in the vaieties noted above. must left shift things by 58, then or them as in other assembly code. right shift things after that
    //by 6. need (at most) size-24 blocks of the highest limbs for each case. i.e we start at (res->_mp_d)[24] for each case and only consider up to 
    //(res->_mp_d)[size-1] at most
   
  //  mpz_set(tmp1,res);
  //  mpz_set(tmp2,kF);
    
    switch(res->_mp_size){
            
        case 36:
            
            
            asm(
                "shrq $0x6, %0;" 
                
                "movq 200(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %0;"
                "shrq $0x6, %1;" 
                
                "movq 208(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %1;"
                "shrq $0x6, %2;" 
                
                "movq 216(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %2;"
                "shrq $0x6, %3;" 
                
                "movq 224(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %3;"
                "shrq $0x6, %4;" 
                
                "movq 232(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %4;"
                "shrq $0x6, %5;" 
                
                "movq 240(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %5;"
                "shrq $0x6, %6;" 
                
                "movq 248(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %6;"
                "shrq $0x6, %7;" 
                
                "movq 256(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %7;"
                "shrq $0x6, %8;" 
                
                "movq 264(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %8;"
                "shrq $0x6, %9;" 
                
                "movq 272(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %9;"
                "shrq $0x6, %10;" 
                
                "movq 280(%12), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %10;"
                "shrq $0x6, %11;" 
            
                :"=g" ((res->_mp_d)[24]), "=g" ((res->_mp_d)[25]), "=g" ((res->_mp_d)[26]), "=g" ((res->_mp_d)[27]), "=g" ((res->_mp_d)[28]), "=g" ((res->_mp_d)[29]), "=g" ((res->_mp_d)[30]), "=g" ((res->_mp_d)[31]), "=g" ((res->_mp_d)[32]), "=g" ((res->_mp_d)[33]), "=g" ((res->_mp_d)[34]), "=g" ((res->_mp_d)[35])
                :"r" (res->_mp_d)
                : "%r8"
                
            );
            
       //     int *point = &(res->_mp_d[24]); 
            
        //    res->_mp_d =  &(res->_mp_d[24]); 
            
            res->_mp_d[0] =  res->_mp_d[24]; 
            res->_mp_d[1] =  res->_mp_d[25]; 
            res->_mp_d[2] =  res->_mp_d[26];  
            res->_mp_d[3] =  res->_mp_d[27]; 
            res->_mp_d[4] =  res->_mp_d[28]; 
            res->_mp_d[5] =  res->_mp_d[29];  
            res->_mp_d[6] =  res->_mp_d[30]; 
            res->_mp_d[7] =  res->_mp_d[31]; 
            res->_mp_d[8] =  res->_mp_d[32]; 
            res->_mp_d[9] =  res->_mp_d[33]; 
            res->_mp_d[10] =  res->_mp_d[34]; 
            res->_mp_d[11] =  res->_mp_d[35];  
            
            
            res->_mp_size = 12;
            
         //   printf("\nFLAG 36\n");
            
            break;
            
        case 37:
            
         //    printf("\nFLAG 37\n");
            
            asm (  "shrq $0x6, %0;" 
                 
                 
                 "movq 200(%1), %%r8;"
                 "shlq $0x3a, %%r8;" 
                 "orq %%r8, %0;"
                 
                 :"=g" ((res->_mp_d)[24])
                 : "r" (res->_mp_d)
                 :"%r8"
                 );
            
            asm (
                
                "shrq $0x6, %1;" 
                
                "movq 208(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %1;"
                "shrq $0x6, %2;" 
                
                "movq 216(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %2;"
                "shrq $0x6, %3;" 
                
                "movq 224(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %3;"
                "shrq $0x6, %4;" 
                
                "movq 232(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %4;"
                "shrq $0x6, %5;" 
                
                "movq 240(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %5;"
                "shrq $0x6, %6;" 
                
                "movq 248(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %6;"
                "shrq $0x6, %7;" 
                
                "movq 256(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %7;"
                "shrq $0x6, %8;" 
                
                "movq 264(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %8;"
                "shrq $0x6, %9;" 
                
                "movq 272(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %9;"
                "shrq $0x6, %10;" 
                
                "movq 280(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %10;"
                "shrq $0x6, %11;" 
                
                "movq 288(%13), %%r8;"
                "shlq $0x3a, %%r8;" 
                "orq %%r8, %11;"
                "shrq $0x6, %12;" 
                
                
                :"=g" ((res->_mp_d)[24]), "=g" ((res->_mp_d)[25]), "=g" ((res->_mp_d)[26]), "=g" ((res->_mp_d)[27]), "=g" ((res->_mp_d)[28]), "=g" ((res->_mp_d)[29]), "=g" ((res->_mp_d)[30]), "=g" ((res->_mp_d)[31]), "=g" ((res->_mp_d)[32]), "=g" ((res->_mp_d)[33]), "=g" ((res->_mp_d)[34]), "=g" ((res->_mp_d)[35]), "=g" ((res->_mp_d)[36])
                :"r" (res->_mp_d)
                : "%r8"
                );
            

            
         //   res->_mp_d =  &(res->_mp_d[24]); 
            
            //              printf("\n res 24: %lx\n", (res->_mp_d)[24]);
            
            
             res->_mp_d[0] =  res->_mp_d[24]; 
             res->_mp_d[1] =  res->_mp_d[25]; 
             res->_mp_d[2] =  res->_mp_d[26];  
             res->_mp_d[3] =  res->_mp_d[27]; 
             res->_mp_d[4] =  res->_mp_d[28]; 
             res->_mp_d[5] =  res->_mp_d[29];  
             res->_mp_d[6] =  res->_mp_d[30]; 
             res->_mp_d[7] =  res->_mp_d[31]; 
             res->_mp_d[8] =  res->_mp_d[32]; 
             res->_mp_d[9] =  res->_mp_d[33]; 
             res->_mp_d[10] =  res->_mp_d[34]; 
             res->_mp_d[11] =  res->_mp_d[35];
             res->_mp_d[12] =  res->_mp_d[36];
            
            //   res->_mp_size = 13;
            break;
  
        default:
            res->_mp_size=0;
            
        
    }
     //  mpz_fdiv_q(tmp3,tmp1,tmp2); 
      
 /*  
    gmp_printf("\n ass res: %Zd\n", res);
    printf("\n ass res size: %d\n", res->_mp_size);
    printf("\n ass res (limbs): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx \n\n", (res->_mp_d)[13], (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);

    gmp_printf("\n mpz res: %Zd\n", tmp3);
    printf("\n mpz res size: %d\n", tmp3->_mp_size);
    printf("\n mpz res (limbs): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx \n\n", (tmp3->_mp_d)[13], (tmp3->_mp_d)[12], (tmp3->_mp_d)[11], (tmp3->_mp_d)[10], (tmp3->_mp_d)[9], (tmp3->_mp_d)[8], (tmp3->_mp_d)[7], (tmp3->_mp_d)[6], (tmp3->_mp_d)[5], (tmp3->_mp_d)[4], (tmp3->_mp_d)[3], (tmp3->_mp_d)[2], (tmp3->_mp_d)[1], (tmp3->_mp_d)[0]);

    printf("\n*****************************\n", countBarrett);
*/
    
  //  mpz_t tmp;    
   // mpz_init2(tmp,13*64);// mpz_init(tmp2); mpz_init(tmp3);
/*    
    mpz_set(tmp1,res);
    mpz_set(tmp2,prime);
    */
/*    int flag=0;    
 //   if( ((res)->_mp_size!=0 && (res)->_mp_size!=12 && (res)->_mp_size!=13) ){
        printf("\n*************** %d **********\n", countBarrett);
        gmp_printf("\n r: %Zd\n", res);
        gmp_printf("\n p: %Zd\n", prime);
        printf("\n r size: %d\n", res->_mp_size);
        printf("\n p size: %d\n", prime->_mp_size);
        flag=1;
 //   }
  */  
 //   mpz_mul(res, res, prime);
    mpz_set(tmp,res);
    mpz_mul_x86_4(res, tmp);
    
  /*  if( mpz_cmp(tmp3,res) != 0 ){
        printf("\nFAILED\n");
        gmp_printf("\nmpz_res: %Zd \n", tmp3);
        printf("\n mpz_res limbs: %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (tmp3->_mp_d)[25], (tmp3->_mp_d)[24], (tmp3->_mp_d)[23], (tmp3->_mp_d)[22], (tmp3->_mp_d)[21], (tmp3->_mp_d)[20], (tmp3->_mp_d)[19], (tmp3->_mp_d)[18], (tmp3->_mp_d)[17], (tmp3->_mp_d)[16], (tmp3->_mp_d)[15], (tmp3->_mp_d)[14], (tmp3->_mp_d)[13], (tmp3->_mp_d)[12], (tmp3->_mp_d)[11], (tmp3->_mp_d)[10], (tmp3->_mp_d)[9], (tmp3->_mp_d)[8], (tmp3->_mp_d)[7], (tmp3->_mp_d)[6], (tmp3->_mp_d)[5], (tmp3->_mp_d)[4], (tmp3->_mp_d)[3], (tmp3->_mp_d)[2], (tmp3->_mp_d)[1], (tmp3->_mp_d)[0]);

    }
    */
 //   printf("\n res (limbs): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx \n\n", (res->_mp_d)[36], (res->_mp_d)[35], (res->_mp_d)[24], (res->_mp_d)[34], (res->_mp_d)[33], (res->_mp_d)[32], (res->_mp_d)[31], (res->_mp_d)[30], (res->_mp_d)[29], (res->_mp_d)[28], (res->_mp_d)[27], (res->_mp_d)[26], (res->_mp_d)[25], (res->_mp_d)[24], (res->_mp_d)[23], (res->_mp_d)[22], (res->_mp_d)[21], (res->_mp_d)[20], (res->_mp_d)[19], (res->_mp_d)[18], (res->_mp_d)[17], (res->_mp_d)[16], (res->_mp_d)[15], (res->_mp_d)[14], (res->_mp_d)[13], (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
  //  if(flag || (res->_mp_size != 25 && res->_mp_size != 0 && res->_mp_size != 24)){
    /*    printf("\n res limbs: %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[25], (res->_mp_d)[24], (res->_mp_d)[23], (res->_mp_d)[22], (res->_mp_d)[21], (res->_mp_d)[20], (res->_mp_d)[19], (res->_mp_d)[18], (res->_mp_d)[17], (res->_mp_d)[16], (res->_mp_d)[15], (res->_mp_d)[14], (res->_mp_d)[13], (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
        gmp_printf("\n res (h): %Zd\n", res);
        printf("\n res size: %d\n", res->_mp_size);
        printf("\n*****************************\n", countBarrett);
   // }
    */    
    mpz_sub(res, a, res);
    
    if( mpz_cmp(res, prime) > 0 ){
        
 //   mpz_sub(res,res,prime); }

        
            /*   printf("\n*************** %d **********\n", countBarrett);
             gmp_printf("\nres: %Zd \n", res);
             printf("\nres: %d \n", res->_mp_size);
             */ 
        
     /*       mpz_t tmp1;
            mpz_init(tmp1);
            mpz_set(tmp1,res);
            
            mpz_t tmp33, tmp14, tmp2;
            mpz_init(tmp2); mpz_init(tmp33); mpz_init(tmp14);
            mpz_set(tmp1,res);
            mpz_set(tmp2,prime);
       */ 
          
        
          
     
            /*    if(res->_mp_size !=12 && res->_mp_size !=13){
             gmp_printf("\nres: %Zd \n", res);
             printf("\nres: %d \n", res->_mp_size);
             printf("\n*****************************\n", countBarrett);
             } 
                     
            /*     if(mpz_cmp(res,prime)==0){
             mpz_set_ui(res,0);
             countBarrett++;
             return;
             }
             */   
          
        
 /*       
            printf("\n*************** %d **************\n", countBarrett);
        
            gmp_printf("\nres: %Zd \n", res);
            printf("\nres size: %d \n", res->_mp_size);
            gmp_printf("\nprime: %Zd \n", prime);
   */     
          //  if( (res->_mp_size)==12 )
       //         (res->_mp_d)[12]=0;
        
            asm(
            
                        "movq	0(%13), %%r8;"
                         "movq	$0xffffffffffffffff, %%r11;"
                         "subq	%%r11, %%r8;"
                         "movq %%r8, %0;"
            
                         "movq	8(%13), %%r8;"
                         //     "movq	$0xffffffffffffffff, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %1;"
            
                         "movq	16(%13), %%r8;"
                         //  "movq	$0xffffffffffffffff, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %2;"
                         
                         "movq	24(%13), %%r8;"
                         //   "movq	$0xffffffffffffffff, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %3;"
            
                         "movq	32(%13), %%r8;"
                         //    "movq	$0xffffffffffffffff, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %4;"
            
                         "movq	40(%13), %%r8;"
                         //   "movq	$0xffffffffffffffff, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %5;"
      //              
                         "movq	48(%13), %%r8;"
                         "movq	$0xf007669a5ce89647, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %6;"
                     
                         "movq	56(%13), %%r8;"
                         "movq	$0xade00d91484504f9, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %7;"
         //            
                         "movq	64(%13), %%r8;"
                         "movq	$0x0979d570c24486e3, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %8;"
                     
                         "movq	72(%13), %%r8;"
                         "movq	$0x8bbae3679a4c7025, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %9;"
           //          
                         "movq	80(%13), %%r8;"
                         "movq	$0xa06a805a9f6808b4, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %10;"
                    
                         "movq	88(%13), %%r8;"
                         "movq	$0xe69ebefa87fabdfa, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %11;"
                     
                         //  "xor %%r11, %%r11;"
                         //  "xor %%r8, %%r8;"
                     
                         "movq	96(%13), %%r8;"
                         "movq	$0x5, %%r11;"
                         "sbbq	%%r11, %%r8;"
                         "movq %%r8, %12;"
                     

                         : "=&c" ((res->_mp_d)[0]), "=&d" ((res->_mp_d)[1]), "=&a" ((res->_mp_d)[2]), "=&m" ((res->_mp_d)[3]), "=&m" ((res->_mp_d)[4]), "=&m" ((res->_mp_d)[5]), "=&m" ((res->_mp_d)[6]), "=&m" ((res->_mp_d)[7]), "=&m" ((res->_mp_d)[8]), "=&m" ((res->_mp_d)[9]), "=&m" ((res->_mp_d)[10]), "=&m" ((res->_mp_d)[11]), "=b" ((res->_mp_d)[12])
                         : "r" ((res->_mp_d))
                         : "%r8", "%r11"//, "memory", "cc"
                         );
        
                        if( (res->_mp_d)[12]==0 )
                             res->_mp_size = 12;
    
        
    }else if( mpz_cmp(res,prime)==0 )
               mpz_set_ui(res,0);

    
    
/*      mpz_sub(tmp33, tmp1, tmp2);
    
    
           printf("\nass res limbs (h) : %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx \n", (res->_mp_d)[0], (res->_mp_d)[1], (res->_mp_d)[2], (res->_mp_d)[3], (res->_mp_d)[4], (res->_mp_d)[5], (res->_mp_d)[6], (res->_mp_d)[7], (res->_mp_d)[8], (res->_mp_d)[9], (res->_mp_d)[10], (res->_mp_d)[11], (res->_mp_d)[12]);
           printf("\nmpz res limbs (h) : %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx \n", (tmp33->_mp_d)[0], (tmp33->_mp_d)[1], (tmp33->_mp_d)[2], (tmp33->_mp_d)[3], (tmp33->_mp_d)[4], (tmp33->_mp_d)[5], (tmp33->_mp_d)[6], (tmp33->_mp_d)[7], (tmp33->_mp_d)[8], (tmp33->_mp_d)[9], (tmp33->_mp_d)[10], (tmp33->_mp_d)[11], (tmp33->_mp_d)[12]);
    
    
           gmp_printf("\nres: %Zd \n", res);
     printf("\nres size: %d \n", res->_mp_size);
     gmp_printf("\nmpz res: %Zd \n", tmp33);
     printf("\nmpz res size: %d \n", tmp33->_mp_size);
     printf("\n*****************************\n", countBarrett);
     
    */
    
    
    
    
    countBarrett++;

}

//testing SIMD instructions
void mpz_mul_x86_0(mpz_t res, mpz_t ain, mpz_t bin){
    
        asm (
             ";"
            // "movaps 0(%6), %%xmm1;"
            // "movaps %%xmm1, %0;"
             
            /* "movq 0(%6), %%rax;" 
             "mulq 0(%7);"
             "movq %%rax, %0;"    //AO*BO  
             "movq %%rdx, %%r8;"    
            
            //   "xorq %%rdx, %%rdx;"
             
             "movq 0(%6), %%rax;"
             "mulq 8(%7);"              //A0*B1
             "addq %%rax, %%r8;"     
             "movq %%rdx, %%r9;"  
             "adcq $0, %%r9;"
             
             "xorq %%r10, %%r10;"  
             
             "movq 8(%6), %%rax;"  //problem seg fault line
             "mulq 0(%7);"         //A1*B0
             "addq %%rax, %%r8;" 
             "movq %%r8, %1;"  
             "adcq %%rdx,%%r9;"    
             "adcq $0, %%r10;"  
           
             "xorq %%r8, %%r8;"
             "movq 0(%6), %%rax;" 
             "mulq  16(%7);"        //A0*B2
             "addq %%rax, %%r9;"
             "adcq %%rdx, %%r10;"
             "adcq $0, %%r8;"
             
             "movq 8(%6), %%rax;"
             "mulq 8(%7);"              //A1*b1
             "addq %%rax, %%r9;"
             "adcq %%rdx, %%r10;"
             "adcq $0, %%r8;"
             
             "movq 16(%6), %%rax;"
             "mulq 0(%7);"            
             "addq %%rax, %%r9;"      
             "movq %%r9, %2;"       //a2*b0
             "adcq %%rdx, %%r10;"
             "adcq $0, %%r8;"
             
             "xorq %%r9, %%r9;"         
             "movq 8(%6), %%rax;"
             "mulq 16(%7);"            
             "addq %%rax, %%r10;"      //a1*b2 nd
             "adcq %%rdx, %%r8;"
             "adcq $0, %%r9;"
            
             "movq 16(%6), %%rax;"
             "mulq 8(%7);"            
             "addq %%rax, %%r10;"          //a2*b1 nd
             "movq %%r10, %3;"  
             "adcq %%rdx, %%r8;"
             "adcq $0, %%r9;"
             
             "movq 16(%6), %%rax;"
             "mulq 16(%7);"            
             "addq %%rax, %%r8;" //a2*b2 nd
             "movq %%r8, %4;" 
             "adcq %%rdx, %%r9;"         
             "adcq $0, %%r9;"
             "movq %%r9, %5;"
      */

         :  "=&g" ((res->_mp_d)[0]), "=&g" ((res->_mp_d)[1]), "=&g" ((res->_mp_d)[2]), "=&g" ((res->_mp_d)[3]), "=&g" ((res->_mp_d)[4]), "=&g" ((res->_mp_d)[5])
         : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))
         : "%rax", "%rdx", "%r10", "%r9", "%r8", "xmm0", "xmm1"
         );
    
    printf("\n\n\n IN FUNC ass res:  %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);

    
    res->_mp_size = 6;

}
 
 
    
void mpz_mul_x86_1(mpz_t res, mpz_t ain, mpz_t bin){

    negFlag=1;
    unsigned long int ra=0;
    unsigned long int rb=0;
    

   
    switch (ain->_mp_size) {
        
        case 0:
            res->_mp_size=0;
            return;
            
        case 1:
            mpz_set(res,bin);
            return;
        
        case 12:
            (ain->_mp_d)[12]=0;

            break;
            
        case -12:
            (ain->_mp_d)[12]=0;
            
            negFlag=-1;
            break;
        
        case -13:
            
            negFlag=-1;
            break;
            
    }
    
    
    switch (bin->_mp_size) {
       
        case 0:
            res->_mp_size=0;
            return;
            
        case 1:
            mpz_set(res,ain);
            return;
            
        case 12:
            (bin->_mp_d)[12]=0;
            
            
        //    if( res->_mp_alloc<25)
        //        _mpz_realloc(res,25);
            
            break;
            
        case -12:
            (bin->_mp_d)[12]=0;
            
            
        //    if( res->_mp_alloc<25)
        //        _mpz_realloc(res,25);
            
            negFlag=-1*negFlag;
            break;
        
        case -13:
            
          //  if( res->_mp_alloc<25)
          //      _mpz_realloc(res,25);
            
            negFlag=-1*negFlag;
            break;
    }
    
    
    
   //  printf("IN FUNC ain:\n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (ain->_mp_d)[12], (ain->_mp_d)[11], (ain->_mp_d)[10], (ain->_mp_d)[9], (ain->_mp_d)[8], (ain->_mp_d)[7], (ain->_mp_d)[6], (ain->_mp_d)[5], (ain->_mp_d)[4], (ain->_mp_d)[3], (ain->_mp_d)[2], (ain->_mp_d)[1], (ain->_mp_d)[0]);
  //  printf("IN FUNC bin:\n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (bin->_mp_d)[12], (bin->_mp_d)[11], (bin->_mp_d)[10], (bin->_mp_d)[9], (bin->_mp_d)[8], (bin->_mp_d)[7], (bin->_mp_d)[6], (bin->_mp_d)[5], (bin->_mp_d)[4], (bin->_mp_d)[3], (bin->_mp_d)[2], (bin->_mp_d)[1], (bin->_mp_d)[0]);

    
      asm (            
        "movq 0(%9), %%rax;" 
        "mulq 0(%10);"
        "movq %%rax, %0;"    
        "movq %%rdx, %7;"           //A0*B0
                                      //0
        
        "xorq %%r10, %%r10;" 
        
        "movq 8(%9), %%rax;"      //problem line
        "mulq 0(%10);"              
        "addq %%rax, %7;"     
        "movq %%rdx, %8;"  
        "adcq $0, %8;"              //A1*B0
        
        "movq 0(%9), %%rax;"  
        "mulq 8(%10);"         
        "addq %%rax, %7;" 
        "movq %7, %1;"  
        "adcq %%rdx,%8;"    
        "adcq $0, %%r10;"                //A0*B1
                                        //1
        
        "xorq %7, %7;" 
        
        "movq 0(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"           //A0*B2
        
        "movq 8(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A1*B1
           
        "movq 16(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %8;"    
        "movq %8, %2;" 
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"            //A2*B0
                                    //2
        "xorq %8, %8;"  
        
        "movq 24(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"              //A3*B0
        
        "movq 0(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"            //A0*B3
        
        "movq 16(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A2*B1
        
        "movq 8(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %%r10;"   
        "movq %%r10, %3;" 
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A1*B2
                                //3
        "xorq %%r10, %%r10;" 
        
        "movq 32(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"            //A4*B0
        
        "movq 0(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"            //A0*B4
        
        "movq 24(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A3*B1
        
        "movq 8(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A1*B3
        
        "movq 16(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %7;"  
        "movq %7, %4;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A2*B2
                                //4
        "xor %7, %7;"
                    
        "movq 40(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"            //A5*B0
       
        "movq 0(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A0*B5
         
        "movq 32(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A4*B1
        
        "movq 8(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A1*B4
        
        "movq 24(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A3*B2
        
        "movq 16(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %8;"
        "movq %8, %5;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"            //A2*B3
                                    //5
        "xor %8, %8;"

        "movq 48(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A6*B0
         
        "movq 0(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A0*B6
        
        "movq 40(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A5*B1
        
        "movq 8(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A1*B5
      
        "movq 32(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"            //A4*B2
       
        "movq 16(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"            //A2*B4
         
        "movq 24(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %%r10;"
        "movq %%r10, %6;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"  
         
    //    "movq %%r8, %7;" 
     //   "movq %%r9, %8;" 
        
                                //A3*B3
                                //6       //break
        :  "=&g" ((res->_mp_d)[0]), "=&g" ((res->_mp_d)[1]), "=&g" ((res->_mp_d)[2]), "=&g" ((res->_mp_d)[3]), "=&g" ((res->_mp_d)[4]), "=&g" ((res->_mp_d)[5]), "=&g" ((res->_mp_d)[6]), "+r" (ra), "+r" (rb)//,  "=&X" ((res->_mp_d)[2]), "=&X" ((res->_mp_d)[3]), "=&X" ((res->_mp_d)[4]),  "=&X" ((res->_mp_d)[5]), "=&X" ((res->_mp_d)[6]), "=&X" ((res->_mp_d)[7]),  "=&X" ((res->_mp_d)[8]), "=&X" ((res->_mp_d)[9]), "=&X" ((res->_mp_d)[10]), "=&X" ((res->_mp_d)[11]), "=&X" ((res->_mp_d)[12]), "=&X" ((res->_mp_d)[13]),  "=&X" ((res->_mp_d)[14]), "=&X" ((res->_mp_d)[15]), "=&X" ((res->_mp_d)[16]),  "=&X" ((res->_mp_d)[17]), "=&X" ((res->_mp_d)[18]), "=&X" ((res->_mp_d)[19]),  "=&X" ((res->_mp_d)[20]), "=&X" ((res->_mp_d)[21]), "=&X" ((res->_mp_d)[22]), "=&X" ((res->_mp_d)[23]), "=&X" ((res->_mp_d)[24]), "=&X" ((res->_mp_d)[25])
        : "r" ((ain->_mp_d)), "r" ((bin->_mp_d)) // "X" ((ain->_mp_d)[11]), "X" ((bin->_mp_d)[11]), "X" ((ain->_mp_d)[12]), "X" ((bin->_mp_d)[12])
        : "%rax", "%rdx", "%r10"
        );
    
    asm (
                 
     //   "movq %11, %%r8;" 
     //   "movq %12, %%r9;" 
    
                    
        //"xor %%r8, %%r8;"
        //"xor %%r9, %%r9;"
        "xor %%r10, %%r10;"

        
        "movq 56(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A7*B0
         
        "movq 0(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A0*B7
         
        "movq 48(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A6*B1
         
        "movq 8(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A1*B6
        
        "movq 40(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A5*B2
      
        "movq 16(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"        //A2*B5
        
        "movq 32(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A4*B3
        
        "movq 24(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %7;"  
        "movq %7, %0;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A3*B4
                                //7
        "xor %7, %7;"

        "movq 64(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"       //A8*B0
          
        "movq 0(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"             //A0*B8
        
        "movq 56(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A7*B1
         
        "movq 8(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A1*B7
        
        "movq 48(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A6*B2
         
        "movq 16(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A2*B6
         
        "movq 40(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A5*B3
       
        "movq 24(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %8;"            
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"            //A3*B5
         
        "movq 32(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %8;" 
        "movq %8, %1;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"            //A4*B4
                                    //8 
        
        "xor %8, %8;"
        
        "movq 72(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"         //A9*B0
         
        "movq 0(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A0*B9
        
        "movq 64(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"          //A8*B1
        
        "movq 8(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"         //A1*B8
          
        "movq 56(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"  //A7*B2
        
        "movq 16(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"         //A2*B7
        
        "movq 48(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A6*B3
        
        "movq 24(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A3*B6
         
        "movq 40(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %%r10;"            
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A5*B4
         
        "movq 32(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %%r10;"
        "movq %%r10, %2;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A4*B5
                                //9  
        "xor %%r10, %%r10;"

        "movq 80(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"         //A10*B0
        
        "movq 0(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"        //A0*B10
        
        "movq 72(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"         //A9*B1
        
        "movq 8(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"         //A1*B9
         
        "movq 64(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A8*B2
        
        "movq 16(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A2*B8
       
        "movq 56(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"          //A7*B3
         
        "movq 24(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A3*B7
        
        "movq 48(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A6*B4
         
        "movq 32(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %7;"            
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"       //A4*B6
      
        "movq 40(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %7;"
        "movq %7, %3;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"         //A5*B5
                                //10
        "xor %7, %7;"

        "movq 88(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"       //A11*B0
        
        "movq 0(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A0*B11
       
        "movq 80(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A10*B1
        
        "movq 8(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"            //A1*B10
        
        "movq 72(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"            //A9*B2
        
        "movq 16(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A2*B9
        
        "movq 64(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A8*B3
        
        "movq 24(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A3*B8
       
        "movq 56(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A7*B4
      
        "movq 32(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A4*B7
        
        "movq 48(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A6*B5
       
        "movq 40(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %8;"
        "movq %8, %4;"
        "adcq %%rdx, %%r10;"
        "adcq $0, %7;"        //A5*B6
                                //11
        "xor %8, %8;"
        
        "movq 96(%9), %%rax;"
        "mulq 0(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"    //A12*B0
        
        "movq 0(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"    //A0*B12
        
        
        "movq 88(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"    //A11*B1
        
        
        "movq 8(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A1*B11
        
        "movq 80(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A10*B2
        
        
        "movq 16(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A2*B10
        
        
        "movq 72(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A9*B3
        
        "movq 24(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A3*B9
       
        
        "movq 64(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A8*B4
        
        "movq 32(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A4*B8
        
        
        "movq 56(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A7*B5
        
        
        "movq 40(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %%r10;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A5*B7
        
        
        "movq 48(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %%r10;"
        "movq %%r10, %5;"
        "adcq %%rdx, %7;"
        "adcq $0, %8;"        //A6*B6
                                //12
        "xor %%r10, %%r10;"
        
        "movq 96(%9), %%rax;"
        "mulq 8(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A12*B1
       
        "movq 8(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A1*B12
       
        "movq 88(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A11*B2
        
        "movq 16(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A2*B11
        
        "movq 80(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A10*B3
       
        "movq 24(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A3*B10
       
        "movq 72(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A9*B4
        
        "movq 32(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A4*B9
        
        "movq 64(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A8*B5
        
        "movq 40(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A5*B8
        
        "movq 56(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;"   //A7*B6
       
        "movq 48(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %7;"
        "movq %7, %6;"
        "adcq %%rdx, %8;"
        "adcq $0, %%r10;" //A6*B7
                        //13                  break
                 
     //   "movq %8, %7;"
     //   "movq %%r10, %8;"
         
     //    "movq %8, %7;"
         "movq %%r10, %7;"
                    
        :  "=&g" ((res->_mp_d)[7]), "=&g" ((res->_mp_d)[8]), "=&g" ((res->_mp_d)[9]), "=&g" ((res->_mp_d)[10]), "=&g" ((res->_mp_d)[11]), "=&g" ((res->_mp_d)[12]), "=&g" ((res->_mp_d)[13]), "+r" (ra), "+r" (rb)//,  "=&X" ((res->_mp_d)[2]), "=&X" ((res->_mp_d)[3]), "=&X" ((res->_mp_d)[4]),  "=&X" ((res->_mp_d)[5]), "=&X" ((res->_mp_d)[6]), "=&X" ((res->_mp_d)[7]),  "=&X" ((res->_mp_d)[8]), "=&X" ((res->_mp_d)[9]), "=&X" ((res->_mp_d)[10]), "=&X" ((res->_mp_d)[11]), "=&X" ((res->_mp_d)[12]), "=&X" ((res->_mp_d)[13]),  "=&X" ((res->_mp_d)[14]), "=&X" ((res->_mp_d)[15]), "=&X" ((res->_mp_d)[16]),  "=&X" ((res->_mp_d)[17]), "=&X" ((res->_mp_d)[18]), "=&X" ((res->_mp_d)[19]),  "=&X" ((res->_mp_d)[20]), "=&X" ((res->_mp_d)[21]), "=&X" ((res->_mp_d)[22]), "=&X" ((res->_mp_d)[23]), "=&X" ((res->_mp_d)[24]), "=&X" ((res->_mp_d)[25])
        : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))//, "X" ((ain->_mp_d)[1]), "X" ((bin->_mp_d)[1]), "X" ((ain->_mp_d)[2]), "X" ((bin->_mp_d)[2]), "X" ((ain->_mp_d)[3]), "X" ((bin->_mp_d)[3]), "X" ((ain->_mp_d)[4]), "X" ((bin->_mp_d)[4]), "X" ((ain->_mp_d)[5]), "X" ((bin->_mp_d)[5]), "X" ((ain->_mp_d)[6]), "X" ((bin->_mp_d)[6]), "X" ((ain->_mp_d)[7]), "X" ((bin->_mp_d)[7]), "X" ((ain->_mp_d)[8]), "X" ((bin->_mp_d)[8]), "X" ((ain->_mp_d)[9]), "X" ((bin->_mp_d)[9]), "X" ((ain->_mp_d)[10]), "X" ((bin->_mp_d)[10]), "X" ((ain->_mp_d)[11]), "X" ((bin->_mp_d)[11]), "X" ((ain->_mp_d)[12]), "X" ((bin->_mp_d)[12])
        : "%rax", "%rdx", "%r8"
        );
    
        asm(
          
     //   "movq %8, %%r9;"
     //   "movq %7, %%r10;"             
        
        "xor %%r8, %%r8;"
        
        "movq 96(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B2
        
        "movq 16(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A2*B12
      
        "movq 88(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B3
        
        "movq 24(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A3*B11
        
        "movq 80(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B4
        
        "movq 32(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A4*B10
        
        "movq 72(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B5
       
        "movq 40(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A5*B9
        
        "movq 64(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B6
        
        "movq 48(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A6*B8
       
        "movq 56(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %8;"
        "movq %8, %0;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A7*B7
                            //14 
        "xor %8, %8;"
        
        "movq 96(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A12*B3
        
        "movq 24(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A3*B12
       
        "movq 88(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A11*B4
        
        "movq 32(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A4*B11
       
        "movq 80(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A10*B5
        
        "movq 40(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A5*B10
        
        "movq 72(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A9*B6
        
        "movq 48(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A6*B9
        
        "movq 64(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A8*B7
       
        "movq 56(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %7;"
        "movq %7, %1;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;" //A7*B8
                        //15 
        "xor %7, %7;"
        
        "movq 96(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A12*B4
        
        "movq 32(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A4*B12
        
        "movq 88(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A11*B5
       
        "movq 40(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A5*B11
        
        "movq 80(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A10*B6
        
        "movq 48(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A6*B10
     
        "movq 72(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A9*B7
        
        "movq 56(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A7*B9
       
        "movq 64(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %%r8;"
        "movq %%r8, %2;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A8*B8
                            //16
        "xor %%r8, %%r8;"

        "movq 96(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B5
        
        "movq 40(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A5*B12
        
        "movq 88(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B6
        
        "movq 48(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A6*B11
        
        "movq 80(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B7
        
        "movq 56(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A7*B10
        
        "movq 72(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B8
        
        "movq 64(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %8;"
        "movq %8, %3;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B9
                            //17 
        "xor %8, %8;"
        
        "movq 96(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A12*B6
        
        "movq 48(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A6*B12
        
        "movq 88(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A11*B7
       
        "movq 56(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A7*B11
       
        "movq 80(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A10*B8
        
        "movq 64(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A8*B10
       
        "movq 72(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %7;"
        "movq %7, %4;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A9*B9
                            //18
        "xor %7, %7;"
        
        "movq 96(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A12*B7
        
        "movq 56(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A7*B12
        
        "movq 88(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A11*B8
        
        "movq 64(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A8*B11
       
        "movq 80(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A10*B9
        
        "movq 72(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %%r8;"
        "movq %%r8, %5;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A9*B10
                            //19
        "xor %%r8, %%r8;"
        

        "movq 96(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B8
        
        "movq 64(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B12
        
        "movq 88(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B9
        
        "movq 72(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B11
        
        "movq 80(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %8;"
        "movq %8, %6;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B10
                            //20           break
                     
       // "movq %7, %8;"
        "movq %%r8, %8;"
                    
        :  "=&g" ((res->_mp_d)[14]), "=&g" ((res->_mp_d)[15]), "=&g" ((res->_mp_d)[16]), "=&g" ((res->_mp_d)[17]), "=&g" ((res->_mp_d)[18]), "=&g" ((res->_mp_d)[19]), "=&g" ((res->_mp_d)[20]), "+r" (ra), "+r" (rb)//,  "=&X" ((res->_mp_d)[2]), "=&X" ((res->_mp_d)[3]), "=&X" ((res->_mp_d)[4]),  "=&X" ((res->_mp_d)[5]), "=&X" ((res->_mp_d)[6]), "=&X" ((res->_mp_d)[7]),  "=&X" ((res->_mp_d)[8]), "=&X" ((res->_mp_d)[9]), "=&X" ((res->_mp_d)[10]), "=&X" ((res->_mp_d)[11]), "=&X" ((res->_mp_d)[12]), "=&X" ((res->_mp_d)[13]),  "=&X" ((res->_mp_d)[14]), "=&X" ((res->_mp_d)[15]), "=&X" ((res->_mp_d)[16]),  "=&X" ((res->_mp_d)[17]), "=&X" ((res->_mp_d)[18]), "=&X" ((res->_mp_d)[19]),  "=&X" ((res->_mp_d)[20]), "=&X" ((res->_mp_d)[21]), "=&X" ((res->_mp_d)[22]), "=&X" ((res->_mp_d)[23]), "=&X" ((res->_mp_d)[24]), "=&X" ((res->_mp_d)[25])
        : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))//, "X" ((ain->_mp_d)[1]), "X" ((bin->_mp_d)[1]), "X" ((ain->_mp_d)[2]), "X" ((bin->_mp_d)[2]), "X" ((ain->_mp_d)[3]), "X" ((bin->_mp_d)[3]), "X" ((ain->_mp_d)[4]), "X" ((bin->_mp_d)[4]), "X" ((ain->_mp_d)[5]), "X" ((bin->_mp_d)[5]), "X" ((ain->_mp_d)[6]), "X" ((bin->_mp_d)[6]), "X" ((ain->_mp_d)[7]), "X" ((bin->_mp_d)[7]), "X" ((ain->_mp_d)[8]), "X" ((bin->_mp_d)[8]), "X" ((ain->_mp_d)[9]), "X" ((bin->_mp_d)[9]), "X" ((ain->_mp_d)[10]), "X" ((bin->_mp_d)[10]), "X" ((ain->_mp_d)[11]), "X" ((bin->_mp_d)[11]), "X" ((ain->_mp_d)[12]), "X" ((bin->_mp_d)[12])
        : "%rax", "%rdx", "%r8"
        );
    
    asm (
                 
     //   "movq %5, %%r8;"
     //   "movq %6, %%r10;"
    
                    
        "xor %%r9, %%r9;"


        "movq 96(%6), %%rax;"
        "mulq 72(%7);"            
        "addq %%rax, %5;"
        "adcq %%rdx, %4;"
        "adcq $0, %%r9;"    //A12*B9
        
        "movq 72(%6), %%rax;"
        "mulq 96(%7);"            
        "addq %%rax, %5;"
        "adcq %%rdx, %4;"
        "adcq $0, %%r9;"    //A9*B12
       
        "movq 88(%6), %%rax;"
        "mulq 80(%7);"            
        "addq %%rax, %5;"
        "adcq %%rdx, %4;"
        "adcq $0, %%r9;"    //A11*B10
        
        "movq 80(%6), %%rax;"
        "mulq 88(%7);"            
        "addq %%rax, %5;"
        "movq %5, %0;"
        "adcq %%rdx, %4;"
        "adcq $0, %%r9;"    //A10*B11
                            //21
        "xor %5, %5;"
        

        "movq 96(%6), %%rax;"
        "mulq 80(%7);"            
        "addq %%rax, %4;"
        "adcq %%rdx, %%r9;"
        "adcq $0, %5;"   //A12*B10
        
        "movq 80(%6), %%rax;"
        "mulq 96(%7);"            
        "addq %%rax, %4;"
        "adcq %%rdx, %%r9;"
        "adcq $0, %5;"   //A10*B12
        
        "movq 88(%6), %%rax;"
        "mulq 88(%7);"            
        "addq %%rax, %4;"
        "movq %4, %1;"
        "adcq %%rdx, %%r9;"
        "adcq $0, %5;"   //A11*B11
                            //22
     //  "xor %4, %4;"
        
        "movq 96(%6), %%rax;"
        "mulq 88(%7);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %5;"
     //   "adcq $0, %4;"    //A12*B11
       
        "movq 88(%6), %%rax;"
        "mulq 96(%7);"            
        "addq %%rax, %%r9;"
        "movq %%r9, %2;"
        "adcq %%rdx, %5;"
     //   "adcq $0, %4;"    //A11*B12
                            //23
     //   "xor %%r9, %%r9;"
        
        "movq 96(%6), %%rax;"
        "mulq 96(%7);"            
        "addq %%rax, %5;"
        "movq %5, %3;"          //A12*B12
  //      "adcq %%rdx, %4;"     //24
  //      "adcq $0, %4;"    
                            
       
   //     "movq %5, %4;" //25

        
        :  "=&g" ((res->_mp_d)[21]), "=&g" ((res->_mp_d)[22]), "=&g" ((res->_mp_d)[23]), "=&g" ((res->_mp_d)[24]), "+r" (rb), "+r" (ra)//,  "=&X" ((res->_mp_d)[2]), "=&X" ((res->_mp_d)[3]), "=&X" ((res->_mp_d)[4]),  "=&X" ((res->_mp_d)[5]), "=&X" ((res->_mp_d)[6]), "=&X" ((res->_mp_d)[7]),  "=&X" ((res->_mp_d)[8]), "=&X" ((res->_mp_d)[9]), "=&X" ((res->_mp_d)[10]), "=&X" ((res->_mp_d)[11]), "=&X" ((res->_mp_d)[12]), "=&X" ((res->_mp_d)[13]),  "=&X" ((res->_mp_d)[14]), "=&X" ((res->_mp_d)[15]), "=&X" ((res->_mp_d)[16]),  "=&X" ((res->_mp_d)[17]), "=&X" ((res->_mp_d)[18]), "=&X" ((res->_mp_d)[19]),  "=&X" ((res->_mp_d)[20]), "=&X" ((res->_mp_d)[21]), "=&X" ((res->_mp_d)[22]), "=&X" ((res->_mp_d)[23]), "=&X" ((res->_mp_d)[24]), "=&X" ((res->_mp_d)[25])
        : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))//, "r" (ra)//, "X" ((ain->_mp_d)[1]), "X" ((bin->_mp_d)[1]), "X" ((ain->_mp_d)[2]), "X" ((bin->_mp_d)[2]), "X" ((ain->_mp_d)[3]), "X" ((bin->_mp_d)[3]), "X" ((ain->_mp_d)[4]), "X" ((bin->_mp_d)[4]), "X" ((ain->_mp_d)[5]), "X" ((bin->_mp_d)[5]), "X" ((ain->_mp_d)[6]), "X" ((bin->_mp_d)[6]), "X" ((ain->_mp_d)[7]), "X" ((bin->_mp_d)[7]), "X" ((ain->_mp_d)[8]), "X" ((bin->_mp_d)[8]), "X" ((ain->_mp_d)[9]), "X" ((bin->_mp_d)[9]), "X" ((ain->_mp_d)[10]), "X" ((bin->_mp_d)[10]), "X" ((ain->_mp_d)[11]), "X" ((bin->_mp_d)[11]), "X" ((ain->_mp_d)[12]), "X" ((bin->_mp_d)[12])
        : "%rax", "%rdx", "%r9"
        );
    
   // printf("\nIN FUN res: %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (resp)[25], (resp)[24], (resp)[23], (resp)[22], (resp)[21], (resp)[20], (resp)[19], (resp)[18], (resp)[17], (resp)[16], (resp)[15], (resp)[14], (resp)[13], (resp)[12], (resp)[11], (resp)[10], (resp)[9], (resp)[8], (resp)[7], (resp)[6], (resp)[5], (resp)[4], (resp)[3], (resp)[2], (resp)[1], (resp)[0]);

    if( (res->_mp_d)[24]==0 )
        res->_mp_size = 24*negFlag;
    else
        res->_mp_size = 25*negFlag;

}

void mpz_mul_x86_2(mpz_t res, mpz_t ain, mpz_t bin){
    
    unsigned long int ra=0;
    unsigned long int rb=0;

    
  //  mpz_set(tmp1, ain);
    //mpz_set(tmp2, bin);
    
    switch (ain->_mp_size) {
            
        case 0:
            res->_mp_size=0;
            return;
            
   /*     case 1:
            mpz_set(res,bin);
            return;
     */       
        case 12:
          
        //    if( ain->_mp_alloc<13)
        //        _mpz_realloc(ain,13);

         
            (ain->_mp_d)[12]=0;
            
           // flag++;
            break;
    }
    
    
    switch (bin->_mp_size) {
            
        case 0:
            res->_mp_size=0;
            return;
            
        case 1:
            mpz_set(res,ain);
            return;
            
        case 12:

      //      if( bin->_mp_alloc<13)
      //           _mpz_realloc(bin,13);
                
            (bin->_mp_d)[12]=0;
            
           // flag++;
        //    if( res->_mp_alloc<25)
        //        _mpz_realloc(res,25);
            
            break;
            
    /*    case 13:
            
     //       if( res->_mp_alloc<25)
     //           _mpz_realloc(res,25);
            
            break;
      */      

    }
    
       
    asm (            
         "movq 0(%9), %%rax;" 
         "mulq 0(%10);"
         "movq %%rax, %0;"    
         "movq %%rdx, %7;"           //A0*B0
         //0
         
         "xorq %%r10, %%r10;" 
         
         "movq 8(%9), %%rax;"      //problem line
         "mulq 0(%10);"              
         "addq %%rax, %7;"     
         "movq %%rdx, %8;"  
         "adcq $0, %8;"              //A1*B0
         
         "movq 0(%9), %%rax;"  
         "mulq 8(%10);"         
         "addq %%rax, %7;" 
         "movq %7, %1;"  
         "adcq %%rdx,%8;"    
         "adcq $0, %%r10;"                //A0*B1
         //1
         
         "xorq %7, %7;" 
         
         "movq 0(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"           //A0*B2
         
         "movq 8(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A1*B1
         
         "movq 16(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %8;"    
         "movq %8, %2;" 
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A2*B0
         //2
         "xorq %8, %8;"  
         
         "movq 24(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"              //A3*B0
         
         "movq 0(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"            //A0*B3
         
         "movq 16(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A2*B1
         
         "movq 8(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %%r10;"   
         "movq %%r10, %3;" 
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A1*B2
         //3
         "xorq %%r10, %%r10;" 
         
         "movq 32(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"            //A4*B0
         
         "movq 0(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"            //A0*B4
         
         "movq 24(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A3*B1
         
         "movq 8(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A1*B3
         
         "movq 16(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %7;"  
         "movq %7, %4;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A2*B2
         //4
         "xor %7, %7;"
         
         "movq 40(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A5*B0
         
         "movq 0(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A0*B5
         
         "movq 32(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A4*B1
         
         "movq 8(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A1*B4
         
         "movq 24(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A3*B2
         
         "movq 16(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %8;"
         "movq %8, %5;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A2*B3
         //5
         "xor %8, %8;"
         
         "movq 48(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A6*B0
         
         "movq 0(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A0*B6
         
         "movq 40(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A5*B1
         
         "movq 8(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A1*B5
         
         "movq 32(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"            //A4*B2
         
         "movq 16(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"            //A2*B4
         
         "movq 24(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %%r10;"
         "movq %%r10, %6;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"  
         
         //    "movq %%r8, %7;" 
         //   "movq %%r9, %8;" 
         
         //A3*B3
         //6       //break
         :  "=&g" ((res->_mp_d)[0]), "=&g" ((res->_mp_d)[1]), "=&g" ((res->_mp_d)[2]), "=&g" ((res->_mp_d)[3]), "=&g" ((res->_mp_d)[4]), "=&g" ((res->_mp_d)[5]), "=&g" ((res->_mp_d)[6]), "+r" (ra), "+r" (rb)
         : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))
         : "%rax", "%rdx", "%r10"
         );
    
    asm (
         
         //   "movq %11, %%r8;" 
         //   "movq %12, %%r9;" 
         
         
         //"xor %%r8, %%r8;"
         //"xor %%r9, %%r9;"
         "xor %%r10, %%r10;"
         
         
         "movq 56(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A7*B0
         
         "movq 0(%9), %%rax;"
         "mulq 56(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A0*B7
         
         "movq 48(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A6*B1
         
         "movq 8(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A1*B6
         
         "movq 40(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A5*B2
         
         "movq 16(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"        //A2*B5
         
         "movq 32(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A4*B3
         
         "movq 24(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %7;"  
         "movq %7, %0;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A3*B4
         //7
         "xor %7, %7;"
         
         "movq 64(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"       //A8*B0
         
         "movq 0(%9), %%rax;"
         "mulq 64(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"             //A0*B8
         
         "movq 56(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A7*B1
         
         "movq 8(%9), %%rax;"
         "mulq 56(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A1*B7
         
         "movq 48(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A6*B2
         
         "movq 16(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A2*B6
         
         "movq 40(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A5*B3
         
         "movq 24(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A3*B5
         
         "movq 32(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %8;" 
         "movq %8, %1;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A4*B4
         //8 
         
         "xor %8, %8;"
         
         "movq 72(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"         //A9*B0
         
         "movq 0(%9), %%rax;"
         "mulq 72(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A0*B9
         
         "movq 64(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"          //A8*B1
         
         "movq 8(%9), %%rax;"
         "mulq 64(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"         //A1*B8
         
         "movq 56(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"  //A7*B2
         
         "movq 16(%9), %%rax;"
         "mulq 56(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"         //A2*B7
         
         "movq 48(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A6*B3
         
         "movq 24(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A3*B6
         
         "movq 40(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A5*B4
         
         "movq 32(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %%r10;"
         "movq %%r10, %2;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A4*B5
         //9  
         "xor %%r10, %%r10;"
         
         "movq 80(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A10*B0
         
         "movq 0(%9), %%rax;"
         "mulq 80(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"        //A0*B10
         
         "movq 72(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A9*B1
         
         "movq 8(%9), %%rax;"
         "mulq 72(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A1*B9
         
         "movq 64(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A8*B2
         
         "movq 16(%9), %%rax;"
         "mulq 64(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A2*B8
         
         "movq 56(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"          //A7*B3
         
         "movq 24(%9), %%rax;"
         "mulq 56(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A3*B7
         
         "movq 48(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A6*B4
         
         "movq 32(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A4*B6
         
         "movq 40(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %7;"
         "movq %7, %3;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A5*B5
         //10
         "xor %7, %7;"
         
         "movq 88(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"       //A11*B0
         
         "movq 0(%9), %%rax;"
         "mulq 88(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A0*B11
         
         "movq 80(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A10*B1
         
         "movq 8(%9), %%rax;"
         "mulq 80(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A1*B10
         
         "movq 72(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A9*B2
         
         "movq 16(%9), %%rax;"
         "mulq 72(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A2*B9
         
         "movq 64(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A8*B3
         
         "movq 24(%9), %%rax;"
         "mulq 64(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A3*B8
         
         "movq 56(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A7*B4
         
         "movq 32(%9), %%rax;"
         "mulq 56(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A4*B7
         
         "movq 48(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A6*B5
         
         "movq 40(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %8;"
         "movq %8, %4;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A5*B6
         //11
         "xor %8, %8;"
         
         "movq 96(%9), %%rax;"
         "mulq 0(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"    //A12*B0
         
         "movq 0(%9), %%rax;"
         "mulq 96(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"    //A0*B12
         
         
         "movq 88(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"    //A11*B1
         
         
         "movq 8(%9), %%rax;"
         "mulq 88(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A1*B11
         
         "movq 80(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A10*B2
         
         
         "movq 16(%9), %%rax;"
         "mulq 80(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A2*B10
         
         
         "movq 72(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A9*B3
         
         "movq 24(%9), %%rax;"
         "mulq 72(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A3*B9
         
         
         "movq 64(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A8*B4
         
         "movq 32(%9), %%rax;"
         "mulq 64(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A4*B8
         
         
         "movq 56(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A7*B5
         
         
         "movq 40(%9), %%rax;"
         "mulq 56(%10);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A5*B7
         
         
         "movq 48(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %%r10;"
         "movq %%r10, %5;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A6*B6
         //12
         "xor %%r10, %%r10;"
         
         "movq 96(%9), %%rax;"
         "mulq 8(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A12*B1
         
         "movq 8(%9), %%rax;"
         "mulq 96(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A1*B12
         
         "movq 88(%9), %%rax;"
         "mulq 16(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A11*B2
         
         "movq 16(%9), %%rax;"
         "mulq 88(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A2*B11
         
         "movq 80(%9), %%rax;"
         "mulq 24(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A10*B3
         
         "movq 24(%9), %%rax;"
         "mulq 80(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A3*B10
         
         "movq 72(%9), %%rax;"
         "mulq 32(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A9*B4
         
         "movq 32(%9), %%rax;"
         "mulq 72(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A4*B9
         
         "movq 64(%9), %%rax;"
         "mulq 40(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A8*B5
         
         "movq 40(%9), %%rax;"
         "mulq 64(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A5*B8
         
         "movq 56(%9), %%rax;"
         "mulq 48(%10);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A7*B6
         
         "movq 48(%9), %%rax;"
         "mulq 56(%10);"            
         "addq %%rax, %7;"
         "movq %7, %6;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;" //A6*B7
         //13                  break
         
         //   "movq %8, %7;"
         //   "movq %%r10, %8;"
         
         //    "movq %8, %7;"
         "movq %%r10, %7;"
         
         :  "=&g" ((res->_mp_d)[7]), "=&g" ((res->_mp_d)[8]), "=&g" ((res->_mp_d)[9]), "=&g" ((res->_mp_d)[10]), "=&g" ((res->_mp_d)[11]), "=&g" ((res->_mp_d)[12]), "=&g" ((res->_mp_d)[13]), "+r" (ra), "+r" (rb)
         : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))
         : "%rax", "%rdx", "%r8"
         );
    
    asm(
        
        //   "movq %8, %%r9;"
        //   "movq %7, %%r10;"             
        
        "xor %%r8, %%r8;"
        
        "movq 96(%9), %%rax;"
        "mulq 16(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B2
        
        "movq 16(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A2*B12
        
        "movq 88(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B3
        
        "movq 24(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A3*B11
        
        "movq 80(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B4
        
        "movq 32(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A4*B10
        
        "movq 72(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B5
        
        "movq 40(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A5*B9
        
        "movq 64(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B6
        
        "movq 48(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A6*B8
        
        "movq 56(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %8;"
        "movq %8, %0;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A7*B7
        //14 
        "xor %8, %8;"
        
        "movq 96(%9), %%rax;"
        "mulq 24(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A12*B3
        
        "movq 24(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A3*B12
        
        "movq 88(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A11*B4
        
        "movq 32(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A4*B11
        
        "movq 80(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A10*B5
        
        "movq 40(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A5*B10
        
        "movq 72(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A9*B6
        
        "movq 48(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A6*B9
        
        "movq 64(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A8*B7
        
        "movq 56(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %7;"
        "movq %7, %1;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;" //A7*B8
        //15 
        "xor %7, %7;"
        
        "movq 96(%9), %%rax;"
        "mulq 32(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A12*B4
        
        "movq 32(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A4*B12
        
        "movq 88(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A11*B5
        
        "movq 40(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A5*B11
        
        "movq 80(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A10*B6
        
        "movq 48(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A6*B10
        
        "movq 72(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A9*B7
        
        "movq 56(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A7*B9
        
        "movq 64(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %%r8;"
        "movq %%r8, %2;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A8*B8
        //16
        "xor %%r8, %%r8;"
        
        "movq 96(%9), %%rax;"
        "mulq 40(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B5
        
        "movq 40(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A5*B12
        
        "movq 88(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B6
        
        "movq 48(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A6*B11
        
        "movq 80(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B7
        
        "movq 56(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A7*B10
        
        "movq 72(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B8
        
        "movq 64(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %8;"
        "movq %8, %3;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B9
        //17 
        "xor %8, %8;"
        
        "movq 96(%9), %%rax;"
        "mulq 48(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A12*B6
        
        "movq 48(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A6*B12
        
        "movq 88(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A11*B7
        
        "movq 56(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A7*B11
        
        "movq 80(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A10*B8
        
        "movq 64(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A8*B10
        
        "movq 72(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %7;"
        "movq %7, %4;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A9*B9
        //18
        "xor %7, %7;"
        
        "movq 96(%9), %%rax;"
        "mulq 56(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A12*B7
        
        "movq 56(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A7*B12
        
        "movq 88(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A11*B8
        
        "movq 64(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A8*B11
        
        "movq 80(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A10*B9
        
        "movq 72(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %%r8;"
        "movq %%r8, %5;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A9*B10
        //19
        "xor %%r8, %%r8;"
        
        
        "movq 96(%9), %%rax;"
        "mulq 64(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B8
        
        "movq 64(%9), %%rax;"
        "mulq 96(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B12
        
        "movq 88(%9), %%rax;"
        "mulq 72(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B9
        
        "movq 72(%9), %%rax;"
        "mulq 88(%10);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B11
        
        "movq 80(%9), %%rax;"
        "mulq 80(%10);"            
        "addq %%rax, %8;"
        "movq %8, %6;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B10
        //20           break
        
        // "movq %7, %8;"
        "movq %%r8, %8;"
        
        :  "=&g" ((res->_mp_d)[14]), "=&g" ((res->_mp_d)[15]), "=&g" ((res->_mp_d)[16]), "=&g" ((res->_mp_d)[17]), "=&g" ((res->_mp_d)[18]), "=&g" ((res->_mp_d)[19]), "=&g" ((res->_mp_d)[20]), "+r" (ra), "+r" (rb)
        : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))
        : "%rax", "%rdx", "%r8"
        );
    
    asm (
         
         //   "movq %5, %%r8;"
         //   "movq %6, %%r10;"
         
         
         "xor %%r9, %%r9;"
         
         
         "movq 96(%6), %%rax;"
         "mulq 72(%7);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A12*B9
         
         "movq 72(%6), %%rax;"
         "mulq 96(%7);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A9*B12
         
         "movq 88(%6), %%rax;"
         "mulq 80(%7);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A11*B10
         
         "movq 80(%6), %%rax;"
         "mulq 88(%7);"            
         "addq %%rax, %5;"
         "movq %5, %0;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A10*B11
         //21
         "xor %5, %5;"
         
         
         "movq 96(%6), %%rax;"
         "mulq 80(%7);"            
         "addq %%rax, %4;"
         "adcq %%rdx, %%r9;"
         "adcq $0, %5;"   //A12*B10
         
         "movq 80(%6), %%rax;"
         "mulq 96(%7);"            
         "addq %%rax, %4;"
         "adcq %%rdx, %%r9;"
         "adcq $0, %5;"   //A10*B12
         
         "movq 88(%6), %%rax;"
         "mulq 88(%7);"            
         "addq %%rax, %4;"
         "movq %4, %1;"
         "adcq %%rdx, %%r9;"
         "adcq $0, %5;"   //A11*B11
         //22
         //  "xor %4, %4;"
         
         "movq 96(%6), %%rax;"
         "mulq 88(%7);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         //   "adcq $0, %4;"    //A12*B11
         
         "movq 88(%6), %%rax;"
         "mulq 96(%7);"            
         "addq %%rax, %%r9;"
         "movq %%r9, %2;"
         "adcq %%rdx, %5;"
         //   "adcq $0, %4;"    //A11*B12
         //23
         //   "xor %%r9, %%r9;"
         
         "movq 96(%6), %%rax;"
         "mulq 96(%7);"            
         "addq %%rax, %5;"
         "movq %5, %3;"          //A12*B12
         //      "adcq %%rdx, %4;"     //24
         //      "adcq $0, %4;"    
         
         
         //     "movq %5, %4;" //25
         
         
         :  "=&g" ((res->_mp_d)[21]), "=&g" ((res->_mp_d)[22]), "=&g" ((res->_mp_d)[23]), "=&g" ((res->_mp_d)[24]), "+r" (rb), "+r" (ra)
         : "r" ((ain->_mp_d)), "r" ((bin->_mp_d))
         : "%rax", "%rdx", "%r9"
         );
    
    
    if( (res->_mp_d)[24]==0)
        res->_mp_size = 24;
    else
        res->_mp_size = 25;
    
}

//special purpose multiplication for Barrett reduction. ain is the fixed 13 limb number and bin is the 24 or 25 limb number
void mpz_mul_x86_3(mpz_t res, mpz_t bin){
    
    unsigned long int ra=0;
    unsigned long int rb=0;
    
  //  printf("\nMUL ENTERED\n");
    
    
    switch (bin->_mp_size) {
        
        case 0:
            mpz_set_ui(res,0);
            return;
            
        case 24:
             (bin->_mp_d)[24]=0;
            break;
        
        case 13:
            (bin->_mp_d)[13]=0;
            (bin->_mp_d)[14]=0;
            (bin->_mp_d)[15]=0;
            (bin->_mp_d)[16]=0;
            (bin->_mp_d)[17]=0;
            (bin->_mp_d)[18]=0;
            (bin->_mp_d)[19]=0;
            (bin->_mp_d)[20]=0;
            (bin->_mp_d)[21]=0;
            (bin->_mp_d)[22]=0;
            (bin->_mp_d)[23]=0;
            (bin->_mp_d)[24]=0;
            break;
            
        case 12:
            (bin->_mp_d)[12]=0;
            (bin->_mp_d)[13]=0;
            (bin->_mp_d)[14]=0;
            (bin->_mp_d)[15]=0;
            (bin->_mp_d)[16]=0;
            (bin->_mp_d)[17]=0;
            (bin->_mp_d)[18]=0;
            (bin->_mp_d)[19]=0;
            (bin->_mp_d)[20]=0;
            (bin->_mp_d)[21]=0;
            (bin->_mp_d)[22]=0;
            (bin->_mp_d)[23]=0;
            (bin->_mp_d)[24]=0;
            break;
            
    }
    
    
  
/*    
 a
 d88b672f67f5f8ff
 5ef4b87df2399756
 d528570ce8895ea5
 a177d0ac7b179728
 599b9401c7d4aa88
 3c3b41cfa08954ce
 3c8f093f6421c1ec
 f8e9a6dddd32146f
 8affc7554c74326c
 9314b528ab0b9bec
 fde79c768b4a3230
 b0300f1dbecda98c
    */
    
      asm (
         
         "xor %4, %4;"
    //     
         "movq $0xa, %%rax;"
         "mulq 88(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
            "adcq $0, %4;"    //A12*B11
         
         "movq $0xd88b672f67f5f8ff, %%rax;"
         "mulq 96(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
          "adcq $0, %4;"    //A11*B12
                            
         "movq $0xb0300f1dbecda98c, %%rax;"
         "mulq 184(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                 //A0*B23  n
         
         "movq $0xfde79c768b4a3230, %%rax;"
         "mulq 176(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A1*B22  n
         
         "movq $0x9314b528ab0b9bec, %%rax;"
         "mulq 168(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A2*B21  n
         
         "movq $0x8affc7554c74326c, %%rax;"
         "mulq 160(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A3*B20  n
         
         "movq $0xf8e9a6dddd32146f, %%rax;"
         "mulq 152(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A4*B19  n
   //      
         "movq $0x3c8f093f6421c1ec, %%rax;"
         "mulq 144(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A5*B18  n
         
         "movq $0x3c3b41cfa08954ce, %%rax;"
         "mulq 136(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A6*B17  n
         
         "movq $0x599b9401c7d4aa88, %%rax;"
         "mulq 128(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A7*B16  n
         
         "movq $0xa177d0ac7b179728, %%rax;"
         "mulq 120(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A8*B15  n
         
         "movq $0xd528570ce8895ea5, %%rax;"
         "mulq 112(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A9*B14  n 
         
         "movq $0x5ef4b87df2399756, %%rax;"
         "mulq 104(%6);"            
         "addq %%rax, %%r9;"
       //  "movq %%r9, %2;"
         "adcq %%rdx, %5;"
         "adcq $0, %4;"                   //A10*B13  n
    //                            //23
         "xor %%r9, %%r9;"
         
         "movq $0xa, %%rax;"
         "mulq 96(%6);"            
         "addq %%rax, %5;"         //A12*B12
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"         
                            
         "movq $0xb0300f1dbecda98c, %%rax;"
         "mulq 192(%6);"            
         "addq %%rax, %5;"       
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A0*B24  n
         
         "movq $0xfde79c768b4a3230, %%rax;"
         "mulq 184(%6);"            
         "addq %%rax, %5;"       
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A1*B23  n
         
         "movq $0x9314b528ab0b9bec, %%rax;"
         "mulq 176(%6);"            
         "addq %%rax, %5;"       
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A2*B22  n
         
         "movq $0x8affc7554c74326c, %%rax;"
         "mulq 168(%6);"            
         "addq %%rax, %5;"        
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A3*B21 n
    //     
         "movq $0xf8e9a6dddd32146f, %%rax;"
         "mulq 160(%6);"            
         "addq %%rax, %5;"       
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A4*B20  n
         
         "movq $0x3c8f093f6421c1ec, %%rax;"
         "mulq 152(%6);"            
         "addq %%rax, %5;"        
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A5*B19  n
         
         "movq $0x3c3b41cfa08954ce, %%rax;"
         "mulq 144(%6);"            
         "addq %%rax, %5;"       
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A6*B18  n
         
         "movq $0x599b9401c7d4aa88, %%rax;"
         "mulq 136(%6);"            
         "addq %%rax, %5;"     
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A7*B17  n
         
         "movq $0xa177d0ac7b179728, %%rax;"
         "mulq 128(%6);"            
         "addq %%rax, %5;"        
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A8*B16  n
         
         "movq $0xd528570ce8895ea5, %%rax;"
         "mulq 120(%6);"            
         "addq %%rax, %5;"       
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A9*B15  n
         
         "movq $0x5ef4b87df2399756, %%rax;"
         "mulq 112(%6);"            
         "addq %%rax, %5;"    
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A10*B14  n
  //       
         "movq $0xd88b672f67f5f8ff, %%rax;"
         "mulq 104(%6);"            
         "addq %%rax, %5;"
         "movq %5, %0;"         
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                   //A11*B13  n
                                //24
          "xor %5, %5;"
         
         "movq $0xfde79c768b4a3230, %%rax;"
         "mulq 192(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                 //A1*B24  n
         
         "movq $0x9314b528ab0b9bec, %%rax;"
         "mulq 184(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"             //A2*B23  n
         
         "movq $0x8affc7554c74326c, %%rax;"
         "mulq 176(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A3*B22  n
         
         "movq $0xf8e9a6dddd32146f, %%rax;"
         "mulq 168(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A4*B21  n
         
         "movq $0x3c8f093f6421c1ec, %%rax;"
         "mulq 160(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A5*B20  n
         
         "movq $0x3c3b41cfa08954ce, %%rax;"
         "mulq 152(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A6*B19  n
   //      
         "movq $0x599b9401c7d4aa88, %%rax;"
         "mulq 144(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A7*B18  n
         
         "movq $0xa177d0ac7b179728, %%rax;"
         "mulq 136(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A8*B17  n
         
         "movq $0xd528570ce8895ea5, %%rax;"
         "mulq 128(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A9*B16  n
         
         "movq $0x5ef4b87df2399756, %%rax;"
         "mulq 120(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A10*B15  n
         
         "movq $0xd88b672f67f5f8ff, %%rax;"
         "mulq 112(%6);"            
         "addq %%rax, %4;"     
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A11*B14  n
         
         "movq $0xa, %%rax;"
         "mulq 104(%6);"            
         "addq %%rax, %4;"
         "movq %4, %1;" 
         "adcq %%rdx, %%r9;"     
         "adcq $0, %5;"                   //A12*B13  n
    //                        //25
                            
         "xor %4, %4;"
         
         "movq $0x9314b528ab0b9bec, %%rax;"
         "mulq 192(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                 //A2*B24  n
         
         "movq $0x8affc7554c74326c, %%rax;"
         "mulq 184(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A3*B23  n
         
         "movq $0xf8e9a6dddd32146f, %%rax;"
         "mulq 176(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A4*B22  n
         
         "movq $0x3c8f093f6421c1ec, %%rax;"
         "mulq 168(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A5*B21  n
         
         "movq $0x3c3b41cfa08954ce, %%rax;"
         "mulq 160(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A6*B20  n
         
         "movq $0x599b9401c7d4aa88, %%rax;"
         "mulq 152(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A7*B19  n
   //      
         "movq $0xa177d0ac7b179728, %%rax;"
         "mulq 144(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A8*B18  n
         
         "movq $0xd528570ce8895ea5, %%rax;"
         "mulq 136(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A9*B17  n 
         
         "movq $0x5ef4b87df2399756, %%rax;"
         "mulq 128(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A10*B16 n
         
         "movq $0xd88b672f67f5f8ff, %%rax;"
         "mulq 120(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A11*B15  n 
         
         "movq $0xa, %%rax;"
         "mulq 112(%6);"            
         "addq %%rax, %%r9;"
         "movq %%r9, %2;"
         "adcq %%rdx, %5;"     
         "adcq $0, %4;"                   //A12*B14  n
                            //26
         "xor %%r9, %%r9;"
         
         
         "movq $0x8affc7554c74326c, %%rax;"
         "mulq 192(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A3*B24  n
         
       //  
         "movq $0xf8e9a6dddd32146f, %%rax;"
         "mulq 184(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A4*B23  n
         
         
         "movq $0x3c8f093f6421c1ec, %%rax;"
         "mulq 176(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A5*B22  n
         
         
         "movq $0x3c3b41cfa08954ce, %%rax;"
         "mulq 168(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A6*B21  n
         
         
         "movq $0x599b9401c7d4aa88, %%rax;"
         "mulq 160(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A7*B20  n
         
     //    
         "movq $0xa177d0ac7b179728, %%rax;"
         "mulq 152(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A8*B19  n
         
         
         "movq $0xd528570ce8895ea5, %%rax;"
         "mulq 144(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A9*B18  n
         
         
         "movq $0x5ef4b87df2399756, %%rax;"
         "mulq 136(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A10*B17  n
         
         
         "movq $0xd88b672f67f5f8ff, %%rax;"
         "mulq 128(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A11*B16  n
         
       
         "movq $0xa, %%rax;"
         "mulq 120(%6);"            
         "addq %%rax, %5;"
         "movq %5, %3;"
         "adcq %%rdx, %4;"     
         "adcq $0, %%r9;"                    //A12*B15  n
                                            //27
         
         
         "movq %%r9, %5;"
         
         
         :   "=&g" ((res->_mp_d)[24]), "=&g" ((res->_mp_d)[25]), "=&g" ((res->_mp_d)[26]), "=&g" ((res->_mp_d)[27]), "+r" (rb), "+r" (ra)//,  "=&X" ((res->_mp_d)[2]), "=&X" ((res->_mp_d)[3]), "=&X" ((res->_mp_d)[4]),  "=&X" ((res->_mp_d)[5]), "=&X" ((res->_mp_d)[6]), "=&X" ((res->_mp_d)[7]),  "=&X" ((res->_mp_d)[8]), "=&X" ((res->_mp_d)[9]), "=&X" ((res->_mp_d)[10]), "=&X" ((res->_mp_d)[11]), "=&X" ((res->_mp_d)[12]), "=&X" ((res->_mp_d)[13]),  "=&X" ((res->_mp_d)[14]), "=&X" ((res->_mp_d)[15]), "=&X" ((res->_mp_d)[16]),  "=&X" ((res->_mp_d)[17]), "=&X" ((res->_mp_d)[18]), "=&X" ((res->_mp_d)[19]),  "=&X" ((res->_mp_d)[20]), "=&X" ((res->_mp_d)[21]), "=&X" ((res->_mp_d)[22]), "=&X" ((res->_mp_d)[23]), "=&X" ((res->_mp_d)[24]), "=&X" ((res->_mp_d)[25])
         :  "r" ((bin->_mp_d))//, "r" (ra)//, "X" ((ain->_mp_d)[1]), "X" ((bin->_mp_d)[1]), "X" ((ain->_mp_d)[2]), "X" ((bin->_mp_d)[2]), "X" ((ain->_mp_d)[3]), "X" ((bin->_mp_d)[3]), "X" ((ain->_mp_d)[4]), "X" ((bin->_mp_d)[4]), "X" ((ain->_mp_d)[5]), "X" ((bin->_mp_d)[5]), "X" ((ain->_mp_d)[6]), "X" ((bin->_mp_d)[6]), "X" ((ain->_mp_d)[7]), "X" ((bin->_mp_d)[7]), "X" ((ain->_mp_d)[8]), "X" ((bin->_mp_d)[8]), "X" ((ain->_mp_d)[9]), "X" ((bin->_mp_d)[9]), "X" ((ain->_mp_d)[10]), "X" ((bin->_mp_d)[10]), "X" ((ain->_mp_d)[11]), "X" ((bin->_mp_d)[11]), "X" ((ain->_mp_d)[12]), "X" ((bin->_mp_d)[12])
         : "%rax", "%rdx", "%r9"
         );
    
    
    /*    
     a
     
     d88b672f67f5f8ff
     5ef4b87df2399756
     d528570ce8895ea5
     a177d0ac7b179728
     599b9401c7d4aa88
     3c3b41cfa08954ce
     3c8f093f6421c1ec
     f8e9a6dddd32146f
     8affc7554c74326c
     9314b528ab0b9bec
     fde79c768b4a3230
     b0300f1dbecda98c
     */
    
    asm(
        "xor %%r9, %%r9;" // 10 is ra, 11 is rb
        
        "movq $0xf8e9a6dddd32146f, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A4*B24  n
                               
        "movq $0x3c8f093f6421c1ec, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"        //A5*B23  n 
        
        "movq $0x3c3b41cfa08954ce, %%rax;"
        "mulq 176(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A6*B22  n
        
        "movq $0x599b9401c7d4aa88, %%rax;"
        "mulq 168(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A7*B21  n
        
        "movq $0xa177d0ac7b179728, %%rax;"
        "mulq 160(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A8*B20  n
        
        "movq $0xd528570ce8895ea5, %%rax;"
        "mulq 152(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A9*B19  n
        
        "movq $0x5ef4b87df2399756, %%rax;"
        "mulq 144(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A10*B18  n
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 136(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A11*B17  n
        
        "movq $0xa, %%rax;"
        "mulq 128(%12);"            
        "addq %%rax, %11;"
        "movq %11, %0;"        
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A12*B16  n
        
                                //28
        "xor %11, %11;"
        
        "movq $0x3c8f093f6421c1ec, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A5*B24  n
        
        "movq $0x3c3b41cfa08954ce, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A6*B23  n
        
        "movq $0x599b9401c7d4aa88, %%rax;"
        "mulq 176(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A7*B22  n
        
        "movq $0xa177d0ac7b179728, %%rax;"
        "mulq 168(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A8*B21  n  
        
        "movq $0xd528570ce8895ea5, %%rax;"
        "mulq 160(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A9*B20  n
        
        "movq $0x5ef4b87df2399756, %%rax;"
        "mulq 152(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A10*B19  n
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 144(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A11*B18  n
        
        "movq $0xa, %%rax;"
        "mulq 136(%12);"            
        "addq %%rax, %10;"
        "movq %10, %1;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"         //A12*B17  n
                                //29
        "xor %10, %10;"
        
        
        "movq $0x3c3b41cfa08954ce, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"          //A6*B24  n
                           
        
        "movq $0x599b9401c7d4aa88, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"         //A7*B23  n
        
        
        "movq $0xa177d0ac7b179728, %%rax;"
        "mulq 176(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"         //A8*B22  n
        
        
        "movq $0xd528570ce8895ea5, %%rax;"
        "mulq 168(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"         //A9*B21  n
        
        
        "movq $0x5ef4b87df2399756, %%rax;"
        "mulq 160(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"         //A10*B20  n
        
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 152(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"         //A11*B19  n
        
        
        "movq $0xa, %%rax;"
        "mulq 144(%12);"            
        "addq %%rax, %%r9;"
        "movq %%r9, %2;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"         //A12*B18  n
                                    //30
        "xor %%r9, %%r9;"
        
        "movq $0x599b9401c7d4aa88, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"        //A7*B24  n
        
        "movq $0xa177d0ac7b179728, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"          //A8*B23  n
        
        "movq $0xd528570ce8895ea5, %%rax;"
        "mulq 176(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"          //A9*B22  n
        
        "movq $0x5ef4b87df2399756, %%rax;"
        "mulq 168(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"          //A10*B21  n
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 160(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"          //A11*B20  n
        
        "movq $0xa, %%rax;"
        "mulq 152(%12);"            
        "addq %%rax, %11;"
        "movq %11, %3;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"         //A12*B19  n
                                //31
        "xor %11, %11;"
        
        "movq $0xa177d0ac7b179728, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"         //A8*B24  n
        
        "movq $0xd528570ce8895ea5, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"         //A9*B23  n
        
        "movq $0x5ef4b87df2399756, %%rax;"
        "mulq 176(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"         //A10*B22  n
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 168(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"          //A11*B21  n
        
        "movq $0xa, %%rax;"
        "mulq 160(%12);"            
        "addq %%rax, %10;"
        "movq %10, %4;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"         //A12*B20 n
                                //32
        "xor %10, %10;"
            
        "movq $0xd528570ce8895ea5, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"        //A9*B24  n
        
        
        "movq $0x5ef4b87df2399756, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"        //A10*B23  n
        
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 176(%12);"            
        "addq %%rax, %%r9;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"        //A11*B22  n
        
        
        "movq $0xa, %%rax;"
        "mulq 168(%12);"            
        "addq %%rax, %%r9;"
        "movq %%r9, %5;"
        "adcq %%rdx, %11;"     
        "adcq $0, %10;"        //A12*B21  n
                                //33
        "xor %%r9, %%r9;"
        
        "movq $0x5ef4b87df2399756, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"          //A10*B24  n
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %11;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"        //A11*B23  n
        
        "movq $0xa, %%rax;"
        "mulq 176(%12);"            
        "addq %%rax, %11;"
        "movq %11, %6;"
        "adcq %%rdx, %10;"     
        "adcq $0, %%r9;"        //A12*B22  n
                                //34
        "xor %11, %11;"
        
        "movq $0xd88b672f67f5f8ff, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %10;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"         //A11*B24  n
        
        
        "movq $0xa, %%rax;"
        "mulq 184(%12);"            
        "addq %%rax, %10;"
        "movq %10, %7;"
        "adcq %%rdx, %%r9;"     
        "adcq $0, %11;"        //A12*B23  n
                            //35
        "xor %10, %10;"
        
        "movq $0xa, %%rax;"
        "mulq 192(%12);"            
        "addq %%rax, %%r9;"
        "movq %%r9, %8;"
      //  "adcq %%rdx, %11;"     
       // "adcq $0, %11;"        //A12*B24  n
                                //36 
        
    //    "movq %11, %9;"        //37

        
        :  "=&g" ((res->_mp_d)[28]), "=&g" ((res->_mp_d)[29]), "=&g" ((res->_mp_d)[30]), "=&g" ((res->_mp_d)[31]), "=&g" ((res->_mp_d)[32]), "=&g" ((res->_mp_d)[33]), "=&g" ((res->_mp_d)[34]), "=&g" ((res->_mp_d)[35]), "=&g" ((res->_mp_d)[36]), "=&g" ((res->_mp_d)[37]), "+r" (ra), "+r" (rb)//,  "=&X" ((res->_mp_d)[2]), "=&X" ((res->_mp_d)[3]), "=&X" ((res->_mp_d)[4]),  "=&X" ((res->_mp_d)[5]), "=&X" ((res->_mp_d)[6]), "=&X" ((res->_mp_d)[7]),  "=&X" ((res->_mp_d)[8]), "=&X" ((res->_mp_d)[9]), "=&X" ((res->_mp_d)[10]), "=&X" ((res->_mp_d)[11]), "=&X" ((res->_mp_d)[12]), "=&X" ((res->_mp_d)[13]),  "=&X" ((res->_mp_d)[14]), "=&X" ((res->_mp_d)[15]), "=&X" ((res->_mp_d)[16]),  "=&X" ((res->_mp_d)[17]), "=&X" ((res->_mp_d)[18]), "=&X" ((res->_mp_d)[19]),  "=&X" ((res->_mp_d)[20]), "=&X" ((res->_mp_d)[21]), "=&X" ((res->_mp_d)[22]), "=&X" ((res->_mp_d)[23]), "=&X" ((res->_mp_d)[24]), "=&X" ((res->_mp_d)[25])
        : "r" ((bin->_mp_d))//, "X" ((ain->_mp_d)[1]), "X" ((bin->_mp_d)[1]), "X" ((ain->_mp_d)[2]), "X" ((bin->_mp_d)[2]), "X" ((ain->_mp_d)[3]), "X" ((bin->_mp_d)[3]), "X" ((ain->_mp_d)[4]), "X" ((bin->_mp_d)[4]), "X" ((ain->_mp_d)[5]), "X" ((bin->_mp_d)[5]), "X" ((ain->_mp_d)[6]), "X" ((bin->_mp_d)[6]), "X" ((ain->_mp_d)[7]), "X" ((bin->_mp_d)[7]), "X" ((ain->_mp_d)[8]), "X" ((bin->_mp_d)[8]), "X" ((ain->_mp_d)[9]), "X" ((bin->_mp_d)[9]), "X" ((ain->_mp_d)[10]), "X" ((bin->_mp_d)[10]), "X" ((ain->_mp_d)[11]), "X" ((bin->_mp_d)[11]), "X" ((ain->_mp_d)[12]), "X" ((bin->_mp_d)[12])
        : "%rax", "%rdx", "%r9"
        );
    
    // printf("\nIN FUN res: %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[37], (res->_mp_d)[36], (res->_mp_d)[35], (res->_mp_d)[24], (res->_mp_d)[34], (res->_mp_d)[33], (res->_mp_d)[32], (res->_mp_d)[31], (res->_mp_d)[30], (res->_mp_d)[29], (res->_mp_d)[28], (res->_mp_d)[27], (res->_mp_d)[26], (res->_mp_d)[25], (res->_mp_d)[24], (res->_mp_d)[23], (res->_mp_d)[22], (res->_mp_d)[21], (res->_mp_d)[20], (res->_mp_d)[19], (res->_mp_d)[18], (res->_mp_d)[17], (res->_mp_d)[16], (res->_mp_d)[15], (res->_mp_d)[14], (res->_mp_d)[13], (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
    
     switch (bin->_mp_size) {
             
       case 24:
             if( (res->_mp_d)[36]==0 )
                  res->_mp_size = 36;
             else
                  res->_mp_size = 37;
         
         break;
             
        case 25:
         //    if( (res->_mp_d)[36]==0 )
               //  ;//res->_mp_size = 36;
           //  else
                 res->_mp_size = 37;
             
             break;
             
    /*     case 13:
             if( (res->_mp_d)[25]==0 )
                 res->_mp_size = 25;
             else 
                 res->_mp_size = 26;
     */        
    //     break;
             
     /*    case 12:
             if( (res->_mp_d)[24]==0 )
                 res->_mp_size = 24;
             else 
                 res->_mp_size = 25;
       */      
      //   break;
     }
    
}



//this assembly implementation of addittion assumes the size of your inputs (in limbs) will be 12, 13 or 0
void mpz_addm_x86(mpz_t res, mpz_t ain, mpz_t bin){
    
    
    if( ain->_mp_size==0 && bin->_mp_size==0){
        mpz_set_ui(res,0);
        return;
    }    
    
    mpz_add(res,ain,bin);
    
    //  if( res->_mp_alloc < 13 )
    //      _mpz_realloc(res,13);
    
    // printf("\nIN FUNC: res size: %d\n", res->_mp_size);
    //  gmp_printf("\nIN FUNC: res : %Zd\n", res);
    
    if( mpz_cmp(res,prime) >= 0){
        
        asm   (  
               "movq	$0xffffffffffffffff, %%r8;"
               "subq	%%r8, %0;"
               
               "movq	$0xffffffffffffffff, %%r8;"     //problem line
               "sbbq	%%r8, %1;"   
               
               "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %2;"
               
               "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %3;"
               
               "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %4;"
               
               "movq	$0xffffffffffffffff, %%r8;"
               "sbbq	%%r8, %5;"
               
               "movq   $0xf007669a5ce89647, %%r8;"
               "sbbq	%%r8, %6;"
               
               "movq	$0xade00d91484504f9, %%r8;"
               "sbbq	%%r8, %7;"
               
               "movq	$0x0979d570c24486e3, %%r8;"
               "sbbq	%%r8, %8;"
               
               "movq	$0x8bbae3679a4c7025, %%r8;"
               "sbbq	%%r8, %9;"
               
               "movq	$0xa06a805a9f6808b4, %%r8;"
               "sbbq	%%r8, %10;"
               
               "movq	$0xe69ebefa87fabdfa, %%r8;"
               "sbbq	%%r8, %11;"
               
               "movq	$0x5, %%r8;"
               "sbbq   %%r8, %12;"
               
               
               : "+g" ((res->_mp_d)[0]), "+g" ((res->_mp_d)[1]), "+g" ((res->_mp_d)[2]), "+g" ((res->_mp_d)[3]), "+g" ((res->_mp_d)[4]), "+g" ((res->_mp_d)[5]), "+g" ((res->_mp_d)[6]), "+g" ((res->_mp_d)[7]), "+g" ((res->_mp_d)[8]), "+g" ((res->_mp_d)[9]), "+g" ((res->_mp_d)[10]), "+g" ((res->_mp_d)[11]), "+g" ((res->_mp_d)[12])
               :
               :"%r8"
               );

        
           printf("IN FUNC RESULT:\n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
    }
    
    if( (res->_mp_d)[24]==0 )
        res->_mp_size = 24;
    else
        res->_mp_size = 13;
    
    
    //res->_mp_size = 13;
    
    
    
}

/*
//this assembly implementation of addittion assumes the size of your inputs (in limbs) will be 12, 13 or 0. Experimental one but faster than mpz_addm_x86.
void mpz_addm_x86_t(mpz_t res, mpz_t ain, mpz_t bin, mpz_t tmp1, mpz_t tmp2 ){
    
    mpz_set(tmp1, ain);
    mpz_set(tmp2, bin);
    
    if( tmp2->_mp_size==0 || tmp1->_mp_size==0){
        
        if( tmp1->_mp_size==0 &&  tmp2->_mp_size!=0){
           mpz_set(res,tmp2);
           return;
        }
        
        if( tmp2->_mp_size==0 &&  tmp1->_mp_size!=0){
            mpz_set(res,tmp1);
            return;
        }
        
        //if we get here both tmp1 and tmp2 are exactly 0
        mpz_set_ui(res,0);
        return;
    }
    
    res->_mp_size = 13;
    
//    if( res->_mp_alloc < 13 )
//        _mpz_realloc(res,13);
    
 //   if( tmp1->_mp_alloc < 13 )
 //       _mpz_realloc(tmp1,13);
    
//    if( tmp2->_mp_alloc < 13 )
//        _mpz_realloc(tmp2,13);
    
    if( tmp1->_mp_size==12 ){
     //   _mpz_realloc(tmp1,13);
        (tmp1->_mp_d)[12]=0;
    }
    
    if( tmp2->_mp_size==12 )
        (tmp2->_mp_d)[12]=0;
    
/*    mpz_t tmp;
    mpz_init(&tmp);
    
      mpz_add(tmp,tmp1,tmp2);
     if( mpz_cmp(tmp,prime) >= 0)
         printf("\nMOD\n");
  */  
/*   // printf("\nIN FUNC: res size: %d\n", res->_mp_size);
  //  gmp_printf("\nIN FUNC: res : %Zd\n", res);
    
  //  if( mpz_cmp(res,prime) >= 0){
        
       // printf("\nMOD\n");

       asm   (  
     //    "subq    $208,%%rsp;"
              
         "movq %13, %%r8;"
         "addq %14, %%r8;"
         "push %%r8;" 
             
         "movq %15, %%r8;"
         "adcq %16, %%r8;" 
         "push %%r8;"    
                     
         "movq %17, %%r8;"
         "adcq %18, %%r8;"   
         "push %%r8;"  
                     
          "movq %19, %%r8;"
          "adcq %20, %%r8;"   
          "push %%r8;" 
                     
          "movq %21, %%r8;"
          "adcq %22, %%r8;"   
          "push %%r8;"  
                     
          "movq %23, %%r8;"
          "adcq %24, %%r8;"  
          "push %%r8;"  
                     
          "movq %25, %%r8;"
          "adcq %26, %%r8;"
          "push %%r8;"  
                     
          "movq %27, %%r8;"
          "adcq %28, %%r8;"
          "push %%r8;"  
             
          "movq %29, %%r8;"
          "adcq %30, %%r8;"   
          "push %%r8;"  
                     
          "movq %31, %%r8;"
          "adcq %32, %%r8;"   
          "push %%r8;"  
                     
          "movq %33, %%r8;"
          "adcq %34, %%r8;"
          "push %%r8;"  
                     
          "movq %35, %%r8;"
          "adcq %36, %%r8;" 
          "push %%r8;"   
                     
          "movq %37, %%r8;"
          "adcq %38, %%r8;" 
          "push %%r8;"                     
         
         "push 96(%%rsp);"    
       //  "movq	$0xffffffffffffffff, %%r8;"
         "subq	$0xffffffffffffffff, 104(%%rsp);"
    
         "push 96(%%rsp);" 
         //"movq	$0xffffffffffffffff, %%r8;" 
         "sbbq	$0xffffffffffffffff, 104(%%rsp);"   
   
         "push 96(%%rsp);" 
        // "movq	$0xffffffffffffffff, %%r8;"
         "sbbq	$0xffffffffffffffff, 104(%%rsp);"
    
         "push 96(%%rsp);"      
       //  "movq	$0xffffffffffffffff, %%r8;"
         "sbbq	$0xffffffffffffffff, 104(%%rsp);"
      
         "push 96(%%rsp);" 
        // "movq	$0xffffffffffffffff, %%r8;"
         "sbbq	$0xffffffffffffffff, 104(%%rsp);"
         
         "push 96(%%rsp);"     
       //  "movq	$0xffffffffffffffff, %%r8;"
         "sbbq	$0xffffffffffffffff, 104(%%rsp);"
 
        "push 96(%%rsp);"     
        "movq   $0xf007669a5ce89647, %%r8;"
        "sbbq	%%r8, 104(%%rsp);"
        
        "push 96(%%rsp);"      
        "movq	$0xade00d91484504f9, %%r8;"
        "sbbq	%%r8, 104(%%rsp);"
        
        "push 96(%%rsp);"       
        "movq	$0x0979d570c24486e3, %%r8;"
        "sbbq	%%r8, 104(%%rsp);"

        "push 96(%%rsp);" 
        "movq	$0x8bbae3679a4c7025, %%r8;"
        "sbbq	%%r8, 104(%%rsp);"

        "push 96(%%rsp);" 
        "movq	$0xa06a805a9f6808b4, %%r8;"
        "sbbq	%%r8, 104(%%rsp);"
        
        "push 96(%%rsp);"      
        "movq	$0xe69ebefa87fabdfa, %%r8;"
        "sbbq	%%r8, 104(%%rsp);"
   
        "push 96(%%rsp);" 
     //   "movq	$0x5, %%r8;"
        "sbbq   $0x5, 104(%%rsp);"
           
        "jc end;"
        "addq $0x68, %%rsp;"   //gets to the part of the stack we want
        "end:"
              
        "pop %%r8;"
        "movq %%r8, %12;"
              
        "pop %%r8;"
        "movq %%r8, %11;"
              
        "pop %%r8;"
        "movq %%r8, %10;"
              
        "pop %%r8;"
        "movq %%r8, %9;"
              
        "pop %%r8;"
        "movq %%r8, %8;"
              
        "pop %%r8;"
        "movq %%r8, %7;"
              
        "pop %%r8;"
        "movq %%r8, %6;"
              
        "pop %%r8;"
        "movq %%r8, %5;"
        
        "pop %%r8;"
        "movq %%r8, %4;"
        
        "pop %%r8;"
        "movq %%r8, %3;"
              
        "pop %%r8;"
        "movq %%r8, %2;"
              
        "pop %%r8;"
        "movq %%r8, %1;"
              
        "pop %%r8;"
        "movq %%r8, %0;"
              
        " jnc end2;"  
        "addq $0x68, %%rsp;"  //resets stack
        "end2:"
        
     //   "addq    $208,%%rsp;"
        
       : "=g" ((res->_mp_d)[0]), "=g" ((res->_mp_d)[1]), "=g" ((res->_mp_d)[2]), "=g" ((res->_mp_d)[3]), "=g" ((res->_mp_d)[4]), "=g" ((res->_mp_d)[5]), "=g" ((res->_mp_d)[6]), "=g" ((res->_mp_d)[7]), "=g" ((res->_mp_d)[8]), "=g" ((res->_mp_d)[9]), "=g" ((res->_mp_d)[10]), "=g" ((res->_mp_d)[11]), "=g" ((res->_mp_d)[12])
       : "g" ((tmp1->_mp_d)[0]), "g" ((tmp2->_mp_d)[0]), "g" ((tmp1->_mp_d)[1]), "g" ((tmp2->_mp_d)[1]), "g" ((tmp1->_mp_d)[2]), "g" ((tmp2->_mp_d)[2]), "g" ((tmp1->_mp_d)[3]), "g" ((tmp2->_mp_d)[3]), "g" ((tmp1->_mp_d)[4]), "g" ((tmp2->_mp_d)[4]), "g" ((tmp1->_mp_d)[5]), "g" ((tmp2->_mp_d)[5]), "g" ((tmp1->_mp_d)[6]), "g" ((tmp2->_mp_d)[6]), "g" ((tmp1->_mp_d)[7]), "g" ((tmp2->_mp_d)[7]), "g" ((tmp1->_mp_d)[8]), "g" ((tmp2->_mp_d)[8]), "g" ((tmp1->_mp_d)[9]), "g" ((tmp2->_mp_d)[9]), "g" ((tmp1->_mp_d)[10]), "g" ((tmp2->_mp_d)[10]), "g" ((tmp1->_mp_d)[11]), "g" ((tmp2->_mp_d)[11]), "g" ((tmp1->_mp_d)[12]), "g" ((tmp2->_mp_d)[12])
       :"%r8"
      );
        
    //    printf("IN FUNC RESULT:\n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
 // }
    
    if( (res->_mp_d)[12]==0 )
        res->_mp_size = 12;
    else
        res->_mp_size = 13;
    
    
     //res->_mp_size = 13;
    
   
}
*/
void mpz_mul_x86_4(mpz_t res, mpz_t ain){
    
    unsigned long int ra=0;
    unsigned long int rb=0;
    
 //   mpz_t tmp;
 //   mpz_init(tmp);
  //  mpz_set(tmp,ain);

    
    switch (ain->_mp_size) {
            
        case 0:
            res->_mp_size=0;
            return;
                 
        case 12:
            (ain->_mp_d)[12]=0;
            break;
    }
    
     
    asm (            
         "movq $0xffffffffffffffff, %%rax;" 
         "mulq 0(%9);"
         "movq %%rax, %0;"    
         "movq %%rdx, %7;"           //A0*B0
         //0
         
         "xorq %%r10, %%r10;" 
         
         "movq $0xffffffffffffffff, %%rax;"      //problem line
         "mulq 0(%9);"              
         "addq %%rax, %7;"     
         "movq %%rdx, %8;"  
         "adcq $0, %8;"              //A1*B0
         
         "movq $0xffffffffffffffff, %%rax;"  
         "mulq 8(%9);"         
         "addq %%rax, %7;" 
         "movq %7, %1;"  
         "adcq %%rdx,%8;"    
         "adcq $0, %%r10;"                //A0*B1
         //1
         
         "xorq %7, %7;" 
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"           //A0*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A1*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %8;"    
         "movq %8, %2;" 
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A2*B0
         //2
         "xorq %8, %8;"  
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"              //A3*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"            //A0*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A2*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %%r10;"   
         "movq %%r10, %3;" 
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A1*B2
         //3
         "xorq %%r10, %%r10;" 
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"            //A4*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"            //A0*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A3*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A1*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %7;"  
         "movq %7, %4;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A2*B2
         //4
         "xor %7, %7;"
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A5*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A0*B5
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A4*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A1*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A3*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %8;"
         "movq %8, %5;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A2*B3
         //5
         "xor %8, %8;"
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A6*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A0*B6
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A5*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A1*B5
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"            //A4*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"            //A2*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %%r10;"
         "movq %%r10, %6;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"  
         
         //    "movq %%r8, %7;" 
         //   "movq %%r9, %8;" 
         
         //A3*B3
         //6       //break
         :  "=&g" ((res->_mp_d)[0]), "=&g" ((res->_mp_d)[1]), "=&g" ((res->_mp_d)[2]), "=&g" ((res->_mp_d)[3]), "=&g" ((res->_mp_d)[4]), "=&g" ((res->_mp_d)[5]), "=&g" ((res->_mp_d)[6]), "+r" (ra), "+r" (rb)
         : "r" ((ain->_mp_d))
         : "%rax", "%rdx", "%r10"
         );
    
    asm (
         
         //   "movq %11, %%r8;" 
         //   "movq %12, %%r9;" 
         
         
         //"xor %%r8, %%r8;"
         //"xor %%r9, %%r9;"
         "xor %%r10, %%r10;"
         
         
         "movq $0xade00d91484504f9, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A7*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 56(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A0*B7
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A6*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A1*B6
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A5*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"        //A2*B5
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A4*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %7;"  
         "movq %7, %0;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A3*B4
         //7
         "xor %7, %7;"
         
         "movq $0x0979d570c24486e3, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"       //A8*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 64(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"             //A0*B8
         
         "movq $0xade00d91484504f9, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A7*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 56(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A1*B7
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A6*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A2*B6
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A5*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %8;"            
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A3*B5
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %8;" 
         "movq %8, %1;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A4*B4
         //8 
         
         "xor %8, %8;"
         
         "movq $0x8bbae3679a4c7025, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"         //A9*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 72(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A0*B9
         
         "movq $0x0979d570c24486e3, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"          //A8*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 64(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"         //A1*B8
         
         "movq $0xade00d91484504f9, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"  //A7*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 56(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"         //A2*B7
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A6*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A3*B6
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %%r10;"            
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A5*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %%r10;"
         "movq %%r10, %2;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A4*B5
         //9  
         "xor %%r10, %%r10;"
         
         "movq $0xa06a805a9f6808b4, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A10*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 80(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"        //A0*B10
         
         "movq $0x8bbae3679a4c7025, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A9*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 72(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A1*B9
         
         "movq $0x0979d570c24486e3, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A8*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 64(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A2*B8
         
         "movq $0xade00d91484504f9, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"          //A7*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 56(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A3*B7
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A6*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %7;"            
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"       //A4*B6
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %7;"
         "movq %7, %3;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"         //A5*B5
         //10
         "xor %7, %7;"
         
         "movq $0xe69ebefa87fabdfa, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"       //A11*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 88(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A0*B11
         
         "movq $0xa06a805a9f6808b4, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A10*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 80(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A1*B10
         
         "movq $0x8bbae3679a4c7025, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"            //A9*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 72(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A2*B9
         
         "movq $0x0979d570c24486e3, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A8*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 64(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A3*B8
         
         "movq $0xade00d91484504f9, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A7*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 56(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A4*B7
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %8;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A6*B5
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %8;"
         "movq %8, %4;"
         "adcq %%rdx, %%r10;"
         "adcq $0, %7;"        //A5*B6
         //11
         "xor %8, %8;"
         
         "movq $0x5, %%rax;"
         "mulq 0(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"    //A12*B0
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 96(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"    //A0*B12
         
         
         "movq $0xe69ebefa87fabdfa, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"    //A11*B1
         
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 88(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A1*B11
         
         "movq $0xa06a805a9f6808b4, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A10*B2
         
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 80(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A2*B10
         
         
         "movq $0x8bbae3679a4c7025, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A9*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 72(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A3*B9
         
         
         "movq $0x0979d570c24486e3, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A8*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 64(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A4*B8
         
         
         "movq $0xade00d91484504f9, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A7*B5
         
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 56(%9);"            
         "addq %%rax, %%r10;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A5*B7
         
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %%r10;"
         "movq %%r10, %5;"
         "adcq %%rdx, %7;"
         "adcq $0, %8;"        //A6*B6
         //12
         "xor %%r10, %%r10;"
         
         "movq $0x5, %%rax;"
         "mulq 8(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A12*B1
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 96(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A1*B12
         
         "movq $0xe69ebefa87fabdfa, %%rax;"
         "mulq 16(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A11*B2
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 88(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A2*B11
         
         "movq $0xa06a805a9f6808b4, %%rax;"
         "mulq 24(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A10*B3
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 80(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A3*B10
         
         "movq $0x8bbae3679a4c7025, %%rax;"
         "mulq 32(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A9*B4
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 72(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A4*B9
         
         "movq $0x0979d570c24486e3, %%rax;"
         "mulq 40(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A8*B5
         
         "movq $0xffffffffffffffff, %%rax;"
         "mulq 64(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A5*B8
         
         "movq $0xade00d91484504f9, %%rax;"
         "mulq 48(%9);"            
         "addq %%rax, %7;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;"   //A7*B6
         
         "movq $0xf007669a5ce89647, %%rax;"
         "mulq 56(%9);"            
         "addq %%rax, %7;"
         "movq %7, %6;"
         "adcq %%rdx, %8;"
         "adcq $0, %%r10;" //A6*B7
         //13                  break
         
         //   "movq %8, %7;"
         //   "movq %%r10, %8;"
         
         //    "movq %8, %7;"
         "movq %%r10, %7;"
         
         :  "=&g" ((res->_mp_d)[7]), "=&g" ((res->_mp_d)[8]), "=&g" ((res->_mp_d)[9]), "=&g" ((res->_mp_d)[10]), "=&g" ((res->_mp_d)[11]), "=&g" ((res->_mp_d)[12]), "=&g" ((res->_mp_d)[13]), "+r" (ra), "+r" (rb)
         : "r" ((ain->_mp_d))
         : "%rax", "%rdx", "%r8"
         );
    
    asm(
        
        //   "movq %8, %%r9;"
        //   "movq %7, %%r10;"             
        
        "xor %%r8, %%r8;"
        
        "movq $0x5, %%rax;"
        "mulq 16(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B2
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 96(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A2*B12
        
        "movq $0xe69ebefa87fabdfa, %%rax;"
        "mulq 24(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B3
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 88(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A3*B11
        
        "movq $0xa06a805a9f6808b4, %%rax;"
        "mulq 32(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B4
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 80(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A4*B10
        
        "movq $0x8bbae3679a4c7025, %%rax;"
        "mulq 40(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B5
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 72(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A5*B9
        
        "movq $0x0979d570c24486e3, %%rax;"
        "mulq 48(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B6
        
        "movq $0xf007669a5ce89647, %%rax;"
        "mulq 64(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A6*B8
        
        "movq $0xade00d91484504f9, %%rax;"
        "mulq 56(%9);"            
        "addq %%rax, %8;"
        "movq %8, %0;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A7*B7
        //14 
        "xor %8, %8;"
        
        "movq $0x5, %%rax;"
        "mulq 24(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A12*B3
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 96(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A3*B12
        
        "movq $0xe69ebefa87fabdfa, %%rax;"
        "mulq 32(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A11*B4
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 88(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A4*B11
        
        "movq $0xa06a805a9f6808b4, %%rax;"
        "mulq 40(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A10*B5
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 80(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A5*B10
        
        "movq $0x8bbae3679a4c7025, %%rax;"
        "mulq 48(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A9*B6
        
        "movq $0xf007669a5ce89647, %%rax;"
        "mulq 72(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A6*B9
        
        "movq $0x0979d570c24486e3, %%rax;"
        "mulq 56(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A8*B7
        
        "movq $0xade00d91484504f9, %%rax;"
        "mulq 64(%9);"            
        "addq %%rax, %7;"
        "movq %7, %1;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;" //A7*B8
        //15 
        "xor %7, %7;"
        
        "movq $0x5, %%rax;"
        "mulq 32(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A12*B4
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 96(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A4*B12
        
        "movq $0xe69ebefa87fabdfa, %%rax;"
        "mulq 40(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A11*B5
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 88(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A5*B11
        
        "movq $0xa06a805a9f6808b4, %%rax;"
        "mulq 48(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A10*B6
        
        "movq $0xf007669a5ce89647, %%rax;"
        "mulq 80(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A6*B10
        
        "movq $0x8bbae3679a4c7025, %%rax;"
        "mulq 56(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A9*B7
        
        "movq $0xade00d91484504f9, %%rax;"
        "mulq 72(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A7*B9
        
        "movq $0x0979d570c24486e3, %%rax;"
        "mulq 64(%9);"            
        "addq %%rax, %%r8;"
        "movq %%r8, %2;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A8*B8
        //16
        "xor %%r8, %%r8;"
        
        "movq $0x5, %%rax;"
        "mulq 40(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B5
        
        "movq $0xffffffffffffffff, %%rax;"
        "mulq 96(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A5*B12
        
        "movq $0xe69ebefa87fabdfa, %%rax;"
        "mulq 48(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B6
        
        "movq $0xf007669a5ce89647, %%rax;"
        "mulq 88(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A6*B11
        
        "movq $0xa06a805a9f6808b4, %%rax;"
        "mulq 56(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B7
        
        "movq $0xade00d91484504f9, %%rax;"
        "mulq 80(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A7*B10
        
        "movq $0x8bbae3679a4c7025, %%rax;"
        "mulq 64(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B8
        
        "movq $0x0979d570c24486e3, %%rax;"
        "mulq 72(%9);"            
        "addq %%rax, %8;"
        "movq %8, %3;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B9
        //17 
        "xor %8, %8;"
        
        "movq $0x5, %%rax;"
        "mulq 48(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A12*B6
        
        "movq $0xf007669a5ce89647, %%rax;"
        "mulq 96(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A6*B12
        
        "movq $0xe69ebefa87fabdfa, %%rax;"
        "mulq 56(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A11*B7
        
        "movq $0xade00d91484504f9, %%rax;"
        "mulq 88(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A7*B11
        
        "movq $0xa06a805a9f6808b4, %%rax;"
        "mulq 64(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A10*B8
        
        "movq $0x0979d570c24486e3, %%rax;"
        "mulq 80(%9);"            
        "addq %%rax, %7;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A8*B10
        
        "movq $0x8bbae3679a4c7025, %%rax;"
        "mulq 72(%9);"            
        "addq %%rax, %7;"
        "movq %7, %4;"
        "adcq %%rdx, %%r8;"
        "adcq $0, %8;"    //A9*B9
        //18
        "xor %7, %7;"
        
        "movq $0x5, %%rax;"
        "mulq 56(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A12*B7
        
        "movq $0xade00d91484504f9, %%rax;"
        "mulq 96(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A7*B12
        
        "movq $0xe69ebefa87fabdfa, %%rax;"
        "mulq 64(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A11*B8
        
        "movq $0x0979d570c24486e3, %%rax;"
        "mulq 88(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A8*B11
        
        "movq $0xa06a805a9f6808b4, %%rax;"
        "mulq 72(%9);"            
        "addq %%rax, %%r8;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A10*B9
        
        "movq $0x8bbae3679a4c7025, %%rax;"
        "mulq 80(%9);"            
        "addq %%rax, %%r8;"
        "movq %%r8, %5;"
        "adcq %%rdx, %8;"
        "adcq $0, %7;"   //A9*B10
        //19
        "xor %%r8, %%r8;"
        
        
        "movq $0x5, %%rax;"
        "mulq 64(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A12*B8
        
        "movq $0x0979d570c24486e3, %%rax;"
        "mulq 96(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A8*B12
        
        "movq $0xe69ebefa87fabdfa, %%rax;"
        "mulq 72(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A11*B9
        
        "movq $0x8bbae3679a4c7025, %%rax;"
        "mulq 88(%9);"            
        "addq %%rax, %8;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A9*B11
        
        "movq $0xa06a805a9f6808b4, %%rax;"
        "mulq 80(%9);"            
        "addq %%rax, %8;"
        "movq %8, %6;"
        "adcq %%rdx, %7;"
        "adcq $0, %%r8;"    //A10*B10
        //20           break
        
        // "movq %7, %8;"
        "movq %%r8, %8;"
        
        :  "=&g" ((res->_mp_d)[14]), "=&g" ((res->_mp_d)[15]), "=&g" ((res->_mp_d)[16]), "=&g" ((res->_mp_d)[17]), "=&g" ((res->_mp_d)[18]), "=&g" ((res->_mp_d)[19]), "=&g" ((res->_mp_d)[20]), "+r" (ra), "+r" (rb)
        :  "r" ((ain->_mp_d))
        : "%rax", "%rdx", "%r8"
        );
    
    asm (
         
         //   "movq %5, %%r8;"
         //   "movq %6, %%r10;"
         
         
         "xor %%r9, %%r9;"
         
         
         "movq $0x5, %%rax;"
         "mulq 72(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A12*B9
         
         "movq  $0x8bbae3679a4c7025, %%rax;"
         "mulq 96(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A9*B12
         
         "movq $0xe69ebefa87fabdfa, %%rax;"
         "mulq 80(%6);"            
         "addq %%rax, %5;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A11*B10
         
         "movq $0xa06a805a9f6808b4, %%rax;"
         "mulq 88(%6);"            
         "addq %%rax, %5;"
         "movq %5, %0;"
         "adcq %%rdx, %4;"
         "adcq $0, %%r9;"    //A10*B11
         //21
         "xor %5, %5;"
         
         
         "movq $0x5, %%rax;"
         "mulq 80(%6);"            
         "addq %%rax, %4;"
         "adcq %%rdx, %%r9;"
         "adcq $0, %5;"   //A12*B10
         
         "movq $0xa06a805a9f6808b4, %%rax;"
         "mulq 96(%6);"            
         "addq %%rax, %4;"
         "adcq %%rdx, %%r9;"
         "adcq $0, %5;"   //A10*B12
         
         "movq $0xe69ebefa87fabdfa, %%rax;"
         "mulq 88(%6);"            
         "addq %%rax, %4;"
         "movq %4, %1;"
         "adcq %%rdx, %%r9;"
         "adcq $0, %5;"   //A11*B11
         //22
         //  "xor %4, %4;"
         
         "movq $0x5, %%rax;"
         "mulq 88(%6);"            
         "addq %%rax, %%r9;"
         "adcq %%rdx, %5;"
         //   "adcq $0, %4;"    //A12*B11
         
         "movq $0xe69ebefa87fabdfa, %%rax;"
         "mulq 96(%6);"            
         "addq %%rax, %%r9;"
         "movq %%r9, %2;"
         "adcq %%rdx, %5;"
         //   "adcq $0, %4;"    //A11*B12
         //23
         //   "xor %%r9, %%r9;"
         
         "movq $0x5, %%rax;"
         "mulq 96(%6);"            
         "addq %%rax, %5;"
         "movq %5, %3;"          //A12*B12
         //      "adcq %%rdx, %4;"     //24
         //      "adcq $0, %4;"    
         
         
         //     "movq %5, %4;" //25
         
         
         :  "=&g" ((res->_mp_d)[21]), "=&g" ((res->_mp_d)[22]), "=&g" ((res->_mp_d)[23]), "=&g" ((res->_mp_d)[24]), "+r" (rb), "+r" (ra)
         : "r" ((ain->_mp_d))
         : "%rax", "%rdx", "%r9"
         );
    
    
    if( (res->_mp_d)[24]==0)
        res->_mp_size = 24;
    else
        res->_mp_size = 25;
    
}


//this assembly implementation of addittion assumes the size of your inputs (in limbs) will be 12, 13 or 0 and ain can also be exactly 1
void mpz_addm_x86_2(mpz_t res, mpz_t ain, mpz_t bin){
    
    if( ain->_mp_size==0 && bin->_mp_size==0){
        mpz_set_ui(res,0);
        return;
    }else  if( ain->_mp_size==1 && bin->_mp_size==0){
    
        mpz_set_ui(res,1);
        return;
    }
    
       
    mpz_add(res,ain,bin);
    
  //  if( res->_mp_alloc < 13 )
  //      _mpz_realloc(res,13);
    
    // printf("\nIN FUNC: res size: %d\n", res->_mp_size);
    //  gmp_printf("\nIN FUNC: res : %Zd\n", res);
    
    if( mpz_cmp(res,prime) >= 0){
        
        asm   (  
                       "movq	$0xffffffffffffffff, %%r8;"
                       "subq	%%r8, %0;"
                 
                       "movq	$0xffffffffffffffff, %%r8;"     //problem line
                       "sbbq	%%r8, %1;"   
                      
                       "movq	$0xffffffffffffffff, %%r8;"
                       "sbbq	%%r8, %2;"
                       
                        "movq	$0xffffffffffffffff, %%r8;"
                       "sbbq	%%r8, %3;"
                       
                        "movq	$0xffffffffffffffff, %%r8;"
                       "sbbq	%%r8, %4;"
                
                       "movq	$0xffffffffffffffff, %%r8;"
                       "sbbq	%%r8, %5;"
                       
                        "movq   $0xf007669a5ce89647, %%r8;"
                       "sbbq	%%r8, %6;"
                     
                       "movq	$0xade00d91484504f9, %%r8;"
                       "sbbq	%%r8, %7;"
                       
                       "movq	$0x0979d570c24486e3, %%r8;"
                       "sbbq	%%r8, %8;"

                       "movq	$0x8bbae3679a4c7025, %%r8;"
                       "sbbq	%%r8, %9;"

                       "movq	$0xa06a805a9f6808b4, %%r8;"
                       "sbbq	%%r8, %10;"

                       "movq	$0xe69ebefa87fabdfa, %%r8;"
                       "sbbq	%%r8, %11;"

                       "movq	$0x5, %%r8;"
                       "sbbq   %%r8, %12;"

                       
                       : "+g" ((res->_mp_d)[0]), "+g" ((res->_mp_d)[1]), "+g" ((res->_mp_d)[2]), "+g" ((res->_mp_d)[3]), "+g" ((res->_mp_d)[4]), "+g" ((res->_mp_d)[5]), "+g" ((res->_mp_d)[6]), "+g" ((res->_mp_d)[7]), "+g" ((res->_mp_d)[8]), "+g" ((res->_mp_d)[9]), "+g" ((res->_mp_d)[10]), "+g" ((res->_mp_d)[11]), "+g" ((res->_mp_d)[12])
                       :
                       :"%r8"
                       );
        
        //    printf("IN FUNC RESULT:\n %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res->_mp_d)[12], (res->_mp_d)[11], (res->_mp_d)[10], (res->_mp_d)[9], (res->_mp_d)[8], (res->_mp_d)[7], (res->_mp_d)[6], (res->_mp_d)[5], (res->_mp_d)[4], (res->_mp_d)[3], (res->_mp_d)[2], (res->_mp_d)[1], (res->_mp_d)[0]);
    }
    
    if( (res->_mp_d)[12]==0 )
        res->_mp_size = 12;
    else
        res->_mp_size = 13;
    
    
    //res->_mp_size = 13;
    
    
}


/***** ALGORITHMS RELATED TO ECC **********/

// One step of Montgomery ladder
void mont_ladder(GF *res1x, GF *res1z,
		 GF *res2x, GF *res2z,
		 const GF x1, const GF z1,
		 const GF x2, const GF z2,
		 const GF dx, const GF dz,
		 const GF A24) {
  GF* tmp = x1.parent->GFtmp;

  add_GF(&tmp[4], x1, z1);         // a = (self.x + self.z)
  sub_GF(&tmp[5], x1, z1);         // b = (self.x - self.z)
  sub_GF(&tmp[6], x2, z2);
  add_GF(&tmp[7], x2, z2);
  sqr_GF(&tmp[1], tmp[4]);         // aa = a.square()
  sqr_GF(&tmp[3], tmp[5]);         // bb = b.square()
  sub_GF(&tmp[0], tmp[1], tmp[3]); // e = aa - bb
  mul_GF(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
  mul_GF(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b

  add_GF(&tmp[2], tmp[6], tmp[7]);
  sqr_GF(&tmp[2], tmp[2]);
  mul_GF(&tmp[2], tmp[2], dz);     // x2 = diff.z*(da + cb).square()

  sub_GF(&tmp[8], tmp[6], tmp[7]);
  sqr_GF(&tmp[8], tmp[8]);
  mul_GF(&tmp[8], tmp[8], dx);     // z2 = diff.x*(da - cb).square()

  mul_GF(res1z, A24, tmp[0]);
  add_GF(res1z, *res1z, tmp[3]);
  mul_GF(res1z, *res1z, tmp[0]);   // z1 = e*(bb + self.curve.A24*e))
  mul_GF(res1x, tmp[1], tmp[3]);  // x1 = aa*bb

  copy_GF(res2x, tmp[2]);
  copy_GF(res2z, tmp[8]);
}

/* Montgomery point doubling */
void mont_double(GF *resx, GF *resz,
		 const GF x, const GF z,
		 const GF A24) {

  GF* tmp = x.parent->GFtmp;
  add_GF(&tmp[0], x, z);           // a = (x + z)
  sqr_GF(&tmp[1], tmp[0]);         // aa = a^2
  sub_GF(&tmp[2], x, z);           // b = (x - z)
  sqr_GF(&tmp[3], tmp[2]);         // bb = b^2
  sub_GF(&tmp[4], tmp[1], tmp[3]); // c = aa - bb
  mul_GF(resz, A24, tmp[4]);   
  add_GF(resz, *resz, tmp[3]);			
  mul_GF(resz, *resz, tmp[4]);   // Z = c (bb + A24 c))	
  mul_GF(resx, tmp[1], tmp[3]);  // X = aa bb  
}

/* Montgomery point tripling */
void mont_triple(GF *resx, GF *resz,
		 const GF x, const GF z,
		 const GF A24) {
  GF* tmp = x.parent->GFtmp;

  // Very dirty function, assuming that mont_double uses
  // registers 0 and 2 to store resp. x+z and x-z
  mont_double(&tmp[5], &tmp[6], x, z, A24);

  sub_GF(&tmp[7], tmp[5], tmp[6]);
  add_GF(&tmp[8], tmp[5], tmp[6]);
  mul_GF(&tmp[5], tmp[7], tmp[0]); // da = (x2 - z2)*a
  mul_GF(&tmp[6], tmp[8], tmp[2]); // cb = (x2 + z2)*b

  add_GF(&tmp[7], tmp[5], tmp[6]);
  sqr_GF(&tmp[7], tmp[7]);
  mul_GF(&tmp[7], tmp[7], z);     // X = z*(da + cb)^2

  sub_GF(&tmp[8], tmp[5], tmp[6]);
  sqr_GF(&tmp[8], tmp[8]);
  mul_GF(resz, tmp[8], x);     // Z = x*(da - cb)^2

  copy_GF(resx, tmp[7]);
}

/*
  Converts a Montgomery point to an Edwards point.
  Guarantee: avoids the temporary registers 0-4.
*/
void mont_to_ed(GF* Rx, GF* Ry,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  /*
    X = x(x+z) / y(x+z)
    Y = y(x-z) / y(x+z)
  */
  add_GF(&tmp[5], x, z);
  sub_GF(&tmp[6], x, z);
  mul_GF(&tmp[7], y, tmp[5]);
  inv_GF(&tmp[8], tmp[7]);
  mul_GF(&tmp[7], tmp[5], tmp[8]);
  mul_GF(Rx, x, tmp[7]);
  mul_GF(&tmp[7], tmp[6], tmp[8]);
  mul_GF(Ry, y, tmp[7]);
}

typedef struct {
    GF tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, dz, dx, x, z;
} thread_p1;

/*void *p3(void *thr){

    thread_p1 *t = (thread_p1 *) thr;
    
   // print_GF(t->tmp2,"tmp2");
    
    /* P3 */
  //  mul_GF(&(t->tmp2), t->tmp0, t->tmp4); // da = (P.x - P.z)*a
/*    mul_GF(&(t->tmp3), t->tmp1, t->tmp5); // cb = (P.x + P.z)*b
    
    add_GF(&(t->tmp8), t->tmp2, t->tmp3);
    sqr_GF(&(t->tmp0), t->tmp8);
    mul_GF(&(t->x), t->tmp0, t->dz);     // x2 = diff.z*(da + cb).square()
   
  /*  sub_GF(&(t->tmp8), t->tmp2, t->tmp3);
    sqr_GF(&(t->tmp0), t->tmp8);
    mul_GF(&(t->z), t->tmp0, t->dx);     // z2 = diff.x*(da - cb).square()
  */  
    
 //  pthread_exit(0);
//}

// Three-point ladder addition step:
//   P1 = 2 P1
//   P2 = dadd(P1, P2, D2)
//   P3 = dadd(P1, P3, D3)
void mont_tradd(GF *x1, GF *z1,
		GF *x2, GF *z2,
		GF *x3, GF *z3,
		const GF dx2, const GF dz2,
		const GF dx3, const GF dz3,
		const GF A24) {
        dieter++;
    
    /*   GF_params *parent1, *parent2, *parent3;
     char * prime = malloc(sizeof(char)*5000);
     mpz_get_str (prime, 10, (x1->parent)->p);
     parent1 = malloc(sizeof(GF_params));
     parent2 = malloc(sizeof(GF_params));
     parent3 = malloc(sizeof(GF_params));
     setup_GF(parent1, prime);
     setup_GF(parent2, prime);
     setup_GF(parent3, prime);
     
     GF* tmpB = parent1->GFtmp;
     GF* tmpC = parent2->GFtmp;*/
    GF* tmp = x1->parent->GFtmp;
    
   /* if(dieter==version){ 
        print_GF(*x3,"\n\nx3\n\n"); 
        print_GF(*z3,"\n\nz3\n\n"); 
        print_GF(*x1,"\n\nx1\n\n"); 
        print_GF(*z1,"\n\nz1\n\n");
        print_GF(*x2,"\n\nx2\n\n"); 
        print_GF(*z2,"\n\nz2\n\n");
    }
    */
    
    // The use of temporary registers in this function is tailored
    // so that it is compatible with
    //       {x1, x2} == {&tmp[0], &tmp[2]}
    //       {z1, z2} == {&tmp[1], &tmp[3]}
    // so that mont_3ladder may safely call it.
    // Be careful when you change indices!
    
    add_GF(&tmp[4], *x1, *z1);         // a = (self.x + self.z)
    sub_GF(&tmp[5], *x1, *z1);         // b = (self.x - self.z)
    sub_GF(&tmp[6], *x2, *z2);
    add_GF(&tmp[7], *x2, *z2);
    sub_GF(&tmp[0], *x3, *z3);
    add_GF(&tmp[1], *x3, *z3); 
    
   /*if(count1==1600){  
        
        print_GF(tmp[6],"\n\ntmp6");
        print_GF(tmp[4],"\n\ntmp4");
        print_GF(tmp[7],"\n\ntmp7");
        print_GF(tmp[5],"\n\ntmp5");
    }*/
    
    // print_GF(tmp[1], "\n\ntmp1");
    
  //  if(dieter==version)
  //      print_GF((tmp[4]),"\n\ntmp4_bt\n\n"); 
    
    /*  copy_GF(&tmpB[4],tmpA[4]);
     copy_GF(&tmpB[5],tmpA[5]);
     
     copy_GF(&tmpC[4],tmpA[4]);
     copy_GF(&tmpC[5],tmpA[5]);
     
     pthread_t thread1;
     thread_p1 thread_params = {tmpA[0], tmpA[1], tmpA[2], tmpA[3], tmpA[4], tmpA[5], tmpA[6], tmpA[7], tmpA[8], dz3, dx3, *x3, *z3}; 
     int rc = pthread_create(&thread1, NULL, p3, (void *) &thread_params); */ 
    
    /* P3 */
    
   /* if(dieter==version){   
        print_GF((tmp[0]),"\n\ntmp0\n\n");
        print_GF((tmp[4]),"\n\ntmp4\n\n");
        print_GF((tmp[1]),"\n\ntmp1\n\n");
        print_GF((tmp[5]),"\n\ntmp5\n\n");
    }*/
    
    mul_GF(&tmp[2], tmp[0], tmp[4]); // da = (P.x - P.z)*a
    mul_GF(&tmp[3], tmp[1], tmp[5]); // cb = (P.x + P.z)*b
    
/*    if(dieter==version){   
        print_GF(tmp[2],"\n\ntmp2\n\n");
        print_GF(tmp[3],"\n\ntmp3\n\n");
    }  */
    
    add_GF(&tmp[8], tmp[2], tmp[3]);
    sqr_GF(&tmp[0], tmp[8]);
    mul_GF(x3, tmp[0], dz3);     // x2 = diff.z*(da + cb).square()
    
    sub_GF(&tmp[8], tmp[2], tmp[3]);
    sqr_GF(&tmp[0], tmp[8]);
    mul_GF(z3, tmp[0], dx3);     // z2 = diff.x*(da - cb).square()
    
    /* P2 */
    mul_GF(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
    
  //  print_GF(tmp[7],"\nx");
  //  print_GF(tmp[5],"\ny");
    
    mul_GF(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b  //6th execution of mul_GF
    
    add_GF(&tmp[2], tmp[6], tmp[7]);
    sqr_GF(&tmp[3], tmp[2]);
    mul_GF(x2, tmp[3], dz2);     // x2 = diff.z*(da + cb).square()
    
    sub_GF(&tmp[3], tmp[6], tmp[7]);
    sqr_GF(&tmp[8], tmp[3]);
    mul_GF(z2, tmp[8], dx2);     // z2 = diff.x*(da - cb).square()
    
    
    /* P1 */
    sqr_GF(&tmp[6], tmp[4]);         // aa = a.square()
    sqr_GF(&tmp[7], tmp[5]);         // bb = b.square()
    sub_GF(&tmp[8], tmp[6], tmp[7]); // e = aa - bb
    mul_GF(&tmp[4], A24, tmp[8]);
    add_GF(&tmp[5], tmp[4], tmp[7]);
    mul_GF(z1, tmp[5], tmp[8]);      // z1 = e*(bb + self.curve.A24*e))
    mul_GF(x1, tmp[6], tmp[7]);      // x1 = aa*bb
    
     count1++;
    
 /*   if(count1==1000){ 
         printf("\n********VALUES TO CHECK***********%d\n",count1);
        print_GF(*(x1), "\n\nt.x1");
        print_GF(*(z1), "\n\nt.z1");
        print_GF(*(x2), "\n\nt.x2");
        print_GF(*(z2), "\n\nt.z2");
        print_GF(*x3,"\n\nt.x3\n\n");
        print_GF(*z3,"\n\nt.z3\n\n"); 
        printf("\n*****************************\n");

        
    }*/
   
   
}


// 3-point ladder to compute P + [t]Q
// Inputs: t, P, Q, Q - P
void mont_3ladder(GF* Rx, GF* Rz,
                  const mpz_t t,
                  const GF Px, const GF Pz,
                  const GF Qx, const GF Qz,
                  const GF QPx, const GF QPz,
                  const GF A24) {
    GF* tmp = Px.parent->GFtmp;
    set_GF(&tmp[0], "0", "1"); 
    set_GF(&tmp[1], "0", "0");
    copy_GF(&tmp[2], Qx); copy_GF(&tmp[3], Qz);
    copy_GF(Rx, Px); copy_GF(Rz, Pz); 
    int bit = mpz_sizeinbase(t, 2) - 1;
    for ( ; bit >=0 ; bit--) {
        if (mpz_tstbit(t, bit) == 0) {
            mont_tradd(&tmp[0], &tmp[1], &tmp[2], &tmp[3], Rx, Rz,
                       Qx, Qz, Px, Pz, A24);
        } else {
            mont_tradd(&tmp[2], &tmp[3], &tmp[0], &tmp[1], Rx, Rz,
                       Qx, Qz, QPx, QPz, A24);
        }
    }
}




/*void *EdwardsCompute(void *thread){
    
    thread_params *thread1;
    thread1 = (thread_params *) thread; 
    
    init_GF(thread1->aPx, (thread1->field)); init_GF(thread1->aPy, (thread1->field));
    mont_to_ed(thread1->aPx, thread1->aPy, *(thread1->Px), *(thread1->Py), *(thread1->Pz));
    init_GF(thread1->aQx, (thread1->field)); init_GF(thread1->aQy, (thread1->field));
    mont_to_ed(thread1->aQx, thread1->aQy, *(thread1->Qx), *(thread1->Qy), *(thread1->Qz));
    
    mul_GF((thread1->tmp4), *(thread1->aPx), *(thread1->aQx)); // tmp4 = C = aPx * aQx
    mul_GF((thread1->tmp5), *(thread1->aPy), *(thread1->aQy)); // tmp5 = D = aPy * aQy
    add_GF((thread1->tmp0), *(thread1->aPx), *(thread1->aPy)); // tmp0 = A = aPx + aPy
    add_GF((thread1->tmp2), *(thread1->aQx), *(thread1->aQy)); // tmp2 = B = aQx + aQy
    mul_GF((thread1->tmp7), *(thread1->tmp4), *(thread1->tmp5));
    mul_GF((thread1->tmp6), *(thread1->d), *(thread1->tmp7)); // tmp6 = E = d * aPx * aQx * aPy * aQy
    sqr_GF((thread1->tmp8), *(thread1->tmp6));
    neg_GF((thread1->tmp7), *(thread1->tmp8));
  /*  add_GF_ui(tmp[7], tmp[7], 1);
    inv_GF(tmp[8], tmp[7]); // tmp8 = 1 / (1-E^2)
    add_GF_ui(tmp[6], tmp[6], 1);
    mul_GF(tmp[7], a, tmp[4]);
    sub_GF(tmp[1], tmp[5], tmp[7]);
    mul_GF(tmp[7], tmp[6], tmp[1]);
    init_GF(PQy, field);
    mul_GF(PQy, tmp[7], tmp[8]); // PQy = (1+E)(D - a C) / (1-E^2)
    neg_GF(tmp[6], tmp[6]);
    add_GF_ui(tmp[6], tmp[6], 2);
    mul_GF(tmp[1], tmp[0], tmp[2]);
    sub_GF(tmp[3], tmp[1], tmp[4]);
    sub_GF(tmp[1], tmp[3], tmp[5]);
    mul_GF(tmp[7], tmp[6], tmp[1]);
    init_GF(PQx, field);
    mul_GF(PQx, tmp[7], tmp[8]);*/ // PQx = (1-E)(A B - C - D) / (1-E^2)

//}

/*
  Computes [m]P + [n]Q, with P and Q points on the Montgomery curve
  with parameters A,B.  Uses Edwards' coordinates for
  calculations.  */
void shamir(GF* Rx, GF* Ry, GF* Rz,
	    const GF A, const GF B,
	    const GF Px, const GF Py, const GF Pz,
	    const GF Qx, const GF Qy, const GF Qz,
	    const mpz_t m, const mpz_t n) {
  
	
   
  // some temporary registers
  GF_params* field = A.parent;
  GF* tmp = field->GFtmp;
  // some other dynamically allocated registers 
  GF a, d, aPx, aPy, aQx, aQy, PQx, PQy;
  
 //  thread_params thread1 = {&aPx, &aPy, &aQx, &aQy, &Px, &Py, &Pz, &Qx, &Qy, &Qz, &field, &d, &a, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6], &tmp[7], &tmp[8]};
    
 // pthread_t threads;  
 // int rc = pthread_create(&threads, NULL, EdwardsCompute, (void *) &thread1);   
	
  /* 
     Parameters of the Edwards curve equivalent to this one:
       a = (A+2)/B
       d = (A-2)/B
  */
  inv_GF(&tmp[0], B);
  copy_GF(&tmp[1], A);
  add_GF_ui(&tmp[1], tmp[1], 2);
  init_GF(&a, field);
  mul_GF(&a, tmp[1], tmp[0]);
  sub_GF_ui(&tmp[1], tmp[1], 4);
  init_GF(&d, field);
  mul_GF(&d, tmp[1], tmp[0]);

  /*
    Computing the Ewdards coordinates of P and Q:
      aPx, aPy = Edwards(P)
      aQx, aQy = Edwards(Q)
  */
  init_GF(&aPx, field); init_GF(&aPy, field);
  mont_to_ed(&aPx, &aPy, Px, Py, Pz);
  init_GF(&aQx, field); init_GF(&aQy, field);
  mont_to_ed(&aQx, &aQy, Qx, Qy, Qz);

  /*
    Computing P+Q using affine Edwards.
  */
  mul_GF(&tmp[4], aPx, aQx); // tmp4 = C = aPx * aQx
  mul_GF(&tmp[5], aPy, aQy); // tmp5 = D = aPy * aQy
  add_GF(&tmp[0], aPx, aPy); // tmp0 = A = aPx + aPy
  add_GF(&tmp[2], aQx, aQy); // tmp2 = B = aQx + aQy
  mul_GF(&tmp[7], tmp[4], tmp[5]);
  mul_GF(&tmp[6], d, tmp[7]); // tmp6 = E = d * aPx * aQx * aPy * aQy
  sqr_GF(&tmp[8], tmp[6]);
  neg_GF(&tmp[7], tmp[8]);
  add_GF_ui(&tmp[7], tmp[7], 1);
  inv_GF(&tmp[8], tmp[7]); // tmp8 = 1 / (1-E^2)
  add_GF_ui(&tmp[6], tmp[6], 1);
  mul_GF(&tmp[7], a, tmp[4]);
  sub_GF(&tmp[1], tmp[5], tmp[7]);
  mul_GF(&tmp[7], tmp[6], tmp[1]);
  init_GF(&PQy, field);
  mul_GF(&PQy, tmp[7], tmp[8]); // PQy = (1+E)(D - a C) / (1-E^2)
  neg_GF(&tmp[6], tmp[6]);
  add_GF_ui(&tmp[6], tmp[6], 2);
  mul_GF(&tmp[1], tmp[0], tmp[2]);
  sub_GF(&tmp[3], tmp[1], tmp[4]);
  sub_GF(&tmp[1], tmp[3], tmp[5]);
  mul_GF(&tmp[7], tmp[6], tmp[1]);
  init_GF(&PQx, field);
  mul_GF(&PQx, tmp[7], tmp[8]); // PQx = (1-E)(A B - C - D) / (1-E^2)*/
  
  int bit = MAX(mpz_sizeinbase(m, 2), mpz_sizeinbase(n, 2)) - 1;
  mpz_set_ui(Rx->a, 0); mpz_set_ui(Ry->a, 0); mpz_set_ui(Rz->a, 0);   
  mpz_set_ui(Rx->b, 0); mpz_set_ui(Ry->b, 1); mpz_set_ui(Rz->b, 1);
  Rx->parent = Ry->parent = Rz->parent = Px.parent;

  for ( ; bit >=0 ; bit--){
    /* Double, using projective Edwards */
    add_GF(&tmp[1], *Rx, *Ry);
    sqr_GF(&tmp[0], tmp[1]); // tmp0 = B = (Rx + Ry)^2
    sqr_GF(&tmp[1], *Rx); // tmp1 = C = Rx^2
    sqr_GF(&tmp[2], *Ry); // tmp2 = D = Ry^2
    mul_GF(&tmp[3], a, tmp[1]); // tmp3 = E = a C
    add_GF(&tmp[4], tmp[3], tmp[2]); // tmp4 = F = E + D
	sqr_GF(&tmp[5], *Rz); // tmp5 = H = Rz^2
    scalar_GF_si(&tmp[7], tmp[5], 2);
    sub_GF(&tmp[6], tmp[4], tmp[7]); // tmp6 = J = F - 2H
    sub_GF(&tmp[7], tmp[0], tmp[1]);
    sub_GF(&tmp[8], tmp[7], tmp[2]);
    mul_GF(Rx, tmp[8], tmp[6]); // Rx = (B-C-D) J
    sub_GF(&tmp[7], tmp[3], tmp[2]);
    mul_GF(Ry, tmp[7], tmp[4]); // Ry = (E-D) F
    mul_GF(Rz, tmp[4], tmp[6]); // Rz = F J

    /* Double and Add, using projective Edwards */
    int r = mpz_tstbit(m, bit) | (mpz_tstbit(n, bit) << 1);
    if (r) {
      if (r == 1) {
	mul_GF(&tmp[0], *Rx, aPx); // tmp0 = C = Rx aPx
	mul_GF(&tmp[1], *Ry, aPy); // tmp1 = D = Ry aPy
	add_GF(&tmp[2], aPx, aPy);  // tmp2 = H = aPx + aPy
      } else if (r == 2) {
	mul_GF(&tmp[0], *Rx, aQx); // tmp0 = C = Rx aQx
	mul_GF(&tmp[1], *Ry, aQy); // tmp1 = D = Ry aQy
	add_GF(&tmp[2], aQx, aQy);  // tmp2 = H = aQx + aQy
      } else {
	mul_GF(&tmp[0], *Rx, PQx); // tmp0 = C = Rx PQx
	mul_GF(&tmp[1], *Ry, PQy); // tmp1 = D = Ry PQy
	add_GF(&tmp[2], PQx, PQy);  // tmp2 = H = PQx + PQy
      }
      sqr_GF(&tmp[3], *Rz); // tmp3 = B = Rz^2
      mul_GF(&tmp[5], tmp[0], tmp[1]);
      mul_GF(&tmp[4], d, tmp[5]); // tmp4 = E = d C D
      sub_GF(&tmp[5], tmp[3], tmp[4]); // tmp5 = F = B - E
      add_GF(&tmp[6], tmp[3], tmp[4]); // tmp6 = G = B + E
      add_GF(&tmp[7], *Rx, *Ry);
      mul_GF(&tmp[8], tmp[7], tmp[2]);
      sub_GF(&tmp[7], tmp[8], tmp[0]);
      sub_GF(&tmp[8], tmp[7], tmp[1]);
      mul_GF(&tmp[7], tmp[8], tmp[5]);
      mul_GF(Rx, *Rz, tmp[7]); // Rx = Rz F ((Rx+Ry)H - C - D)
      mul_GF(&tmp[7], a, tmp[0]);
      sub_GF(&tmp[8], tmp[1], tmp[7]);
      mul_GF(&tmp[7], tmp[6], tmp[8]);
      mul_GF(Ry, *Rz, tmp[7]); // Ry = Rz G (D - a C)
      mul_GF(Rz, tmp[5], tmp[6]); // Rz = F G
    }
  }

  /* Convert to Montgomery */
  add_GF(&tmp[0], *Rz, *Ry);
  sub_GF(&tmp[1], *Rz, *Ry);
  mul_GF(Ry, tmp[0], *Rz); // Ry = (Rz+Ry)Rz
  mul_GF(Rz, tmp[1], *Rx); // Rz = (Rz-Ry)Rx
  mul_GF(Rx, tmp[0], *Rx); // Rx = (Rz+Ry)Rx

  clear_GF(&a); clear_GF(&d);
  clear_GF(&aPx); clear_GF(&aPy);
  clear_GF(&aQx); clear_GF(&aQy);
  clear_GF(&PQx); clear_GF(&PQy);
}




/************* ISOGENIES ******************/
typedef struct {
  GF u, r;
} iso;

typedef GF iso2;

typedef struct {
  GF p, p2;
} iso3;

void init_iso3( iso3 *iso, GF_params *parent){
 
    init_GF(&(iso->p), parent);
    init_GF(&(iso->p2), parent);
    
}

void copy_iso3( iso3 *res, iso3 iso){
    
    copy_GF(&(res->p), iso.p);
    copy_GF(&(res->p2), iso.p2);
    
}

typedef struct {
  GF Ap2;
} iso4;

/* Utility routine to compute (A+2)/4 */
void a24(GF* A24, const GF A) {
  GF_params* field = A.parent;
  GF* tmp = field->GFtmp;

  add_GF_ui(&tmp[0], A, 2);
  mpz_set_ui(field->tmp1, 4);
  mpz_invert(field->tmp2, field->tmp1, field->p);
  scalar_GF(A24, tmp[0], field->tmp2);
}

/*
  Compute an isomorphism of the montgomery curve
  sending (x,z) to (0,0).
*/
void isom_comp(iso* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B, const GF A24,
	       const GF x, const GF z) {
  GF* tmp = A.parent->GFtmp;
  
  mont_double(&tmp[1], &tmp[2], x, z, A24);
  div_GF(&tmp[0], tmp[1], tmp[2]);  // P2x = x([2]P) / z([2]P)
  neg_GF(&iso->r, tmp[0]);  // r = -P2x
  scalar_GF_si(&tmp[1], tmp[0], 3);
  add_GF(&tmp[1], tmp[1], A);  // a2 = 3 P2x + A
  mul_GF(&tmp[2], iso->r, z);
  add_GF(&tmp[3], tmp[2], x);
  div_GF(&iso->u, z, tmp[3]);  // u = z / (x - z P2x)
  mul_GF(iA, tmp[1], iso->u);  // iA = a2 u
  mul_GF(iB, B, iso->u);  // iB = B u
  a24(iA24, *iA);
}

/* Amply an isomorphism of Montgomery curves */
void isom_apply(GF* X, GF* Y, GF* unused,
		const iso iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], iso.r, z);
  add_GF(&tmp[1], x, tmp[0]);
  if (Y)
    mul_GF(Y, y, iso.u); // Y = y u
  mul_GF(X, tmp[1], iso.u); // X = (x + r z) u
}


/*
  Compute a 2-isogeny of the montgomery curve
  sending (x,z) to (1,...).
*/
void iso2_comp(iso2* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B,
	       const GF x, const GF z) {
  GF* tmp = x.parent->GFtmp;

  sub_GF(&tmp[0], x, z);
  sqr_GF(&tmp[1], tmp[0]);
  inv_GF(&tmp[0], tmp[1]);
  mul_GF(&tmp[1], tmp[0], z);
  mul_GF(iso, tmp[1], x); // iA2 = x z / (x-z)^2
  add_GF_ui(&tmp[0], A, 6);
  mul_GF(iB, B, *iso); // iB = B iA2
  mul_GF(iA, tmp[0], *iso); // iA = (A+6) iA2
  a24(iA24, *iA);
}

/* Apply a 2-isogeny of Montgomery curves */
void iso2_apply(GF* X, GF* Y, GF* Z,
		const iso2 iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;
  
  sub_GF(&tmp[3], x, z);
  sqr_GF(&tmp[4], tmp[3]);
  if (Y) {
    mul_GF(&tmp[4], x, tmp[4]); // ... X = x iA2 (x - z)^2
    sqr_GF(&tmp[0], x); // Px2 = x^2
    sqr_GF(&tmp[1], z);
    sub_GF(&tmp[2], tmp[0], tmp[1]);
    mul_GF(&tmp[1], y, tmp[2]);
    mul_GF(Y, iso, tmp[1]); // Y = iA2 y (x^2 - z^2)
    mul_GF(Z, z, tmp[0]); // Z = z x^2
  } else {
    mul_GF(Z, z, x); // Z = z x
  }
  mul_GF(X, iso, tmp[4]); // X = iA2 (x - z)^2
    
    count2iso++;
}

/*
  Compute a 3-isogeny of the montgomery curve
*/
void iso3_comp(iso3* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B,
	       const GF x, const GF z) {
  GF* tmp = x.parent->GFtmp;

  div_GF(&iso->p, x, z);             // p
  sqr_GF(&iso->p2, iso->p);          // p^2
  
  scalar_GF_si(&tmp[3], iso->p, -6); 
  add_GF(&tmp[4], tmp[3], A);
  mul_GF(&tmp[3], tmp[4], iso->p);
  add_GF_ui(&tmp[4], tmp[3], 6);     // (-6p + A)p + 6

  mul_GF(iB, B, iso->p2);      // iB = B p^2
  mul_GF(iA, tmp[4], iso->p);  // iA = ((-6p + A)p + 6)p

  a24(iA24, *iA);
}

/* Apply a 3-isogeny of Montgomery curves */
void iso3_apply(GF* X, GF* Y, GF* Z,
		const iso3 iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], z, iso.p);
  sub_GF(&tmp[1], x, tmp[0]); // h = x - p z
                              // if zero, P is in the kernel
  mul_GF(&tmp[2], x, iso.p);
  sub_GF(&tmp[0], tmp[2], z); // rh = x p - z
  sqr_GF(&tmp[3], tmp[0]);
  mul_GF(&tmp[2], x, tmp[3]); // X0 = x (x p - z)^2

  if (Y) {
    mul_GF(&tmp[2], tmp[2], tmp[1]); // X0 *= h
    mul_GF(&tmp[3], x, z);
    sub_GF_ui(&tmp[4], iso.p2, 1);
    mul_GF(&tmp[5], tmp[3], tmp[4]);
    scalar_GF_si(&tmp[3], tmp[5], -2);
    mul_GF(&tmp[4], tmp[0], tmp[1]);
    add_GF(&tmp[5], tmp[3], tmp[4]);
    mul_GF(&tmp[3], tmp[5], tmp[0]); // (rh (rh h + 2xz(1-p^2)))
    sqr_GF(&tmp[7], tmp[1]);
    mul_GF(&tmp[8], tmp[7], z);
    mul_GF(Y, y, tmp[3]); // Y = y (rh (rh h + 2xz(1-p^2)))
    mul_GF(Z, tmp[8], tmp[1]); // Z = h^2 h z
  } else {
    sqr_GF(&tmp[3], tmp[1]);
    mul_GF(Z, tmp[3], z); // Z = h^2 z
  }
  copy_GF(X, tmp[2]);
    
  count3iso++;
}



/*
  Compute a 4-isogeny of the Montgomery curve
  sending (1,...) to infinity.
*/
void iso4_comp(iso4* iso, GF* iA, GF* iB, GF* iA24,
	       const GF A, const GF B) {
  GF* tmp = A.parent->GFtmp;

  add_GF_ui(&iso->Ap2, A, 2); // Ap2 = A + 2
  sub_GF_ui(&tmp[0], A, 2);
  neg_GF(&tmp[0], tmp[0]);
  inv_GF(&tmp[2], tmp[0]); // iAm2 = 1 / (2-A)
  add_GF_ui(&tmp[0], A, 6);
  mul_GF(&tmp[1], tmp[0], tmp[2]);
  mul_GF(iB, B, tmp[2]); // iB = B iAm2
  scalar_GF_si(iA, tmp[1], -2); // iA = -2 (A+6) iAm2
  a24(iA24, *iA);
}

/* Apply a 4-isogeny of Montgomery curves */
void iso4_apply(GF* X, GF* Y, GF* Z,
		const iso4 iso,
		const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], x, z); // z1 = x z
  sub_GF(&tmp[2], x, z);
  sqr_GF(&tmp[1], tmp[2]); // x1 = (x - z)^2
  mul_GF(&tmp[2], tmp[0], iso.Ap2); // zA2 = z1 Ap2
  scalar_GF_si(&tmp[3], tmp[0], 4); // fourz = 4 z1
  add_GF(&tmp[6], tmp[1], tmp[2]);
  add_GF(&tmp[4], tmp[1], tmp[3]);
  mul_GF(&tmp[5], tmp[6], tmp[4]); // x0 = (x1+zA2)(x1+fourz)

  if (Y) {
    mul_GF(&tmp[4], x, tmp[1]); // B = x x1
    mul_GF(&tmp[5], tmp[5], tmp[4]); // x0 *= B
    sqr_GF(&tmp[6], z);
    sub_GF(&tmp[7], tmp[0], tmp[6]);
    scalar_GF_si(&tmp[6], tmp[7], 2);
    add_GF(&tmp[7], tmp[6], tmp[1]); // C = x1 + 2(z1 - z^2)
    sqr_GF(&tmp[6], tmp[1]);
    mul_GF(&tmp[8], tmp[2], tmp[3]);
    sub_GF(&tmp[6], tmp[6], tmp[8]); // D = x1^2 - zA2 fourz
    mul_GF(&tmp[8], tmp[7], tmp[6]);
    mul_GF(Y, y, tmp[8]); // Y = y C D
    sqr_GF(&tmp[6], tmp[4]);
    sub_GF_ui(&tmp[7], iso.Ap2, 4);
      
    //The code below implements neg_GF(&tmp[7], tmp[7]). For some reason just calling the function here has unpredictable behavious due to a bug in mpz_sub
      if (mpz_sgn(tmp[7].a) == 0)
          mpz_set(tmp[7].a, tmp[7].a);
      else
        mpz_sub(tmp[7].a, tmp[7].parent->p, tmp[7].a); 
                                    
      if (mpz_sgn(tmp[7].b) == 0)
          mpz_set(tmp[7].b, x.b);
      else
          mpz_sub(tmp[7].b, tmp[7].parent->p, tmp[7].b);
      
    mul_GF(&tmp[8], tmp[6], tmp[7]);
    mul_GF(Z, z, tmp[8]); // Z = z B^2 (4 - Ap2)
  } else {
    sub_GF(&tmp[4], tmp[3], tmp[2]);
    mul_GF(Z, tmp[1], tmp[4]); // Z = x1 (fourz - zA2)
  }
  copy_GF(X, tmp[5]);
}

/* Apply a 4-isogeny of Montgomery curves (hopefully thread safe since we can specify temp variables */
 void iso4_apply_t(GF* X, GF* Y, GF* Z,
                 const iso4 iso,
                 const GF x, const GF y, const GF z, GF * tmp) {
     
     mpz_t tmpA[3];
     mpz_init(tmpA[0]);
     mpz_init(tmpA[1]);
     mpz_init(tmpA[2]);
     
    // dieter=1; 
    mul_GF_t(&tmp[0], x, z, tmpA); // z1 = x z
    sub_GF_t(&tmp[2], x, z, tmpA);
    sqr_GF_t(&tmp[1], tmp[2], tmpA); // x1 = (x - z)^2  //probably the buggy line
     
  /*   print_GF(tmp[1], "tmp1");
     print_GF(*Y, "Y");
     print_GF(*Z, "Z");
     print_GF(x, "x");
     print_GF(y, "y");
     print_GF(z, "z");  */
     
    mul_GF_t(&tmp[2], tmp[0], iso.Ap2, tmpA); // zA2 = z1 Ap2
    scalar_GF_si_t(&tmp[3], tmp[0], 4, tmpA); // fourz = 4 z1
    add_GF_t(&tmp[6], tmp[1], tmp[2], tmpA);
    add_GF_t(&tmp[4], tmp[1], tmp[3], tmpA);
    mul_GF_t(&tmp[5], tmp[6], tmp[4], tmpA); // x0 = (x1+zA2)(x1+fourz)
 
    if (Y) {
        mul_GF_t(&tmp[4], x, tmp[1], tmpA); // B = x x1
        mul_GF_t(&tmp[5], tmp[5], tmp[4], tmpA); // x0 *= B
        sqr_GF_t(&tmp[6], z, tmpA);
        sub_GF_t(&tmp[7], tmp[0], tmp[6], tmpA);
        scalar_GF_si_t(&tmp[6], tmp[7], 2, tmpA);
        add_GF_t(&tmp[7], tmp[6], tmp[1], tmpA); // C = x1 + 2(z1 - z^2)
        sqr_GF_t(&tmp[6], tmp[1], tmpA);
        
      //  printf("\ndieter\n");
        
        mul_GF_t(&tmp[8], tmp[2], tmp[3], tmpA);
        sub_GF_t(&tmp[6], tmp[6], tmp[8], tmpA); // D = x1^2 - zA2 fourz
        mul_GF_t(&tmp[8], tmp[7], tmp[6], tmpA);
        mul_GF_t(Y, y, tmp[8], tmpA); // Y = y C D
        sqr_GF_t(&tmp[6], tmp[4], tmpA);
        sub_GF_ui_t(&tmp[7], iso.Ap2, 4, tmpA);
 
        //The code below implements neg_GF(&tmp[7], tmp[7]). For some reason just calling the function here has unpredictable behavious due to a bug in mpz_sub
        if (mpz_sgn(tmp[7].a) == 0)
            mpz_set(tmp[7].a, tmp[7].a);
        else
            mpz_sub(tmp[7].a, tmp[7].parent->p, tmp[7].a); 
 
        if (mpz_sgn(tmp[7].b) == 0)
            mpz_set(tmp[7].b, x.b);
        else
            mpz_sub(tmp[7].b, tmp[7].parent->p, tmp[7].b);
 
            mul_GF_t(&tmp[8], tmp[6], tmp[7], tmpA);
            mul_GF_t(Z, z, tmp[8], tmpA); // Z = z B^2 (4 - Ap2)
        } else {
            sub_GF_t(&tmp[4], tmp[3], tmp[2], tmpA);
            mul_GF_t(Z, tmp[1], tmp[4], tmpA); // Z = x1 (fourz - zA2)
          }
        copy_GF(X, tmp[5]);
      //  dieter=0;
     
    // mpz_clear(tmpA[0]); mpz_clear(tmpA[1]); mpz_clear(tmpA[2]); 
    }
 
 
/******* COMPOSITE ISOGENIES **************/

/* Implementation of a queue */

typedef struct queue_point {
    GF x, z;
    int h;
    struct queue_point *next, *prev;
};


void copy_QP(queue_point *res, queue_point qp){
   // printf("dieter00\n");
    copy_GF(&(res->x), qp.x); //printf("dieter--\n");
    copy_GF(&(res->z), qp.z); //printf("dieter++\n");
    res->h = qp.h; // printf("dieter3\n");
    res->next = qp.next; //printf("dieter4\n");
    res->prev = qp.prev; //printf("dieter5\n");
    
}

void print_QP(queue_point qp, char * str){
    
    printf("\n%s\n",str);
    print_GF(qp.x, "x");
    print_GF(qp.z, "z");
    printf("\nh: %d", qp.h);
    
}

void init_QP( queue_point *qp, GF_params *parent ){
    //printf("dieter++\n");
    init_GF(&(qp->x), parent);
    init_GF(&(qp->z), parent);
    qp->h=0; //printf("dieter2\n");
    qp->next = malloc(sizeof(queue_point));
    qp->prev = malloc(sizeof(queue_point)); //printf("dieter3\n");
}

typedef struct {
    GF *X, *Y, *Z, x, y, z;
    iso4  d;
} thread_p2;

void * apply_t(void * thr){
    
    thread_p2 *t = (thread_p2 *) thr;
    
    int i=0;
    GF tmp[GF_TMP_REGS];// = malloc(sizeof(GF)*GF_TMP_REGS);
    for (i = 0 ; i < GF_TMP_REGS ; i++)
        init_GF(&tmp[i], t->X->parent);
   
    iso4_apply_t(t->X, t->Y, t->Z, t->d, t->x, t->y, t->z, tmp);
}

#define Q_INIT(q,field) do {	     \
    q = malloc(sizeof(queue_point)); \
    if (q) {			     \
      q->next = q->prev = NULL;	     \
      init_GF(&q->x, field);	     \
      init_GF(&q->z, field);	     \
      q->h = 0;			     \
    }				     \
  } while(0)
#define Q_CLEAR(q) do {   \
      clear_GF(&q->x);	  \
      clear_GF(&q->z);	  \
      free(q);		  \
    } while(0)
#define Q_PUSH(tail,q) do {			\
    tail->next = q;				\
    q->prev = tail;				\
    tail = q;					\
  } while(0)
#define Q_POP(tail,q) do {				\
    q = tail;						\
    tail = tail->prev;					\
    if (tail) {						\
      tail->next = NULL;				\
    }							\
  } while(0)
#define Q_NEXT(q) (q->next)
#define Q_PREV(q) (q->prev)
#define Q_ISHEAD(q) (q->prev==NULL)
#define Q_ISTAIL(q) (q->next==NULL)

// These bits of code are almost identical for 1, 2, 3, 4
// isogenies, thus we "template" them.
#define APPLY_ISOG(apply,obj,lower) do {		\
    for ( tmp = tail ; tmp ; tmp = Q_PREV(tmp)) {	\
      apply(&tmp->x, NULL, &tmp->z, obj,		\
	    tmp->x, tmp->x, tmp->z);			\
      tmp->h = tmp->h - lower;				\
    }							\
    if (Px && Py && Pz){  \
      apply(Px, Py, Pz, obj, *Px, *Py, *Pz);\
     }                           \
    if (Qx && Qy && Qz)					\
      apply(Qx, Qy, Qz, obj, *Qx, *Qy, *Qz);		\
  } while (0)              
#define COMP_ISOG(comp,obj) do {			\
    Q_POP(tail, tmp);					\
    comp(&obj, A, B, A24, *A, *B, tmp->x, tmp->z);	\
    Q_CLEAR(tmp);					\
  } while (0)


union isogenies {
  struct {
    iso d1;
    iso2 d2;
    iso4 d4;
  };
  iso3 d3;
};




void *p2(){
    
}

void *p3(){
    
}
    

/* Push (Px, Py, Pz) and (Qx, Qy, Qz) through the isogeny of kernel
 generated by (Rx, Rz) using the given strategy. */
void push_through_iso(GF *A, GF *B, GF *A24,
                      const GF Rx, const GF Rz,
                      const int ell, int *strategy, int h,
                      GF *Px, GF *Py, GF *Pz,
                      GF *Qx, GF *Qy, GF *Qz, int e) {
	
    GF_params* field = A->parent;
    int split, i, first = 1, first2=1;
    union isogenies phi;
    queue_point *tail, *tmp;
  /*  count2iso = 0;
    count2mul = 0;
    count3iso = 0;
    count3mul = 0;
    */
    if (ell == 2) {
        init_GF(&phi.d1.u, field);
        init_GF(&phi.d1.r, field);
        init_GF(&phi.d2, field);
        init_GF(&phi.d4.Ap2, field);
    } else {
        init_GF(&phi.d3.p, field);
        init_GF(&phi.d3.p2, field);
    }
    
    Q_INIT(tail, field);
    copy_GF(&tail->x, Rx);
    copy_GF(&tail->z, Rz);
    tail->h = h;
    while (tail) {
        h = tail->h;
        split = strategy[h];
        // Descend to the floor
        while (h > 1) {
            Q_INIT(tmp, field);
            copy_GF(&tmp->x, tail->x);
            copy_GF(&tmp->z, tail->z);
            
            for ( i=0 ; i < h - split ; i++) {	  
                if (ell == 2){
                    mont_double(&tmp->x, &tmp->z,
                                tmp->x, tmp->z, *A24); 
                    
                   // count2mul++;
                }
                else{ 
                    mont_triple(&tmp->x, &tmp->z,
                                tmp->x, tmp->z, *A24);
                   // count3mul++;
                }
                
            }
            tmp->h = split;
            
            Q_PUSH(tail, tmp);
            h = split;
            split = strategy[h];
        }
        // For ell=2, at the first iteration, bring the
        // 2-torsion point to (0,0)
        if (ell == 2 && first) {
            first = 0;
            Q_INIT(tmp, field); // slight abuse
            mont_double(&tmp->x, &tmp->z, tail->x, tail->z, *A24);
            isom_comp(&phi.d1, A, B, A24,
                      *A, *B, *A24, tmp->x, tmp->z);
            Q_CLEAR(tmp);
            APPLY_ISOG(isom_apply, phi.d1, 0);
        }
        
        
        // Compute and apply the isogeny
        if (ell == 2) {

            COMP_ISOG(iso2_comp, phi.d2);
            APPLY_ISOG(iso2_apply, phi.d2, 1);
                
        } else {
            COMP_ISOG(iso3_comp, phi.d3);
            
           // print_GF(phi.d3.p, "phi.d3.p"); 
           // print_GF(phi.d3.p2, "phi.d3.p2"); 
            
            APPLY_ISOG(iso3_apply, phi.d3, 1);
        }
       // count3++;   
        
    }
    // For ell=2 there is still a 4-isogeny to apply
    if (ell == 2) {
        iso4_comp(&phi.d4, A, B, A24, *A, *B);
        // This works because the queue is empty
        APPLY_ISOG(iso4_apply, phi.d4, 2);

    }
    
    if (ell == 2) {
        clear_GF(&phi.d1.u);
        clear_GF(&phi.d1.r);
        clear_GF(&phi.d2);
        clear_GF(&phi.d4.Ap2);
    } else {
        clear_GF(&phi.d3.p);
        clear_GF(&phi.d3.p2);
    } 
    
}





////////////////////dieter////////////////////////////////////////////////

//threading stuff
    

/*    
void *add_ui_t(void *arithInfo){
    add_GF_ui(arithInfo->a, arithInfo->(*x), arithInfo->(*u));
}

void *sub_t(void *arithInfo){
     sub_GF(arithInfo->a, arithInfo->(*x), arithInfo->(*y));
}
    
void *sub_ui(void *arithInfo){
    sub_GF_ui(arithInfo->a, arithInfo->(*x), arithInfo->(*u));
}

void *neg(void *arithInfo){
    neg_GF(arithInfo->a,  arithInfo->(*x));
        
}
    
void *scalar(void *arithInfo){
    ;
        
}

void *scalar_si(void *arithInfo){
    ;   
        
}
    
void *mul(void *arithInfo){
    ;   
        
}
    
void *sqr(void *arithInfo){
    ;   
        
}

void *inv(void *arithInfo){
    ;   
        
}
    
void *div(void *arithInfo){
    ;   
        
}
 
 */
    
//Montgomery Curve
typedef struct {
    GF A, B, A24;
} MC;

void copy_MC(MC *res, MC curve){
    copy_GF( &(*res).A, curve.A );
    copy_GF( &(*res).B, curve.B );
    copy_GF( &(*res).A24, curve.A24 );
}

void init_MC(MC* x) {
    GF_params* parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    init_GF( &(*x).A, parent);
    init_GF( &(*x).B, parent);
    init_GF( &(*x).A24, parent);
}


void clear_MC(MC *x) {
    clear_GF(&(*x).A);
    clear_GF(&(*x).B);
    clear_GF(&(*x).A24);
}

void print_Curve(MC *x){
    print_GF( (*x).A,"A" );
    print_GF( (*x).B,"B" );
    print_GF( (*x).A24,"A24" ); 
}

void set_Curve( MC *curve, GF A, GF B, GF A24 ){
    (*curve).A = A;
    (*curve).B = B;
    (*curve).A24 = A24;
}


void j_invariant(GF *final, MC *curve){
    GF tmp1, tmp2, tmp3, denom, num;
    
    GF_params *parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    init_GF( &tmp1, parent );
    init_GF( &tmp2, parent );
    init_GF( &tmp3, parent );
    init_GF( &denom, parent );
    init_GF( &num, parent );
    sqr_GF(&tmp1, (*curve).A);
    sub_GF_ui(&denom, tmp1, 4);
    sub_GF_ui(&tmp2, tmp1, 3 );
    sqr_GF( &tmp3, tmp2);
    mul_GF( &tmp1, tmp3, tmp2);
    scalar_GF_si(&num, tmp1, 256);
    div_GF(final, num, denom);
}


//Montgomery Point
typedef struct {
    GF x,y,z;
    MC curve;
} MP;

void init_MP(MP* a) {
    GF_params *parent;
    MC *curve;
    curve = malloc(sizeof(MC));
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    
    init_GF( &(*a).x, parent );
    init_GF( &(*a).y, parent );
    init_GF( &(*a).z, parent );
    init_MC(curve);
}
    
void copy_MP( MP *res, MP P ){
    copy_GF(&(*res).x,P.x);
    copy_GF(&(*res).y,P.y);
    copy_GF(&(*res).z,P.z);  
    copy_MC(&(*res).curve,P.curve);
}
    
void set_Curve_MC( MP *res, MC curve ){
    (*res).curve.A = curve.A;
    (*res).curve.B = curve.B;
    (*res).curve.A24 = curve.A24;
}

void set_MP( MP *res, GF x, GF y, GF z, MC curve ){
    (*res).x = x;
    (*res).y = y;
    (*res).z = z;
    
    (*res).curve.A = curve.A;
    (*res).curve.B = curve.B;
    (*res).curve.A24 = curve.A24;
}


void clear_MP(MP *a) {
    clear_GF(&(*a).x);
    clear_GF(&(*a).y);
    clear_GF(&(*a).z);
}

void print_MP(MP *a, char * string){
    printf("%s: \n",string);
    print_GF( (*a).x,"x" );
    print_GF( (*a).y,"y" );
    print_GF( (*a).z,"z" );
    print_Curve( &((*a).curve) );
}


//P+Q make sure both curves are the same!!!!
void  add(MP *res, MP Q, MP P){

    GF *x, *y, *z;
    
    GF* pxtmp = P.x.parent->GFtmp;
    GF* pytmp = P.y.parent->GFtmp;
    GF* pztmp = P.y.parent->GFtmp;

    x = malloc(sizeof(GF));
    z = malloc(sizeof(GF));
    y = malloc(sizeof(GF));
    
    init_GF(x, P.x.parent);
    init_GF(z, P.x.parent);
    init_GF(y, P.x.parent);
    
    mul_GF( &pxtmp[1], Q.y, P.z );
    mul_GF( &pxtmp[2] , Q.x, P.z );
    
 //   printf("dietergeorge-1\n");
  //  print_GF(pxtmp[2],"\npx2");
    
    mul_GF( &pxtmp[3] , Q.z, P.z );

    mul_GF( &pxtmp[4], P.y ,Q.z );
    
   // printf("dietergeorge-3\n");
   // print_GF(pxtmp[4],"\npx4");
    
    
    sub_GF( &pxtmp[5], pxtmp[4], pxtmp[1]);
    sqr_GF( &pxtmp[6], pxtmp[5] );

    mul_GF( &pxtmp[4], P.x, Q.z );
    sub_GF( &pxtmp[7], pxtmp[4], pxtmp[2] );
    
    sqr_GF( &pxtmp[8], pxtmp[7] );
    
  /*  printf("dietergeorge\n");
    print_GF(pxtmp[7],"\nx");
    print_GF(pxtmp[8],"\ny");
    */
    mul_GF( &pxtmp[0], pxtmp[7], pxtmp[8] );
    mul_GF( &pytmp[2], pxtmp[8], pxtmp[2] );

    mul_GF( &pxtmp[4], Q.curve.B, pxtmp[6] );
    mul_GF( &pytmp[1], Q.curve.A, pxtmp[8] );
    sub_GF( &pytmp[8], pxtmp[4], pytmp[1] );
    mul_GF( &pztmp[1], pytmp[8], pxtmp[3] );
    sub_GF( &pxtmp[4], pytmp[1], pxtmp[0] );
    scalar_GF_si( &pytmp[1], pytmp[2], 2 );
    sub_GF( &pytmp[4], pxtmp[4], pytmp[1] ); 

    mul_GF(x, pxtmp[7], pytmp[4] );
    mul_GF(z, pxtmp[0], pxtmp[3] );

    sub_GF( &pxtmp[4], pytmp[2], pytmp[4] );
    mul_GF( &pytmp[1], pxtmp[5], pxtmp[4] );   
    mul_GF( &pytmp[6], pxtmp[0], pxtmp[1] );
    sub_GF( y, pytmp[1], pytmp[6] );
    
    (*res).x = *x;
    (*res).y = *y;
    (*res).z = *z;
    (*res).curve = P.curve;
    
}

//P-Q make sure both curves are the same!!
void  subtract(MP *res, MP P, MP Q){
    MP *point;
    point = malloc(sizeof(MP));
    init_MP(point);
    
    neg(point,Q);
    add(res,P,*point);
    free(point);

}


//Returns the negative of P
void neg(MP *res, MP P){
    GF *yN;
    
    yN = malloc(sizeof(GF));
    init_GF( yN, (P).x.parent);
    neg_GF( yN , (P).y );

    (*res).x = (P).x;
    (*res).y = *yN;
    (*res).z = (P).z;
    (*res).curve = P.curve;
    free(yN);
}


//returns a random integer from 0 to m (inclusive)
void rand_range(mpz_t *num, mpz_t m){
    
    unsigned long int bytes;
    mpz_t *tmp1, *tmp2;
    tmp1 = malloc(sizeof(mpz_t));
    tmp2 = malloc(sizeof(mpz_t));
    mpz_init(tmp1);
    mpz_init(tmp2);
    bytes = ceil( mpz_sizeinbase(m, 2)/8.0 ) + 20;
    unsigned char buf[bytes];
    RAND_bytes(buf, bytes);
    
    unsigned char buf1[bytes];
    char* buf2 = buf1;

    int k=0;
    for (k; k < bytes; k++)
        buf2 += sprintf(buf2, "%02X", buf[k]);
    
    sprintf(buf2,"\n");
    *(buf2 + 1) = '\0';
    mpz_set_str(tmp1, buf1, 16);
    mpz_add_ui(tmp2, m, 1);
    
    mpz_mod(num, tmp1, tmp2); 
    free(tmp1);free(tmp2);
}



void rand_subgroup(mpz_t *m, mpz_t *n, char * l, char * e){
  
    mpz_t *num, *_l, *_e, *tmp1, *le, *le1, *l1;
    num = malloc(sizeof(mpz_t));
    _l = malloc(sizeof(mpz_t));
    _e = malloc(sizeof(mpz_t));
    tmp1 = malloc(sizeof(mpz_t));
    le = malloc(sizeof(mpz_t));
    le1 = malloc(sizeof(mpz_t));
    l1 = malloc(sizeof(mpz_t));
     
    long int e_int = 0;
    e_int = atoi(e);
    mpz_init(le);
    mpz_init(le1);
    mpz_init(l1);
    mpz_init(tmp1);
    mpz_init(num);
    
    mpz_init_set_str(_l, l, 10);
    mpz_init_set_str(_e, e, 10);
    mpz_add_ui(l1, _l, 1);
    rand_range(num, *l1);
        
    if ( mpz_cmp_ui(*num, 1)==0 ){
        mpz_set_ui(m, 1);
        mpz_pow_ui(le, _l, e_int);
        rand_range(n, *le);
    }else{
        mpz_pow_ui(le1, _l, e_int-1);
        rand_range(tmp1, *le1);
        mpz_mul(m, tmp1, _l);
        mpz_set_ui(n, 1);
    }
    
    //for real execution, one should comment out the block below and then uncomment the above if/else block
  /*  while( mpz_cmp_ui(*num, 1)==0 ){
       rand_range(num, *l1);
        
       // mpz_set_ui(m, 1);
       // mpz_pow_ui(le, _l, e_int);
       // rand_range(n, *le);
    }
    

    mpz_pow_ui(le1, _l, e_int-1);
    rand_range(tmp1, *le1);
    mpz_mul(m, tmp1, _l);
    mpz_set_ui(n, 1);
*/
    
    
    free(num);free(_l);free(_e);free(tmp1);free(le);free(le1);free(l1);
    
}

void keygen_c_dfc1(MP *Pother2, MP *Qother2, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, MP Pother, MP Qother, MP QP, MP PQ, int e){

    
    GF *Px, *Py, *Pz, *Qx, *Qy, *Qz, *Rx, *Rz, *tmp, *A, *B, *A24;
   // MP *QP, *PQ; 
    MC *E;
  //  struct timeval tv1, tv2, diff;

    Rx = malloc(sizeof(GF));
    init_GF(Rx,P.x.parent);
    Rz = malloc(sizeof(GF));
    init_GF(Rz,P.x.parent);
    tmp = malloc(sizeof(GF));
    init_GF(tmp,P.x.parent);
    Px = malloc(sizeof(GF));
    init_GF(Px,P.x.parent);
    Py = malloc(sizeof(GF));
    init_GF(Py,P.x.parent);
    Pz = malloc(sizeof(GF));
    init_GF(Pz,P.x.parent);
    Qx = malloc(sizeof(GF));
    init_GF(Qx,P.x.parent);
    Qy = malloc(sizeof(GF));
    init_GF(Qy,P.x.parent);
    Qz = malloc(sizeof(GF));
    init_GF(Qz,P.x.parent);

    if(mpz_cmp_ui(M, 1)==0){
        mont_3ladder(Rx, Rz, N, P.x, P.z, Q.x, Q.z, (QP).x, (QP).z, P.curve.A24); 
        
    }else if(mpz_cmp_ui(N, 1)==0){
        mont_3ladder(Rx, Rz, M, Q.x, Q.z, P.x, P.z, (PQ).x, (PQ).z, P.curve.A24);
    }else{
        
        shamir(Rx, tmp, Rz, P.curve.A, P.curve.B, P.x, (GF)P.y, P.z, Q.x, (GF)Q.y, Q.z, M, N);
        
    }

    E = malloc(sizeof(MC));
    init_MC(E);
    
    A = malloc(sizeof(GF));
    init_GF(A,P.x.parent);
    copy_GF(A,P.curve.A);

    B = malloc(sizeof(GF));
    init_GF(B,P.x.parent);
    copy_GF(B,P.curve.B);
    
    A24 = malloc(sizeof(GF));
    init_GF(A24,P.x.parent);
    copy_GF(A24,P.curve.A24);
    
    set_Curve(E, *A, *B, *A24);
    
    copy_GF(Px, Pother.x);
    copy_GF(Py, Pother.y);
    copy_GF(Pz, Pother.z);

    copy_GF(Qx, Qother.x);
    copy_GF(Qy, Qother.y);
    copy_GF(Qz, Qother.z);
    
  //  gettimeofday(&tv1,NULL);
    push_through_iso( &(*E).A, &(*E).B, &(*E).A24, *Rx, *Rz, ell, strat, len - 1, Px, Py, Pz, Qx, Qy, Qz, e);
   // gettimeofday(&tv2,NULL);
    
    //timersub(&tv2, &tv1, &diff);
    
 /*   printf("\n\ncount2iso: %d\n\n", count2iso);
    printf("\n\ncount2mul: %d\n\n", count2mul);
    printf("\n\ncount3iso: %d\n\n", count3iso);
    printf("\n\ncount3mul: %d\n\n", count3mul);
    printf ("\n\nPush_Through Time:%f \n\n", (double) diff.tv_usec/1000000 + diff.tv_sec);
  */  
    (*Pother2).x = *Px;
    (*Pother2).y = *Py;
    (*Pother2).z = *Pz;
    
    (*Qother2).x = *Qx;
    (*Qother2).y = *Qy;
    (*Qother2).z = *Qz;
    
    set_Curve_MC(Pother2, *E);
    set_Curve_MC(Qother2, *E);
    
    free(Rx);free(Rz);free(tmp);free(Px);free(Py);free(Pz);free(Qx);free(Qy);free(Qz);free(E);free(A);free(B);
    
}


void keygen_c_dfc2(MC *E, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, int e){
    
    GF *Px, *Py, *Pz, *Qx, *Qy, *Qz, *Rx, *Rz, *tmp, *A, *B, *A24;
    MP *QP, *PQ;
   // struct timeval tv1, tv2, diff;
    
    Rx = malloc(sizeof(GF));
    init_GF(Rx,P.x.parent);
    Rz = malloc(sizeof(GF));
    init_GF(Rz,P.x.parent);
    tmp = malloc(sizeof(GF));
    init_GF(tmp,P.x.parent);
    
    Px = malloc(sizeof(GF));
    init_GF(Px,P.x.parent);
    Py = malloc(sizeof(GF));
    init_GF(Py,P.x.parent);
    Pz = malloc(sizeof(GF));
    init_GF(Pz,P.x.parent);
    
    Qx = malloc(sizeof(GF));
    init_GF(Qx,P.x.parent);
    Qy = malloc(sizeof(GF));
    init_GF(Qy,P.x.parent);
    Qz = malloc(sizeof(GF));
    init_GF(Qz,P.x.parent);

    copy_GF(Rx,P.x);
   // set_GF(Rx, "0", "0");
    copy_GF(Rz,P.x);
   // set_GF(Rz, "0", "0");
    copy_GF(tmp,P.x);
   // set_GF(tmp, "0", "0");
    
    
    if(mpz_cmp_ui(M, 1)==0){
        QP = malloc(sizeof(MP));
        init_MP(QP);
        subtract(QP,Q,P);
        mont_3ladder(Rx, Rz, N, P.x, P.z, Q.x, Q.z, (*QP).x, (*QP).z, P.curve.A24); 
        free(QP);

         
    }else if(mpz_cmp_ui(N, 1)==0){
        PQ = malloc(sizeof(MP));
        init_MP(PQ);
        subtract(PQ,P,Q);
        mont_3ladder(Rx, Rz, M, Q.x, Q.z, P.x, P.z, (*PQ).x, (*PQ).z, P.curve.A24);
        free(PQ);
       
    }else{
            
        shamir(Rx, tmp, Rz, P.curve.A, P.curve.B, P.x, (GF)P.y, P.z, Q.x, (GF)Q.y, Q.z, M, N);

    }
    
    A = malloc(sizeof(GF));
    init_GF(A,P.x.parent);
    copy_GF(A,P.curve.A);
    
    B = malloc(sizeof(GF));
    init_GF(B,P.x.parent);
    copy_GF(B,P.curve.B);
    
    A24 = malloc(sizeof(GF));
    init_GF(A24,P.x.parent);
    copy_GF(A24,P.curve.A24);
    
    set_Curve(E, *A, *B, *A24);
    
 //   gettimeofday(&tv1,NULL);
    push_through_iso( &(*E).A, &(*E).B, &(*E).A24, *Rx, *Rz, ell, strat, len - 1, NULL, NULL, NULL, NULL, NULL, NULL, e);
 //   gettimeofday(&tv2,NULL);
    
  //  timersub(&tv2, &tv1, &diff);
    
  /*  printf("\n\ncount2iso: %d\n\n", count2iso);
    printf("\n\ncount2mul: %d\n\n", count2mul);
    printf("\n\ncount3iso: %d\n\n", count3iso);
    printf("\n\ncount3mul: %d\n\n", count3mul);
    printf ("\n\nPush_Through Time:%f \n\n", (double) diff.tv_usec/1000000 + diff.tv_sec);
   */ 
    free(Rx);free(Rz);free(tmp);free(Px);free(Py);free(Pz);free(Qx);free(Qy);free(Qz);
    
    
}
   // returns 0 if shared keys are different, 1 if they are the same
double ss_isogeny_exchange_dfc(double *time, char * eA, char * eB, char * lA_str, char * lB_str, int *strA, int lenA, int *strB, int lenB, MP *PA, MP *QA, MP *PB, MP *QB){
    mpz_t *mA, *nA, *mB, *nB;
    struct timeval tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8, diff1, diff2, diff3, diff4, total;
    int lA, lB;
    int good=0;
    
  
    mA = malloc(sizeof(mpz_t));
    nA = malloc(sizeof(mpz_t));
    mB = malloc(sizeof(mpz_t));
    nB = malloc(sizeof(mpz_t));

    mpz_init2(mA, MPZ_MEMORY*64);
    mpz_init2(nA, MPZ_MEMORY*64);
    mpz_init2(mB, MPZ_MEMORY*64);
    mpz_init2(nB, MPZ_MEMORY*64);
    
   
    printf("\n\n----------------Randomly Generating Secret Keys\n");
    rand_subgroup(mA,nA,lA_str,eA);
    rand_subgroup(mB,nB,lB_str,eB);
    
    lA = atoi(lA_str);
    lB = atoi(lB_str);
    
/*    mpz_set_str(mA, "1133909468128",10);
    mpz_set_str(nA, "1",10);
    mpz_set_str(mB, "1",10);
    mpz_set_str(nB, "262874198694843566480101580581059992241",10);
  */ 
    
   /*  mpz_set_str(mA, "1",10);
     mpz_set_str(nA, "2087036195557",10);
     mpz_set_str(mB, "157156508598187157953832503245610341564",10);
     mpz_set_str(nB, "1",10);
     */
    
  /*  mpz_set_str(mA, "1",10);
    mpz_set_str(nA, "1231352653461",10);
    mpz_set_str(mB, "97338415178644876851956561975635670025",10);
    mpz_set_str(nB, "1",10);
    */
    
   /* mpz_set_str(mA, "1",10);
    mpz_set_str(nA, "699314364328",10);
    mpz_set_str(mB, "325406247161099601785004392803197906816",10);
    mpz_set_str(nB, "1",10);
    */
    
  /*   mpz_set_str(mA, "534649502542",10);
     mpz_set_str(nA, "1",10);
     mpz_set_str(mB, "165066014652889966740329010995050040775",10);
     mpz_set_str(nB, "1",10);
    */ 
    
 /*    mpz_set_str(mA, "912276545892",10);
     mpz_set_str(nA, "1",10);
     mpz_set_str(mB, "1",10);
     mpz_set_str(nB, "383405685953065612689085621797716635602",10);
   */   
    
 /*    mpz_set_str(mA, "1",10);
     mpz_set_str(nA, "73466438175038822138008022909149526646583765505717211199948202029729505640771777877098509272247553096270752136502873",10);
     mpz_set_str(mB, "18785263630172217833641411052393420167565262711490998494449247492328383488174436481302973056280518911742701981486503",10);
     mpz_set_str(nB, "1",10);
   */
   
 /*   mpz_set_str(mA, "1",10);
    mpz_set_str(nA, "20084264679067769000915020459705611909150018420813500797659487000554224441026014298022901782968064536732833647117539",10);
    mpz_set_str(mB, "3711594234558361667389254214433445389228295117387213708461816094530892999965380041842556219662218345205421524153178",10);
    mpz_set_str(nB, "1",10);
 */
    
    MC EsA, EsB;
    MP * points1;
    MP * points2;
    MP phiPB, phiQB, phiPA, phiQA, QPB, PQB, QPA, PQA ;
    
    init_MP(&phiPB);
    init_MP(&phiQB);
    init_MP(&phiPA);
    init_MP(&phiQA);
    
    init_MC(&EsA);
    init_MC(&EsB);
    
    init_MP(&QPB);
    init_MP(&PQB);
    init_MP(&QPA);
    init_MP(&PQA);
    
    t = malloc(sizeof(mt_info));
    Q_INIT(t->qp2, (PA->x).parent);
    Q_INIT(t->qp3, (PA->x).parent);
    
    subtract(&QPB, *QB, *PB);
    subtract(&PQB, *PB, *QB);
    subtract(&QPA, *QA, *PA);
    subtract(&PQA, *PA, *QA);
    
    int eA_num = atoi(eA);
    int eB_num = atoi(eB);

    printf("----------------Generating Alice's Public Data\n");
    gettimeofday(&tv1,NULL);
    keygen_c_dfc1(&phiPB, &phiQB, *mA, *nA, lA, strA, lenA, *PA, *QA, *PB, *QB, QPA, PQA, eA_num);
    gettimeofday(&tv2,NULL);
    printf("----------------Generating Bob's Public Data\n");
    gettimeofday(&tv3,NULL);
    keygen_c_dfc1(&phiPA, &phiQA, *mB, *nB, lB, strB, lenB, *PB, *QB, *PA, *QA, QPB, PQB, eB_num);
    gettimeofday(&tv4,NULL);
    printf("----------------Computing Shared Key on Alice's Side\n");
    gettimeofday(&tv5,NULL);
    keygen_c_dfc2(&EsA, *mA, *nA, lA, strA, lenA, phiPA, phiQA, eA_num);
    gettimeofday(&tv6,NULL);
    printf("----------------Computing Shared Key on Bob's Side\n");
    gettimeofday(&tv7,NULL);
    keygen_c_dfc2(&EsB, *mB, *nB, lB, strB, lenB, phiPB, phiQB, eB_num);
    gettimeofday(&tv8,NULL);
    
    GF EsA_j, EsB_j;
    GF_params *parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    init_GF( &EsA_j, parent );
    init_GF( &EsB_j, parent );
    
    j_invariant(&EsA_j, &EsA);
    j_invariant(&EsB_j, &EsB);
    printf("\n\n");
    
    if ( equals(&EsA_j, &EsB_j)!=1 ){
        gmp_printf("ERROR: the shared keys don't match! Here's the secret keys:\n\tmA = %Zd\n\tnA = %Zd\n\tmB = %Zd\n\tnB = %Zd\n",*mA,*nA,*mB,*nB);
        
        printf("*************EsA*************\n");
        print_Curve(&EsA);  
        print_GF(EsA_j,"EsA J-Invariant");
        printf("\n*************EsB*************\n");
        print_Curve(&EsB); 
        print_GF(EsB_j,"EsB J-Invariant");
    }else{
        
        printf("*************Shared Key(s)*************\n\n");
        printf("*************EsA*************\n");
        print_Curve(&EsA);  
        print_GF(EsA_j,"EsA J-Invariant");
        printf("\n*************EsB*************\n");
        print_Curve(&EsB); 
        print_GF(EsB_j,"EsB J-Invariant");
        printf("\n*************Secret Keys*************\n");
        gmp_printf(" mA: %Zd\n nA: %Zd\n mB: %Zd\n nB: %Zd\n",*mA,*nA,*mB,*nB);
        
        printf("\n*************J-Invariants Match*************\n\n");
        good = 1;
    }
   
    timersub(&tv2, &tv1, &diff1);
    timersub(&tv4, &tv3, &diff2);
    timersub(&tv6, &tv5, &diff3);
    timersub(&tv8, &tv7, &diff4);
    
    timeradd(&diff1, &diff2, &total);
    timeradd(&diff3, &total, &total);
    timeradd(&diff4, &total, &total);
    
    printf("\n\nMontgomery Ladders Total: %d\n\n",  count1);
    
    printf ("Keygen 1 Real Time (sec) :%f \n", (double) diff1.tv_usec/1000000 + diff1.tv_sec);
    printf ("Keygen 2 Real Time (sec) :%f \n", (double) diff2.tv_usec/1000000 + diff2.tv_sec);
    printf ("Keygen 3 Real Time (sec) :%f \n", (double) diff3.tv_usec/1000000 + diff3.tv_sec);
    printf ("Keygen 4 Real Time (sec) :%f \n", (double) diff4.tv_usec/1000000 + diff3.tv_sec);
    printf ("Total Time (sec) :%f \n", (double) total.tv_usec/1000000.0 + total.tv_sec);
    
    *time = total.tv_usec/1000000.0 + total.tv_sec;
    return good;
    
    free(mA);free(mB);free(nA);free(nB);free(parent);free(t);
}

//reads public parameters for use with ss_isogeny_exchange_dfc in from file. Note that this file must match the format and naming conventions
//of the one generated by ss_isogeny_gen_file(). k and m are Barrett precomputations
void params_from_file( char *kF1, int *kB1,  char *mB1, char * p, char *eA, char *eB, char *lA, char *lB, int *strA, int *lenA, int *strB, int *lenB, MP *PA, MP *QA, MP *PB, MP *QB, char * file ){

    FILE *fr; 
    fr = fopen (file, "rt");
    MC curve;
    init_MC(&curve);
    int lineMax = 10000;
    char line[lineMax];

    //taking in p,lA, lB, eA, eB
    fgets(line, lineMax, fr);
    sscanf (line, "%s", p);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", lA);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", lB);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", eA);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", eB);
    
    //taking in the curve
    GF A, B, A24;
    GF_params *parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent, p);
   
    char *a;
    char *b;
    a = malloc(sizeof(char)*lineMax);
    b = malloc(sizeof(char)*lineMax);

    init_GF(&A, parent);
    init_GF(&B, parent);
    init_GF(&A24, parent);
 
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&A, a, b);
   
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&B, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
   
    set_GF(&A24, a, b);

    set_Curve( &curve, A, B, A24);
  
    //taking in PA
    GF x1,y1,z1;
    
    init_GF(&x1, parent);
    init_GF(&y1, parent);
    init_GF(&z1, parent);
   
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x1, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y1, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z1, a, b);
    
    set_MP( PA, x1, y1, z1, curve );
    
    //taking in QA
    GF x2,y2,z2;
    
    init_GF(&x2, parent);
    init_GF(&y2, parent);
    init_GF(&z2, parent);
    
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x2, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y2, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z2, a, b);
    
    set_MP( QA, x2, y2, z2, curve );
    
    //taking in PB
    GF x3,y3,z3;
    
    init_GF(&x3, parent);
    init_GF(&y3, parent);
    init_GF(&z3, parent);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x3, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y3, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z3, a, b);
    
    set_MP( PB, x3, y3, z3, curve );
    
    //taking in QB
    GF x4,y4,z4;
    
    init_GF(&x4, parent);
    init_GF(&y4, parent);
    init_GF(&z4, parent);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x4, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y4, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z4, a, b);
    
    set_MP( QB, x4, y4, z4, curve );
    
    //taking in strA
    fgets(line, lineMax, fr);
    sscanf (line, "%d", lenA);
    
    int k=0;
    for(; k<*lenA; k++){
        fgets(line, lineMax, fr);
        sscanf (line, "%d", &strA[k]);
    }
    
    //taking in strB
    fgets(line, lineMax, fr);
    sscanf (line, "%d", lenB);
    
    int r=0;
    for(; r<*lenB; r++){
        fgets(line, lineMax, fr);
        sscanf (line, "%d", &strB[r]);
    }
    
    //Barrett precomputations
    fgets(line, lineMax, fr);
    sscanf (line, "%d", kB1);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", mB1);
                  
    fgets(line, lineMax, fr);
    sscanf (line, "%s", kF1);
    
    free(a);free(b);
}

//1st argument specifies file with parameters, second # of times to run the key exchange. If second argument is null, only 1 iteration is performed
int main(int argc, char *argv[]) {

 //below is the code used to run the key exchange. It's commented out to test the x86 assembly optimizations. comment those out and uncomment the code
 //below to run the key exchange   
  int iterations;
  MPZ_MEMORY = 37; //sets amount of memory (in limbs) to allocate to mpz objects
   // mpz_init(&prime);
    //mpz_set_str(prime, "bb41c3ca78b2286ca800001", 16);
    
     mpz_init(&prime);
     mpz_set_str(prime, "9161191555982008052298538759697325872858383005444503030763917191888120427263653604739574602371851919945332710234806205297475768266460658683484318356498713773944703702864057467786913144364234277796785269800198817400814717913480036351", 10);
    
    if( argc < 2 ){
        printf("ERROR: Must specify file name in 1st command line agrument.\n");
        
    }else{
        printf("File where Public parameters are: '%s'\n", argv[1]);
    }
    
    if( !argv[2] )
        iterations = 1;
    else{
        iterations = atoi(argv[2]);
        printf("Number of Iterations: %s\n", argv[2]);
    }    
    
    int MAX_LENGTH = 10000;
    int *strA, *strA_t, *strB, *strB_t; 
    int lenA=0;
    int lenB=0;
    kB = malloc(sizeof(int));
    MP *PA;
    MP *QA;
    MP *PB;
    MP *QB;
    char *p, *eA, *eB, *lA, *lB, *mB_str, *kF_str;
    
    p = malloc(sizeof(char)*MAX_LENGTH);
    eA = malloc(sizeof(char)*MAX_LENGTH);
    eB = malloc(sizeof(char)*MAX_LENGTH);
    lA = malloc(sizeof(char)*MAX_LENGTH);
    lB = malloc(sizeof(char)*MAX_LENGTH);
    PA = malloc(sizeof(MP));
    QA = malloc(sizeof(MP));
    PB = malloc(sizeof(MP));
    QB = malloc(sizeof(MP));
    strA_t = malloc(MAX_LENGTH*sizeof(int));
    strB_t = malloc(MAX_LENGTH*sizeof(int));
    
    //initializing a couple of global variables
    mB_str = malloc(sizeof(char)*MAX_LENGTH);
    kF_str = malloc(sizeof(char)*MAX_LENGTH);
    
    mpz_init2(tmp, MPZ_MEMORY*64);
    mpz_init2(mB, MPZ_MEMORY*64);
    mpz_init2(kF, MPZ_MEMORY*64);

    
    params_from_file( kF_str, kB, mB_str, p, eA, eB, lA, lB, strA_t, &lenA, strB_t, &lenB, PA, QA, PB, QB, argv[1] );

    strA = malloc(lenA*sizeof(int));
    strB = malloc(lenB*sizeof(int));
    mpz_set_str(mB, mB_str,10);
    mpz_set_str(kF, kF_str,10);
    
    //Putting strategies into arrays of the exact required length
    int r=0;
    for(r; r<lenA; r++)
        strA[r] = strA_t[r];
        
    int k=0;
    for(k; k<lenB; k++)
        strB[k] = strB_t[k]; 

    //run the key exchange 'iterations' number of times and calculate the arverage time taken
    int i=0;
    int good = 0;
    double totalTime;
    double avgTime;
    int errors=0;
   
    //testing Barrett reduction method
  /*  mpz_t a, res, tmp;
    mpz_init(tmp);
    mpz_init(res);
    mpz_init(a);
    mpz_set_str(a,"13975575381836014193582095525003842844356469803607471",10);
    
    printf("*****BEFORE BARRETT****");
    gmp_printf("\na: %Zd\n", a);
    gmp_printf("\nres: %Zd\n\n", res);
    
    barrett(res, a);
    
    printf("*****AFTER BARRETT****");
    gmp_printf("\na: %Zd\n", a);
    gmp_printf("\nres: %Zd\n\n", res);
   */
 
/*   //below is code for basic testing for the x86 assembly optimizations
   mpz_t a, b1, c, pa, res_ass, res_mpz;
    struct timeval tv1, tv2, tv3, tv4, diff1, diff2, tv5, tv6, diff3;
    
    mpz_init2(res_ass, 64*MPZ_MEMORY);
    mpz_init2(res_mpz, 64*MPZ_MEMORY);
    mpz_init(&a);
    mpz_init(&b1);
    mpz_init(&c);
    mpz_init(&pa);
    
   // _mpz_realloc(res_ass,6);
   // mpz_set_str(pa, "1b5e4d5ffffffffffffffffffffffff", 16);
  
 //   printf("\n\nSIZE: %d\n\n",  mpz_sizeinbase(b,2));
  // gmp_printf("\n%Zd\n",a);
    
 //   mpz_set_str(pa, "91611915559820080522985387596f7325872858383005444503030763917191888120427263653604739574602371851919945332710234806205297475768266460658683484318356498713773944703702864057467786913144364234277796785269800198817400814717913480036351", 10);

  //  mpz_set_str(c, "fffffffffffffffff", 16);
    
     GF x;
     GF_params *parent = malloc(sizeof(GF_params));
     setup_GF(parent, "9161191555982008052298538759697325872858383005444503030763917191888120427263653604739574602371851919945332710234806205297475768266460658683484318356498713773944703702864057467786913144364234277796785269800198817400814717913480036351"); 
     init_GF(&x,parent);

    
   // int i=0;
    int iterations1 = 1;
    
    mpz_set_str(b1, "ffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    mpz_set_str(a,  "ffffffffffffffffffffffffffffffffffffffffffffffff", 16);
    gmp_printf("\nbefore ass a: %Zx\n",a);
    printf("before ass a size: %d\n",a->_mp_size);
     gmp_printf("\nbefore ass b: %Zx\n",b1);
     printf("before ass b size: %d\n",b1->_mp_size);
        
    gettimeofday(&tv1,NULL);
    
    for(i=0; i<iterations1; i++){ 
        mpz_mul_x86_0(res_ass, a, b1);   // 3 limb multiplication with SIMD instructions...small test case example
     //  mpz_addm_x86(res_ass, a, b1);
        
    }

    gettimeofday(&tv2,NULL);

  //  mpz_set_str(b1, "111111111111111111111111111111111111111111111111", 16);
  //  mpz_set_str(a,  "111111111111111111111111111111111111111111111111", 16);    
    
    
    gmp_printf("\nbefore mpz a: %Zx",a);
    printf("\nbefore mpz a size: %d\n",a->_mp_size);
    gmp_printf("\nbefore mpz b: %Zx",b1);
    printf("\nbefore mpz b size: %d\n",b1->_mp_size);
    
   
    gettimeofday(&tv3,NULL);
    
    for(i=0; i<iterations1; i++){
        mpz_mul(res_mpz, a, b1);
        //gmp_printf("\n\nafter add before mpd (mpz): %Zx\n\n", x.parent->tmp1);
        //mpz_mod(res_mpz, x.parent->tmp1, x.parent->p);   
    }
    
        gettimeofday(&tv4,NULL);
    
    mpz_t tmp1, tmp2;
    
    mpz_init(tmp1);
    mpz_init(tmp2);

    
    

    
    timersub(&tv2, &tv1, &diff1);
    timersub(&tv4, &tv3, &diff2);
    timersub(&tv6, &tv5, &diff3);
    
    printf ("\nX86-64 Assembly Version Real Time  for %d iterations(sec) :%f \n", iterations1, (double) diff1.tv_usec/(1000000) + diff1.tv_sec);
    printf ("MPZ Version Real Time  for %d iterations (sec) :%f \n", iterations1,(double) diff2.tv_usec/(1000000) + diff2.tv_sec);
  //  printf ("EXP X86-64 Assembly Version Real Time  for %d iterations(sec) :%f\n", iterations1,(double) diff3.tv_usec/(1000000) + diff3.tv_sec);
  //  gmp_printf("\n\nafter function: \n%Zd\n", c);
    
   
    gmp_printf("\n\nres_ass: %Zd\n", res_ass);
  //  gmp_printf("\nres_ass (hex): %Zx\n", res_ass);
    printf("res_ass size: %d\n", (res_ass)->_mp_size);
   // dont ask for the hex code...everything dies...
    
    gmp_printf("\nres_mpz: %Zd",res_mpz);
    printf("\nres_mpz size: %d",res_mpz->_mp_size);
    
    printf("\n\n\nass res:  %lx %lx %lx %lx %lx %lx\n\n", (res_ass->_mp_d)[5], (res_ass->_mp_d)[4], (res_ass->_mp_d)[3], (res_ass->_mp_d)[2], (res_ass->_mp_d)[1], (res_ass->_mp_d)[0]);
    printf("mpz res: %lx %lx %lx %lx %lx %lx\n\n", (res_mpz->_mp_d)[5], (res_mpz->_mp_d)[4], (res_mpz->_mp_d)[3], (res_mpz->_mp_d)[2], (res_mpz->_mp_d)[1], (res_mpz->_mp_d)[0]);
   
  //  printf("\n\n\nass res: %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res_ass->_mp_d)[37], (res_ass->_mp_d)[36], (res_ass->_mp_d)[35], (res_ass->_mp_d)[24], (res_ass->_mp_d)[34], (res_ass->_mp_d)[33], (res_ass->_mp_d)[32], (res_ass->_mp_d)[31], (res_ass->_mp_d)[30], (res_ass->_mp_d)[29], (res_ass->_mp_d)[28], (res_ass->_mp_d)[27], (res_ass->_mp_d)[26], (res_ass->_mp_d)[25], (res_ass->_mp_d)[24], (res_ass->_mp_d)[23], (res_ass->_mp_d)[22], (res_ass->_mp_d)[21], (res_ass->_mp_d)[20], (res_ass->_mp_d)[19], (res_ass->_mp_d)[18], (res_ass->_mp_d)[17], (res_ass->_mp_d)[16], (res_ass->_mp_d)[15], (res_ass->_mp_d)[14], (res_ass->_mp_d)[13], (res_ass->_mp_d)[12], (res_ass->_mp_d)[11], (res_ass->_mp_d)[10], (res_ass->_mp_d)[9], (res_ass->_mp_d)[8], (res_ass->_mp_d)[7], (res_ass->_mp_d)[6], (res_ass->_mp_d)[5], (res_ass->_mp_d)[4], (res_ass->_mp_d)[3], (res_ass->_mp_d)[2], (res_ass->_mp_d)[1], (res_ass->_mp_d)[0]);
  //  printf("mpz res: %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", (res_mpz->_mp_d)[37], (res_mpz->_mp_d)[36], (res_mpz->_mp_d)[35], (res_mpz->_mp_d)[24], (res_mpz->_mp_d)[34], (res_mpz->_mp_d)[33], (res_mpz->_mp_d)[32], (res_mpz->_mp_d)[31], (res_mpz->_mp_d)[30], (res_mpz->_mp_d)[29], (res_mpz->_mp_d)[28], (res_mpz->_mp_d)[27], (res_mpz->_mp_d)[26], (res_mpz->_mp_d)[25], (res_mpz->_mp_d)[24], (res_mpz->_mp_d)[23], (res_mpz->_mp_d)[22], (res_mpz->_mp_d)[21], (res_mpz->_mp_d)[20], (res_mpz->_mp_d)[19], (res_mpz->_mp_d)[18], (res_mpz->_mp_d)[17], (res_mpz->_mp_d)[16], (res_mpz->_mp_d)[15], (res_mpz->_mp_d)[14], (res_mpz->_mp_d)[13], (res_mpz->_mp_d)[12], (res_mpz->_mp_d)[11], (res_mpz->_mp_d)[10], (res_mpz->_mp_d)[9], (res_mpz->_mp_d)[8], (res_mpz->_mp_d)[7], (res_mpz->_mp_d)[6], (res_mpz->_mp_d)[5], (res_mpz->_mp_d)[4], (res_mpz->_mp_d)[3], (res_mpz->_mp_d)[2], (res_mpz->_mp_d)[1], (res_mpz->_mp_d)[0]);

    
  //  gmp_printf("\n\nres_ass_exp: %Zd", x.parent->tmp2);
  //  gmp_printf("\n\nres_ass_exp (hex): %Zx\n\n", x.parent->tmp2);
   
   
    if( mpz_cmp(res_ass,res_mpz)==0 ){
      //  if( mpz_cmp(res_ass, x.parent->tmp2)==0 )
           printf("\n\n***************PASSED*******************\n\n");
       // else
        //   printf("\n\nFAILED\n\n%d\n\n", mpz_cmp(res_mpz,res_ass));
    
    } else
         printf("\n\nFAILED\n\n%d\n\n", mpz_cmp(res_mpz,res_ass));
    
  /*  unsigned long int fish1 = (a->_mp_d)[0];
    unsigned long int fish2 = (a->_mp_d)[1];
    
    gmp_printf("\n%lx\n", fish1);
    gmp_printf("\n%lx\n", fish2);

    
    
    unsigned long int dieter1 = (a->_mp_d);
    unsigned long int dieter2 = malloc(sizeof(long));
    int size = sizeof(long);
    printf("\nsize: %d\n",size);
    
    asm ("movq 8(%1), %0;"
         : "=c" (dieter2)
         : "c" (dieter1)
         );
     
    printf("\n%lx\n",dieter2);
   
    
    
   
  */
    //rest of key exchange
   double times[iterations];
 
    for(i; i<iterations; i++){
        
      //  pthread_create(&threads[0], NULL, p2, NULL);
      //  pthread_create(&threads[1], NULL, p3, NULL);
        
        double *time;
        time = malloc(sizeof(double));
        good = ss_isogeny_exchange_dfc(time, eA, eB, lA, lB, strA, lenA, strB, lenB, PA, QA, PB, QB);
        totalTime += *time;
        times[i] = *time;
        
        if (!good)
            errors +=1;
        
        printf("\nErrors/Total so far: %d/%d\n", errors,i+1);
    }
    
    
    avgTime = totalTime/(iterations*1.0);
    printf ("\nAverage Time for %d iteration(s) (sec) : %f \n", iterations, avgTime);
    printf ("\Number of key exchanges NOT completed successfully : %d \n", errors);
    printf ("\Number of times mpz_mod had to be subsituted for Barrett reduction : %d \n", barrettCount);
    printf("\nWITH_FIX\n");
    printf("\n# add_GF: %d", countAdd/iterations);
    printf("\nParameters: %s\n",  argv[1]);
    
   // for(i=0; i<iterations; i++)
   //     printf("\n %f", times[i]);

    
 /*  
    arith *info = malloc(sizeof(arith));
    info->a = &((*PA).x); 
    info->x = &((*PA).y);
    info->y = &((*PA).z);
    
    add_GF(&((*QA).x), (*PA).y, (*PA).z);
    print_GF((*QA).x, "res");
    
    pthread_t threads[10];  
    int rc = pthread_create(&threads[0], NULL, add_t, (void *) info); 
    
    print_GF((*PA).x, "3tmp4");  
    printf("\n\n*************\n\n");   
    
    */
    
    //some key exchange information
 /*   printf("errorCount: %d\n", errorCount);
    printf("size: %d\n", size);
     printf("Zerosize: %d\n", zero);
    
    free(kB); free(p); free(eA); free(eB); free(lA); free(lB); free(PA); free(QA); free(PB); free(QB); free(strA_t); free(strB_t); free(mB_str); free(kF_str);
    free(strA); free(strB);
*/
    int color[4] = {1, 2, 3, 4};
    
  }

	
 
						  
			
