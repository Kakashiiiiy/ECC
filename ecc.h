#ifndef ___ECC_FUNCTIONS___
#define ___ECC_FUNCTIONS___
#include <gmp.h>
typedef struct Curve{
	mpz_t p;
	mpz_t a;
	mpz_t b;
	mpz_t x;
	mpz_t y;
	mpz_t q;

}Curve;

typedef struct Point{
	mpz_t x1,y1;
	int infinity;
}Point;

void naf(Point *p, mpz_t mal,Point *returne,Curve *cu);
void initCruve(Curve *a);
void free_curve(Curve *a);
long *calc_naf_representation(mpz_t exponent,int size);
void doubleandadd(Point *p, mpz_t mal,Point *returne,Curve *cu);
void pointoperation(Point *point1,Point *point2,Curve *cu);

#endif

