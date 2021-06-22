#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ecc.h"

#define DEBUG 0

static inline void add(mpz_t a,mpz_t b,mpz_t n)
{
	mpz_add(a,a,b);
	mpz_mod(a,a,n);
}

static inline void sub(mpz_t a,mpz_t b,mpz_t n)
{
	mpz_sub(a,a,b);
	mpz_mod(a,a,n);
}

static inline void mul(mpz_t a,mpz_t b,mpz_t n){
	mpz_mul(a,a,b);
	mpz_mod(a,a,n);
}

static inline void calc_inv(mpz_t a,mpz_t n,mpz_t inverse){
	mpz_invert(inverse,a,n);
}


static void pointaddition(Point *point1,Point *point2,Curve *cu)
{	
	mpz_t s,p1,x3,y3;
	mpz_init(s);
	if (mpz_cmp(point1->x1,point2->x1)==0) //testen ob p1/p2 invers zueinander sind
	{
		mpz_neg(s,point1->y1);
		mpz_mod(s,s,cu->p);
		if (mpz_cmp(s,point2->y1)==0){
			point1->infinity = 0;
			mpz_clear(s);
			return;
		}
		
	}
	mpz_inits(p1,x3,y3,NULL);
	mpz_sub(p1,point2->y1,point1->y1); //p1 = y2-y1
	mpz_sub(s,point2->x1,point1->x1);  //s = x2 -x1
	calc_inv(s, cu->p, s);  	//inv of s
	mul(s, p1, cu->p);        	// s * p1

	mpz_mul(p1,s,s);     		//p1 = s**2
	sub(p1,point1->x1,cu->p);	// p1 - x1
	sub(p1,point2->x1,cu->p);	// p1 - x2
	mpz_set(x3,p1);				//x3 = p1

	mpz_sub(p1,point1->x1,x3);	//x1-x3
	mul(p1,s,cu->p);			//p1 * s = p1
	mpz_sub(y3,p1,point1->y1);	//y3 = p1 -s

	mpz_mod(y3,y3,cu->p);
	mpz_mod(x3,x3,cu->p);

	mpz_set(point1->x1,x3);		  //set Values for return
	mpz_set(point1->y1,y3);
	mpz_clears(s,p1,x3,y3,NULL);  //free mem
}

static void pointdouble(Point *p,Curve *cu) {
	mpz_t s,p1,x3,y3;
	if (mpz_cmp_d(p->y1,0)==0) //Punkte der Ordnung 2 haben y=0 vgl HA 18
	{
		p->infinity = 1;
		return;
	}
	mpz_inits(s,p1,x3,y3,NULL);

	mpz_mul(p1,p->x1,p->x1);
	mpz_mul_ui(p1,p1,3);
	add(p1,cu->a,cu->p);

	mpz_mul_ui(s,p->y1,2);
	calc_inv(s, cu->p, s);
	mul(s,p1,cu->p);
	//compute of s
	mpz_pow_ui(p1,s,2);     //p1 = s**2
	sub(p1,p->x1,cu->p);	// p1 - x1
	sub(p1,p->x1,cu->p);	// p1 - x2
	mpz_set(x3,p1);			//x3 = p1

	mpz_sub(p1,p->x1,x3);	//x1-x3
	mul(p1,s,cu->p);		//p1 * s = p1
	mpz_sub(y3,p1,p->y1);	//y3 = p1 -s

	mpz_mod(y3,y3,cu->p);
	mpz_mod(x3,x3,cu->p);

	mpz_set(p->x1,x3);
	mpz_set(p->y1,y3);
	mpz_clears(s,p1,x3,y3,NULL);
}

void pointoperation(Point *point1,Point *point2,Curve *cu){
	if (point1->infinity)
	{
		if (!point2->infinity)
		{
			mpz_set(point1->x1,point2->x1);
			mpz_set(point1->y1,point2->y1);
			point1->infinity = 0;
		}
		return;
	}
	if (point2->infinity) return;
	if (mpz_cmp(point1->x1,point2->x1)==0 && mpz_cmp(point1->y1,point2->y1)==0) //falls beide x und y gleich wird gedoublt
	{	
		pointdouble(point1, cu);
		return;
	}
	else{
		pointaddition(point1, point2, cu);
	}

}

void doubleandadd(Point *p, mpz_t factor,Point *returne,Curve *cu) {
	if (DEBUG)
	{
		gmp_printf("DAA: Point:\n%ZX\n%ZX\nMal:\n%ZX\n\n",p->x1,p->y1,factor);
	}

	mpz_mod(factor,factor,cu->q);
	int range = mpz_sizeinbase(factor,2);
	Point *dummy_value = malloc(sizeof(Point));
	
	mpz_init_set(dummy_value->x1,p->x1);
	mpz_init_set(dummy_value->y1,p->y1);

	for (int i = range-2; i >=0; --i)
	{
		pointdouble(dummy_value,cu);
		if (mpz_tstbit(factor,i)){
			pointaddition(dummy_value,p,cu);
		}
	}
	mpz_set(returne->x1,dummy_value->x1);
	mpz_set(returne->y1,dummy_value->y1);
	free(dummy_value);
}

long *calc_naf_representation(mpz_t exponent,int size){
	//Hier die Naf berechnen, daran denken, dass die Naf die Länge
	//der Binärrepräsentation um 1 erhöhen kann | Array dynamisch anfordern und returnen
	mpz_t loop_var;
	mpz_init_set(loop_var,exponent);
	long * dummy = malloc((size+1) * sizeof(long));  //Naf String kann 1 Laenger sein
	
	dummy[size] = 2; //Later Naf expansion check
	
	//Berechne NAF wie im Algo im Script
	//----------------------------------------------------------------------------------
	long x;
	for (int i = 0; i <= size; ++i)
	{
		if (mpz_cmp_ui(loop_var, 0) == 0)break;
		x = mpz_fdiv_ui(loop_var,4);
		switch (x) {
		case 1:
			dummy[i] = 2 - x;
			mpz_sub_ui(loop_var,loop_var,1);
			break;
		case 3:
			dummy[i] = 2 - x;
			mpz_add_ui(loop_var,loop_var,1);
			break;
		default:
			dummy[i] = 0;
		}
		mpz_tdiv_q_ui(loop_var,loop_var,2);
	}	
	//----------------------------------------------------------------------------------
	mpz_clear(loop_var);
	return dummy;
}

void naf(Point *p, mpz_t mal,Point *returne,Curve *cu) { //p und returne MUESSEN unterschiedlich sein
	if (DEBUG)
	{
		gmp_printf("NAF: Point:\n%ZX\n%ZX\nMal:\n%ZX\n\n",p->x1,p->y1,mal);
	}
	mpz_mod(mal,mal,cu->q);   //mal mit ordnung der Kurve reduzieren.
	Point inverse_p;
	inverse_p.infinity = 0;
	mpz_init_set(inverse_p.x1,p->x1);
	mpz_init(inverse_p.y1);

	mpz_neg(inverse_p.y1,p->y1);
	mpz_mod(inverse_p.y1,inverse_p.y1,cu->p);
	
	int size = mpz_sizeinbase(mal,2);
	long * dummy = calc_naf_representation(mal,size);
	//testen ob sich binary rep. verlaengert hat durch NAF Bildung
	if (dummy[size]==2)size--;

	mpz_set(returne->x1,p->x1);
	mpz_set(returne->y1,p->y1);
	returne->infinity = 0;
	// Abhängig vom NAF Ergebnis Square mul oder INV Mul
	// -------------------------------------------------------------------------------------
	for (int i = size-1; i >= 0; --i)
	{
		pointoperation(returne, returne, cu);
		switch (dummy[i]) {
		case 1:
			pointoperation(returne,p,cu);
			break;
		case -1:
			pointoperation(returne,&inverse_p,cu);
			break;
		}
	}
	// -------------------------------------------------------------------------------------
	// Free Memory
	free(dummy);
	mpz_clear(inverse_p.x1);
	mpz_clear(inverse_p.y1);

}

void initCruve(Curve *a){
	mpz_init_set_str(a->a,"7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9",16);
	mpz_init_set_str(a->p,"A9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377",16);
	mpz_init_set_str(a->x,"8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262",16);
	mpz_init_set_str(a->y,"547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997",16);
	mpz_init_set_str(a->b,"26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6",16);
	mpz_init_set_str(a->q,"A9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7",16);
}

void free_curve(Curve *a){
	mpz_clears(a->a,a->p,a->x,a->y,a->b,a->q,NULL);
}

int test_eq(Point *p1, Point *p2){
	return !( (mpz_cmp(p1->x1,p2->x1)==0) && (mpz_cmp(p1->y1,p2->y1)==0) && (p1->infinity==p2->infinity) );
}

int main()
{
	mpz_t factor;
	Curve curve;
	int err=0,size;
	long *dummy;
	long compare_value_1_naf[] = {-1l,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	long compare_value_2_naf[] = {-1,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,-1,0,-1,0,0,0,0,0,1,};
	initCruve(&curve);
	mpz_init_set_str(factor,"1111111111111111111111111111111111111111111111",2);
	size = mpz_sizeinbase(factor,2);
	dummy = calc_naf_representation(factor,size);
	err |= memcmp(dummy, &compare_value_1_naf, sizeof(compare_value_1_naf));
	mpz_set_str(factor,"FAFAFAFAFAFAFAFAFAFAFAFAFAF",16);
	size = mpz_sizeinbase(factor,2);
	dummy = calc_naf_representation(factor,size);
	err |= memcmp(dummy, &compare_value_2_naf, sizeof(compare_value_2_naf));
	if(err){
		printf("Error NAF calculation\n");
		exit(-1);
	}else {
		printf("[+] NAF calculation\n");
	}
	Point test,return_value,compare_value;
	test.infinity = 0;
	compare_value.infinity = 0;
	return_value.infinity = 0;
	mpz_init_set_str(test.x1,"8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262",16);
	mpz_init_set_str(test.y1,"547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997",16);
	mpz_init_set_str(compare_value.x1,"609638E7A3679D04AD149F698E58731598C14EA6F84631427346F912D7EFEE97",16);
	mpz_init_set_str(compare_value.y1,"4AF0802EAE475811C3686A89694D1631E14F375E811EE1C70A5D5DF647E5D31F",16);
	mpz_inits(return_value.x1,return_value.y1,NULL);
	naf(&test, factor, &return_value, &curve);
	err = test_eq(&return_value, &compare_value);
	if(err){
		printf("Error Double and Add using NAF representation\n");
		exit(-1);
	}else {
		printf("[+] NAF Point Calc\n");
	}
	doubleandadd(&test, factor, &return_value, &curve);
	err = test_eq(&return_value, &compare_value);
	if(err){
		printf("Error Double and Add \n");
		exit(-1);
	}else {
		printf("[+] Double and Add\n");
	}
	gmp_randstate_t state;
	gmp_randinit_default(state);
	for (int i = 0; i < 1000; ++i)
	{
		mpz_urandomb(factor,state,500);
		naf(&test, factor, &compare_value, &curve);
		doubleandadd(&test, factor, &return_value, &curve);
		err = test_eq(&return_value, &compare_value);
		if (err)
		{
			printf("NAF and Double and Add behave differently at least one is false\n");
			gmp_printf("DAA Result:\n%ZX\n%ZX\n",return_value.x1,return_value.y1);
			puts("\n-------------------------------------------\n");
			gmp_printf("NAF Result:\n%ZX\n%ZX\n",compare_value.x1,compare_value.y1);
			exit(-1);
		}
	}
	printf("---------------[+]All Test passed------------------\n");




	mpz_clears(test.x1,test.y1,NULL);
	free_curve(&curve);
	return 0;
}
