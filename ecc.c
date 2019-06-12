#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

//mpz_t s,p1,x3,y3;
mpz_t x0_t,x1_t,y0_t,y1_t,q_t,t,a1_t,n1_t; //global um performance durch realloc zu sparen

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
	int on_curve;
}Point;

void add(mpz_t a,mpz_t b,mpz_t n)
{
	mpz_add(a,a,b);
	mpz_mod(a,a,n);
}

void sub(mpz_t a,mpz_t b,mpz_t n)
{
	mpz_sub(a,a,b);
	mpz_mod(a,a,n);
}

void mul(mpz_t a,mpz_t b,mpz_t n){
	mpz_mul(a,a,b);
	mpz_mod(a,a,n);
}

void square_and_mul(mpz_t a,mpz_t b,mpz_t n,mpz_t returne){ //a ^ b mod n
	int range = mpz_sizeinbase(b,2)-2;
	//char * z = mpz_get_str(NULL,2,b);
	mpz_t cov;
	mpz_init(cov);
	mpz_add(cov,cov,a);
	mpz_set(returne,a);

	for (; range>=0; --range)
	{
		mul(returne,returne,n);
		if (mpz_tstbit(b,range)){ //49 ascii 1
			mul(returne,cov,n);
		}
	}
}

void inv_euler(mpz_t a,mpz_t n,mpz_t inverse){
	mpz_t h;
	mpz_init(h);
	mpz_sub_ui(h,n,2);
	square_and_mul(a, h, n,inverse);
}

void inv_euklid(mpz_t a,mpz_t n,mpz_t inverse){
	//mpz_t x0_t,x1_t,y0_t,y1_t,q_t,t,a1_t,n1_t; //global um performance durch realloc zu sparen
	mpz_mod(a,a,n);

	mpz_inits(a1_t,n1_t,x0_t,y1_t,q_t,t,NULL);
	//mpz_init(n1_t);
	//mpz_init(x0_t);
	//mpz_init(y1_t);
	//mpz_init(q_t);
	//mpz_init(t);

	mpz_init_set_ui(x1_t,1);
	mpz_init_set_ui(y0_t,1);

	mpz_set(a1_t,a);
	mpz_set(n1_t,n);

	while(mpz_cmp_ui(a1_t,0)) {
	    mpz_fdiv_q(q_t,n1_t,a1_t);
	    mpz_set(t,a1_t);
	   	mpz_mod(a1_t,n1_t,a1_t);
	    mpz_set(n1_t,t);

	    mpz_set(t,y0_t);
	    mpz_set(y0_t,y1_t);
	    mul(y1_t,q_t,n);
	    sub(y1_t,t,n);

	    mpz_set(t,x0_t);
	    mpz_set(x0_t,x1_t);
	    mul(x1_t,q_t,n);
	    mpz_sub(x1_t,t,x1_t);

	}
	mpz_set(inverse,x0_t);
}

void calc_inv(mpz_t a,mpz_t n,mpz_t inverse){
	inv_euklid(a, n, inverse);
	//inv_euler(a, n, inverse);
	//mpz_invert(inverse,a,n);
}


void pointaddition(Point *point1,Point *point2,Curve *cu)
{
	static mpz_t s,p1,x3,y3;
	mpz_inits(s,p1,x3,y3,NULL);

	mpz_sub(p1,point2->y1,point1->y1); //p1 = y2-y1
	mpz_sub(s,point2->x1,point1->x1);  //s = x2 -x1
	calc_inv(s, cu->p, s);  //inv of s
	//gmp_printf("%Zd\n",s);
	mul(s, p1, cu->p);        // s * p1
	//compute of s

	mpz_mul(p1,s,s);     //p1 = s**2
	sub(p1,point1->x1,cu->p);			// p1 - x1
	sub(p1,point2->x1,cu->p);			// p1 - x2
	mpz_set(x3,p1);			//x3 = p1

	mpz_sub(p1,point1->x1,x3);		//x1-x3
	mul(p1,s,cu->p);			//p1 * s = p1
	mpz_sub(y3,p1,point1->y1);		//y3 = p1 -s

	mpz_mod(y3,y3,cu->p);
	mpz_mod(x3,x3,cu->p);

	//gmp_printf("%Zd:%Zd\n",x3,y3);
	mpz_set(point1->x1,x3);
	mpz_set(point1->y1,y3);
	mpz_clear(s);
	mpz_clear(p1);
}

void pointdouble(Point *p,Curve *cu) {
	static mpz_t s,p1,x3,y3;
	mpz_inits(s,p1,x3,y3,NULL);

	mpz_mul(p1,p->x1,p->x1);
	mpz_mul_ui(p1,p1,3);
	add(p1,cu->a,cu->p);

	mpz_mul_ui(s,p->y1,2);
	calc_inv(s, cu->p, s);
	mul(s,p1,cu->p);
	//compute of s
	mpz_pow_ui(p1,s,2);     //p1 = s**2
	sub(p1,p->x1,cu->p);			// p1 - x1
	sub(p1,p->x1,cu->p);			// p1 - x2
	mpz_set(x3,p1);			//x3 = p1

	mpz_sub(p1,p->x1,x3);		//x1-x3
	mul(p1,s,cu->p);			//p1 * s = p1
	mpz_sub(y3,p1,p->y1);		//y3 = p1 -s

	mpz_mod(y3,y3,cu->p);
	mpz_mod(x3,x3,cu->p);

	mpz_set(p->x1,x3);
	mpz_set(p->y1,y3);
	mpz_clear(s);
	mpz_clear(p1);
	//gmp_printf("%Zd:%Zd\n",x3,y3);

}

void doubleandadd(Point *p, mpz_t mal,Point *returne,Curve *cu) {
	mpz_mod(mal,mal,cu->q);
	int range = mpz_sizeinbase(mal,2);
	//char * z = mpz_get_str(NULL,2,mal);
	Point *puff = malloc(sizeof(Point));
	mpz_init_set(puff->x1,p->x1);
	mpz_init_set(puff->y1,p->y1);

	for (int i = range-2; i >=0; --i)
	{
		pointdouble(puff,cu);
		if (mpz_tstbit(mal,i)){
			//if(!mpz_cmp(puff->x1,p->x1) && !mpz_cmp(puff->y1,p->y1)){ //man hätte "mal" auch redzizeren können
			//	pointdouble(puff,cu);
			//}
			//else{
				pointaddition(puff,p,cu);
			//}
		}
	}
	mpz_set(returne->x1,puff->x1);
	mpz_set(returne->y1,puff->y1);
}

void naf(Point *p, mpz_t mal,Point *returne,Curve *cu) {
	mpz_mod(mal,mal,cu->q); //reduzieren, um zu verhindern das p = puff
	//gmp_printf("%Zx\n",mal);
	int i = mpz_sizeinbase(mal,2)-1;
	if (i<=0){
		mpz_set(returne->x1,p->x1);
		mpz_set(returne->y1,p->y1);
		return;
	};
	
	Point puff; //= malloc(sizeof(Point));
	Point inv; //= malloc(sizeof(Point));
	mpz_init_set(puff.x1,p->x1);
	mpz_init_set(puff.y1,p->y1);
	mpz_init_set(inv.x1,p->x1);
	mpz_init(inv.y1);
	mpz_neg(inv.y1,p->y1);

	while(mpz_tstbit(mal,i)) {
	    pointdouble(&puff, cu);
	    i--;
	}
	pointaddition(&puff,&inv,cu);

	while(i >=0) {
	    pointdouble(&puff, cu);
	    if (mpz_tstbit(mal,i-1)){
	    	if (mpz_tstbit(mal,i-2))
	    	{
	    		pointaddition(&puff, p, cu);
	    		i--;
	    		while(mpz_tstbit(mal,i)) {
	    		    pointdouble(&puff, cu);
	    		    i--;
	    		}
	    		pointaddition(&puff, &inv, cu);
	    		continue;
	    	}
	    }
	    if (mpz_tstbit(mal,i))
	    {
	    	pointaddition(&puff, p, cu);
	    }
	    i--;
	}
	mpz_set(returne->x1,puff.x1);
	mpz_set(returne->y1,puff.y1);
}

void initCruve(Curve *a){
	mpz_init_set_str(a->a,"7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9",16);
	mpz_init_set_str(a->p,"A9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377",16);
	mpz_init_set_str(a->x,"8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262",16);
	mpz_init_set_str(a->y,"547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997",16);
	mpz_init_set_str(a->b,"26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6",16);
	mpz_init_set_str(a->q,"A9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7",16);
}

void ecdh(Curve * cu, Point *po,mpz_t ka,mpz_t kb) {
	mpz_init_set_str(ka,"20A5B20E076E77984380CB49173F6ED7FDED87E645747133F63888907245E5D8",16);
	mpz_init_set_str(kb,"63690612179A5742A7DB7003F0545E866CAF9DE086BF272A0E1827165381B399",16);
	mpz_t t,inv,xor1,xor2;

	mpz_inits(t,inv,xor1,xor2,NULL);

	mpz_mul(t,ka,kb);
	mpz_mod(t,t,cu->q);
	Point p1;
	mpz_inits(p1.x1,p1.y1,NULL);
	doubleandadd(po,t,&p1,cu);
	gmp_printf("%ZX\n%ZX\n",p1.x1,p1.y1);

	mpz_init_set_str(inv,"340282366920938463463374607431768211455",10); //2**128-1
	mpz_and(xor1,p1.x1,inv);
	mpz_tdiv_q_2exp(xor2,p1.x1,128); //right shift
	mpz_xor(xor1,xor1,xor2);
	gmp_printf("%Zx\n",xor1);
}

int main()
{
	mpz_t inv,ka,kb;
	mpz_init(inv);
	Curve curve;
	initCruve(&curve);


	Point test;
	mpz_init_set_str(test.x1,"8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262",16);
	mpz_init_set_str(test.y1,"547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997",16);

	//pointdouble(&test, &curve);
	//gmp_printf("%ZX\n%ZX\n",test.x1,test.y1);
	//for (int i = 0; i < 10000; ++i)
		{
			ecdh(&curve, &test, ka, kb);
		}	
	//ecdh(&curve, &test, ka, kb);
	return 0;
}