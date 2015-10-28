/*
 *************************************************************************
 *
 *		R E L E A S E D   I N   P U B L I C   D O M A I N 
 *
 * @(#)       Title : vnkunx!/usr/coert/projects/msdos/pol5.c
 *
 * @(#)     Version : 2.011		Creation date : 31-NOV-1989
 *
 * @(#)    Abstract : Calculates roots for 2nd, 3rd, 4th, or 5th degree
 * @(#)               polynominals.  An approximation is used only for
 * @(#)               the fifth degree.
 *
 *	       Lint :	flags:	-lm
 *		       errors:	none
 *		     warnings:	pol5.c(314): y1 redefinition hides earlier one
 *
 *	    Compile :	flags:	-O -lm
 *		       errors:	none
 *                   warnings:	none
 *
 *	       Link :	flags:	-o pol5
 *		       errors:	none
 *		     warnings:	none
 *
 *
 *          Remarks : SCO XENIX SysV release 2.2.3 native C compiler used.
 *                    Source should be compatible with MS Quick C compiler v1.0.
 *
 *          Authors : Michiel Niemeijer, the Netherlands.     NMR
 *                    Coert Vonk, the Netherlands.            VNK
 *
 *
 *    Modifications :
 *	0.010(HP41)  - HP41 Original (beta release!)      NMR 16-NOV-1986
 *	1.000(HP41)  - When solving P4, the largest real
 *			root should be eliminated to get
 *			and third degree equation         NMR 26-JAN-1986
 *	1.100(HP41)  - Added i/o handling                 VNK 02-FEB-1986
 *	2.010(MSDOS) - MSDOS Original (beta release!)     VNK 29-OCT-1989
 *	2.011(XENIX) - XENIX Original (beta release!)     VNK 21-NOV-1989
 *
 ************************************************************************
 */

/* You may distribute this code and its derivatives so long as you do not
   charge for them; you may charge a fee sufficient to cover distribution
   costs.  You may not include it in a product to be sold.  You may not use
   it as an inducement to buy.

   You may not modify this copyright notice; you may only add to the
   modification log that precedes this copyright notice.

   You may modify the code.  If you modify the code and distribute it, you
   must update the modification log to say who you are and when you made
   the change. */


#include <stdio.h>

#ifndef _COMPLEX_DEFINED
struct complex {
	double re,im;     /* real and imaginary parts */
	} ;
#define _COMPLEX_DEFINED
#endif

#include <math.h>

/* A kludge because of the missing LINELEN, MAX, and M_PI in the includes. */

#ifndef LINELEN
#define LINELEN 80
#endif

#ifndef M_PI
#define M_PI 3.1415926
#endif

#ifndef MAX
#define MAX(a,b) ((a>b) ? a : b)
#endif

#define ITTMAX  100             /* iteration loop max */
#define EPSILON 0.00000001      /* precision needed */
#define MAXORDER 5

#define EXIT(a,b,c)	{ (void)fprintf(stderr,b,c); (void)exit(a); }

/* forward referenced functions */

void    chi();
double  power(double, double);
double  phi(double,double);
void    FindThirdRoot(double * ,struct complex *);
void    FindTwoRoots(double *, struct complex *, struct complex *);
void    FindFourRoots(double *, struct complex *, struct complex *,
		struct complex *, struct complex *);
void    FindFifthRoot(double *, struct complex *);


int
main()
{
	double  coef[MAXORDER+1];
	struct complex  root[MAXORDER+1];        /* calculated roots */
	int     order,i;
	double  tmp1,tmp2;

	chi();
	(void)printf("pol5.c: input order? "); (void)scanf("%d",&order);

	if (order>MAXORDER)
		EXIT(2,"pol5/main: order to high, MAXORDER==%d\n",MAXORDER);

	for (i=0; i<=order; i++) {
		(void)printf("pol5.c: coefficient of x^%0d? ", i);
		(void)scanf("%lf",&coef[i]);
	}

	/* rescale to highest coefficient==1 */
	if (coef[order]!=0)
		for (i=0; i<order; i++)
			coef[i] /= coef[order];

	switch (order) {
	case 5: FindFifthRoot(coef,&root[4]);
		/* eliminate this fifth root */
		tmp1 = coef[3]; coef[3] = coef[4] + root[4].re;
		tmp2 = coef[2]; coef[2] = tmp1 + root[4].re*coef[3];
		tmp1 = coef[1]; coef[1] = tmp2 + root[4].re*coef[2];
		coef[0] = tmp1 + root[4].re*coef[1];

;for later; eliminate 5th root
;  a3' = a4 + root5->re
;  a2' = a3 + root5->re * a3
;  a1' = a2 + root5->re * a2
;  a0' = a1 + root5->re * a1


// a3' = a4 + root[4].re
// a2' = a3 + root[4].re * a3
// a1' = a2 + root[4].re * a2
// a0' = a1 + root[4].re * a1

	case 4: FindFourRoots(coef,&root[0],&root[1],&root[2],&root[3]);
		break;
	case 3: FindThirdRoot(coef,&root[2]);
		/* eliminate this third root */
		coef[0] = -coef[0]/root[2].re;
		coef[1] = coef[2] + root[2].re;
		coef[2] = 0; /* for the simple minded */

// a2' = 0
// a1' = a2 + root[2].re
// a0' = -a0 / root[2].re

	case 2: FindTwoRoots(coef,&root[0],&root[1]);
		break;
	default:EXIT(3,"pol5.c: order not supported (%0d).\n",i);
	}

	for (i=0; i<order; i++)
		(void)printf("pol5.c: calculated root[%d] == (%+lf%+lfj)\n",
			i, root[i].re, root[i].im);
	return(0);
}


void
chi()
{
	(void)puts("pol5.c - Solves polynomials up to the fifth degree v2.011.");
	(void)puts("  All rights reserved by Coert Vonk, November 1989.");
}

double
power(base,pwr)
double base,pwr;
{
	if (base==0)
		return(0);
	if (base<0)
		return(-exp(pwr*log(-base)));
	else
		return(exp(pwr*log(base)));
}

double
phi(re,im)
double re,im;
{
	if (re==0)
		if (im>0)
			return(M_PI/2);
		else
			return(-M_PI/2);
	if (re>0)
		return(atan(im/re));
	else
		return(atan(im/re)+M_PI);
}


/*
 *-----------------------------------------------------------------------
 *
 *             Routine : Find third root in a Third degree polynominal.
 *
 *              Status : beta
 *
 *          References : HP-UN #6, pg. 16/17
 *                       HP-UN #7, pg. 22/23
 *
 *        Return Value : void
 *
 *           Arguments : a, &root3
 *
 *  Implicit Arguments :
 *
 *             Updates : root3
 *
 *        Side effects :
 *
 *             Remarks :
 *
 *----------------------------------------------------------------------
 */

void
FindThirdRoot(a,root3)
struct complex  *root3;
double *a;
{
	double  p,q,tmp;

	root3->im=0; /* always, because real root */
	p = -a[2]*a[2]/9 + a[1]/3;
	q = -a[2]*a[2]*a[2]/27 + a[2]*a[1]/6 - a[0]/2;

	if (q==0) {
	        root3->re=0; // should perhaps be -a[2]/3
		return;
	}

	tmp = p*p*p + q*q;
	if (tmp>0)
		root3->re = power(q+sqrt(tmp),(double)1/3) +
			    power(q-sqrt(tmp),(double)1/3) -
			    a[2]/3;
	else
		root3->re = 2 * power(sqrt(-tmp+q*q),(double)1/3) *
			    cos( phi(q,sqrt(-tmp))/3 ) - a[2]/3;
}


/*
 *-----------------------------------------------------------------------
 *
 *             Routine : Finds two roots in a Second degree polynominal.
 *
 *              Status : beta
 *
 *          References : HP-UN #6, pg. 16/17
 *                       HP-UN #7, pg. 22/23
 *
 *        Return Value : void
 *
 *           Arguments : a, &root1, &root2
 *
 *  Implicit Arguments :
 *
 *             Updates : root1, root2
 *
 *        Side effects :
 *
 *             Remarks :
 *
 *----------------------------------------------------------------------
 */

void
FindTwoRoots(a,root1,root2)
double  *a;
struct complex  *root1,*root2;
{
	double          discr;

	discr = a[1]*a[1]/4 - a[0];
	if (discr>=0) {
		root1->re = -a[1]/2 + sqrt(discr);
		root2->re = -a[1]/2 - sqrt(discr);
		root1->im = root2->im = 0;
	} else {
		root1->re = root2->re = -a[1]/2;
		root2->im = -(root1->im = sqrt(-discr) );
	}
}


/*
 *-----------------------------------------------------------------------
 *
 *             Routine : Finds four roots in a Fifth degree polynomial.
 *
 *              Status : beta
 *
 *          References : HP-UN #6, pg. 16/17
 *                       HP-UN #7, pg. 22/23
 *
 *        Return Value : void
 *
 *           Arguments : a, &root1, &root2, &root3, &root4
 *
 *  Implicit Arguments :
 *
 *             Updates : root1, root2, root3, root4
 *
 *        Side effects :
 *
 *             Remarks :
 *
 *----------------------------------------------------------------------
 */

void
FindFourRoots(a,root1,root2,root3,root4)
double  *a;
struct complex  *root1,*root2,*root3,*root4;
{
	double          b[MAXORDER+1];
	struct complex  y1,y2,y3;
	double          A,B,C,D,LargestRoot;

	b[2] = -a[2];
	b[1] = a[3]*a[1] - 4*a[0];
	b[0] = a[0]*(4*a[2]-a[3]*a[3]) - a[1]*a[1];

	FindThirdRoot(b,&y1);
	b[0] = -b[0] / y1.re;
	b[1] = b[2] + y1.re;
	FindTwoRoots(b,&y2,&y3);

	LargestRoot = MAX( y1.re, MAX(y2.re,y3.re) );

	A = a[3]/2;
	B = LargestRoot/2;
	C = sqrt(A*A-a[2]+LargestRoot);
	D = sqrt(B*B-a[0]);

	if ((A*B-a[1]/2)>=0) {
		b[1] = A+C;
		b[0] = B+D;
		FindTwoRoots(b,root1,root2);
		b[1] = A-C;
		b[0] = B-D;
		FindTwoRoots(b,root3,root4);
	} else {
		b[1] = A-C;
		b[0] = B+D;
		FindTwoRoots(b,root1,root2);
		b[1] = A+C;
		b[0] = B-D;
		FindTwoRoots(b,root3,root4);
	}
	return;
}


/*
 *-----------------------------------------------------------------------
 *
 *             Routine : Find fifth root in a Fifth degree polynomial.
 *
 *              Status : beta
 *
 *          References : HP-UN #6, pg. 16/17
 *                       HP-UN #7, pg. 22/23
 *
 *        Return Value : void
 *
 *           Arguments : a, &root5
 *
 *  Implicit Arguments :
 *
 *             Updates : root5
 *
 *        Side effects :
 *
 *             Remarks :
 *
 *----------------------------------------------------------------------
 */

void
FindFifthRoot(a,root5)
double *a;
struct complex *root5;
{
	double  m,n,o,fn,err;
	unsigned itt=0;

	m=a[0]/100; n=a[0]; o=0;
	do {
		fn = a[0] + n*(a[1] + n*(a[2] + n*(a[3] + n*(a[4] + n))));
		m = m*fn;
		if (o!=fn) m = m/(o-fn);
		o = fn;
		err = n;
		n = n + m;
		err = err - n;
		if (err<0) err=0-err; /* abs() wouldn't do the job */
	} while(err>EPSILON && itt++!=ITTMAX);
	if (itt==ITTMAX)
		EXIT(1,"pol5/FindFifthRoot, iteration divergent.\n",NULL);
	root5->im = 0;
	root5->re = n + o*fn;
}
