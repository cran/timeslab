/* timeslab.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b54 = 1.;
static integer c__1 = 1;
static doublereal c_b225 = 2.;
static integer c__0 = 0;
static doublereal c_b307 = 0.;
static doublereal c_b318 = -2.;
static doublereal c_b375 = -1.;

/* ********************************************************************** */
/* Subroutine */ int arfilt_(alpha, p, rvar, x, n, m, w, sigma, a, ier)
doublereal *alpha;
integer *p;
doublereal *rvar, *x;
integer *n, *m;
doublereal *w, *sigma, *a;
integer *ier;
{
    /* System generated locals */
    integer x_dim1, x_offset, w_dim1, w_offset, a_dim1, a_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k, l;
    static doublereal al, al1;
    static integer ip1;
    static doublereal sr0;


/*   Subroutine to multiply the (n x m) matrix  x by the inverse square */
/*   root of the (n x n) covariance matrix of an AR(p,alpha,rvar) process.
 */

/* ***********************************************************************
 */


/*   Find coefficients and residual variances for orders 1, ... , p: */

    /* Parameter adjustments */
    a_dim1 = *p;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --sigma;
    --alpha;
    w_dim1 = *n;
    w_offset = w_dim1 + 1;
    w -= w_offset;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    sigma[*p] = *rvar;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	a[*p + i__ * a_dim1] = alpha[i__];
    }

    if (*p > 1) {
	for (i__ = *p - 1; i__ >= 1; --i__) {
	    ip1 = i__ + 1;
	    *ier = ip1;
	    al = a[i__ + 1 + (i__ + 1) * a_dim1];
	    if (abs(al) >= (float)1.) {
		return 0;
	    }
	    al1 = 1. - al * al;
	    sigma[i__] = sigma[ip1] / al1;
	    i__1 = i__;
	    for (j = 1; j <= i__1; ++j) {
/* L15: */
		a[i__ + j * a_dim1] = (a[ip1 + j * a_dim1] - al * a[ip1 + (
			ip1 - j) * a_dim1]) / al1;
	    }

/* L20: */
	}

    }

    *ier = 1;
    if ((d__1 = a[a_dim1 + 1], abs(d__1)) >= (float)1.) {
	return 0;
    }
    sr0 = sqrt(sigma[1] / ((float)1. - a[a_dim1 + 1] * a[a_dim1 + 1]));
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	sigma[i__] = sqrt(sigma[i__]);
    }
    *ier = 0;

/*   Multiply: */

    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	w[k * w_dim1 + 1] = x[k * x_dim1 + 1] / sr0;
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
/* Computing MIN */
	    i__3 = i__ - 1;
	    l = min(i__3,*p);

	    c__ = x[i__ + k * x_dim1];
	    i__3 = l;
	    for (j = 1; j <= i__3; ++j) {
/* L110: */
		c__ += a[l + j * a_dim1] * x[i__ - j + k * x_dim1];
	    }

/* L120: */
	    w[i__ + k * w_dim1] = c__ / sigma[l];
	}
/* L100: */
    }


    return 0;
} /* arfilt_ */


/* ******************************************************************** */

/* Subroutine */ int mxpd_(x, alpha, beta, np, nq, rvar, ntf, ntl, nhf, nhl, 
	nr, st, w, a, sh, st1, st2, xp, std)
doublereal *x, *alpha, *beta;
integer *np, *nq;
doublereal *rvar;
integer *ntf, *ntl, *nhf, *nhl, *nr;
doublereal *st, *w, *a, *sh, *st1, *st2, *xp, *std;
{
    /* System generated locals */
    integer st_dim1, st_offset, w_dim1, w_offset, a_dim1, a_offset, sh_dim1, 
	    sh_offset, st1_dim1, st1_offset, st2_dim1, st2_offset, i__1, i__2,
	     i__3, i__4;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int vneg_(), mxma_();
    static integer mxpq;
    static doublereal c__;
    static integer i__, j, k, l;
    static doublereal z__;
    extern /* Subroutine */ int mnvar_();
    static doublereal c1, c2;
    static integer l1, ih, mxpqp1, it;
    static doublereal xh[100], xt[100];
    static integer iptstd;
    extern /* Subroutine */ int wilson_();
    static integer nr2, nr4;
    static doublereal xt1[100], xt2[100];
    extern /* Subroutine */ int movxyr_();
    static integer ier, nr24;
    static doublereal st11;
    static integer npr;


/* ******************************************************************** */


/*   initial conditions: */

    /* Parameter adjustments */
    --x;
    --alpha;
    --beta;
    st2_dim1 = *nr;
    st2_offset = st2_dim1 + 1;
    st2 -= st2_offset;
    st1_dim1 = *nr;
    st1_offset = st1_dim1 + 1;
    st1 -= st1_offset;
    sh_dim1 = *nr;
    sh_offset = sh_dim1 + 1;
    sh -= sh_offset;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    w_dim1 = *nr;
    w_offset = w_dim1 + 1;
    w -= w_offset;
    st_dim1 = *nr;
    st_offset = st_dim1 + 1;
    st -= st_offset;
    --xp;
    --std;

    /* Function Body */
    iptstd = 1;
    mxpq = max(*np,*nq);
    mxpqp1 = mxpq + 1;
    if (*np > 0) {
	vneg_(&alpha[1], np);
    }
    if (*nq > 0) {
	vneg_(&beta[1], nq);
    }

/*   Form Initial State Vector (xh) and Covariance Matrix (sh): */

    wilson_(&alpha[1], np, &beta[1], nq, xt, &mxpqp1, xt1, &mxpqp1, xt2, &
	    mxpq, &ier);
    if (ier == 1) {
	goto L99;
    }
    if (*np > 0) {
	vneg_(&alpha[1], np);
    }
    if (*nq > 0) {
	vneg_(&beta[1], nq);
    }


    i__1 = *nr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xh[i__ - 1] = (float)0.;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    sh[i__ + j * sh_dim1] = xt[i__ - j] * *rvar;
/* L50: */
	    sh[j + i__ * sh_dim1] = sh[i__ + j * sh_dim1];
	}
    }
    if (*nr > 1) {
	mxma_(&alpha[1], &beta[1], np, nq, nr, xt1);
	i__2 = *nr;
	for (j = 2; j <= i__2; ++j) {
	    i__1 = *nr;
	    for (k = j; k <= i__1; ++k) {
		c__ = sh[j + k * sh_dim1];
		i__3 = j - 2;
		for (l = 0; l <= i__3; ++l) {
		    c1 = (float)1.;
		    if (l > 0) {
			c1 = xt1[l - 1];
		    }
		    c2 = 1.;
		    l1 = l + k - j;
		    if (l1 > 0) {
			c2 = xt1[l1 - 1];
		    }
/* L52: */
		    c__ -= c1 * c2 * *rvar;
		}
		sh[j + k * sh_dim1] = c__;
/* L51: */
		sh[k + j * sh_dim1] = c__;
	    }
	}
    }

/*   form a and w: */

    xt1[0] = (float)1.;
    if (*nr > 1) {
	i__1 = *nr;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    c__ = (float)0.;
	    if (i__ - 1 <= *nq) {
		c__ = beta[i__ - 1];
	    }
	    if (*np > 0) {
/* Computing MIN */
		i__3 = *np, i__4 = i__ - 1;
		i__2 = min(i__3,i__4);
		for (j = 1; j <= i__2; ++j) {
/* L54: */
		    c__ -= alpha[j] * xt1[i__ - j - 1];
		}
	    }
/* L55: */
	    xt1[i__ - 1] = c__;
	}
    }
    i__1 = *nr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nr;
	for (j = 1; j <= i__2; ++j) {
/* L56: */
	    w[i__ + j * w_dim1] = *rvar * xt1[i__ - 1] * xt1[j - 1];
	}
    }
    i__2 = *nr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *nr;
	for (j = 1; j <= i__1; ++j) {
/* L57: */
	    a[i__ + j * a_dim1] = (float)0.;
	}
    }
    if (*np > 0) {
	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L58: */
	    a[*nr + (*nr - i__ + 1) * a_dim1] = -alpha[i__];
	}
    }
    if (*nr > 1) {
	i__1 = *nr - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L59: */
	    a[i__ + (i__ + 1) * a_dim1] = (float)1.;
	}
    }

/*   iterate: */

    npr = 0;
    nr2 = *nr * *nr;
    nr24 = nr2 << 2;
    nr4 = *nr << 2;


    i__1 = *ntl;
    for (it = 0; it <= i__1; ++it) {

/*  x tilda and sigma tilda: */

	mnvar_(xh, &a[a_offset], &sh[sh_offset], &w[w_offset], nr, nr, xt, &
		st[st_offset]);

/*   Forecasts: */

	if (it >= *ntf && it <= *ntl) {
	    movxyr_(xt1, xt, nr);
	    movxyr_(&st1[st1_offset], &st[st_offset], &nr2);
	    i__2 = *nhl;
	    for (ih = 1; ih <= i__2; ++ih) {
		if (ih >= *nhf && ih <= *nhl) {
		    ++npr;
		    xp[npr] = xt1[0];
		    if (iptstd == 1) {
			std[npr] = sqrt(st1[st1_dim1 + 1]);
		    }
		}
		if (ih == *nhl) {
		    goto L60;
		}


		mnvar_(xt1, &a[a_offset], &st1[st1_offset], &w[w_offset], nr, 
			nr, xt2, &st2[st2_offset]);
		movxyr_(xt1, xt2, nr);
		movxyr_(&st1[st1_offset], &st2[st2_offset], &nr2);
L60:
		;
	    }
	}
	if (it == *ntl) {
	    goto L100;
	}

/*   Update x hat and sigma hat: */

	st11 = st[st_dim1 + 1];
	i__2 = *nr;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L65: */
	    xt1[i__ - 1] = st[i__ * st_dim1 + 1] / st11;
	}

/*   xhat: */

	z__ = x[it + 1] - xt[0];
	i__2 = *nr;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L70: */
	    xh[i__ - 1] = xt[i__ - 1] + xt1[i__ - 1] * z__;
	}

/*  sigma hat: */

	i__2 = *nr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *nr;
	    for (j = 1; j <= i__3; ++j) {
/* L68: */
		sh[i__ + j * sh_dim1] = st[i__ + j * st_dim1] - xt1[i__ - 1] *
			 st[j * st_dim1 + 1];
	    }
	}
/* 	do 66 i=1,nr */
/* 	do 66 j=1,nr */
/*  66	st1(i,j)=0. */
/* 	do 67 i=1,nr */
/*    	st1(i,i)=1. */
/*  67	st1(i,1)=st1(i,1)-xt1(i) */
/* 	do 68 i=1,nr */
/* 	do 68 j=1,nr */
/* 	c=0. */
/* 		do 69 k=1,nr */
/*  69		c=c+st1(i,k)*st(k,j) */
/*  68	sh(i,j)=c */


L100:
	;
    }



L99:
    return 0;
} /* mxpd_ */



/* ******************************************************************* */

/* Subroutine */ int mnvar_(x, a, sig, w, ndim, n, y, b)
doublereal *x, *a, *sig, *w;
integer *ndim, *n;
doublereal *y, *b;
{
    /* System generated locals */
    integer a_dim1, a_offset, sig_dim1, sig_offset, b_dim1, b_offset, w_dim1, 
	    w_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k, l;
    static doublereal c1;


/*   Subroutine to calculate */

/*   y=a*x     and     b=a*sig*a' + w */

/*   for nx1 x and nxn a and sig  in Kalman Filter for ARMAPRED */

/* ******************************************************************* */



    /* Parameter adjustments */
    b_dim1 = *ndim;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    w_dim1 = *ndim;
    w_offset = w_dim1 + 1;
    w -= w_offset;
    sig_dim1 = *ndim;
    sig_offset = sig_dim1 + 1;
    sig -= sig_offset;
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ < *n) {
	    c__ = x[i__ + 1];
	    goto L20;
	}

	c__ = (float)0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L10: */
	    c__ += a[i__ + j * a_dim1] * x[j];
	}
L20:
	y[i__] = c__;
    }


    if (*n == 1) {
	goto L40;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    c__ = w[i__ + j * w_dim1] + sig[i__ + 1 + (j + 1) * sig_dim1];
	    b[i__ + j * b_dim1] = c__;
/* L30: */
	    b[j + i__ * b_dim1] = c__;
	}
    }

/*   i=n: */

L40:
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
	c__ = w[*n + j * w_dim1];
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    if (j < *n) {
		c1 = sig[k + (j + 1) * sig_dim1];
		goto L60;
	    }
	    c1 = (float)0.;
	    i__3 = *n;
	    for (l = 1; l <= i__3; ++l) {
/* L50: */
		c1 += sig[k + l * sig_dim1] * a[j + l * a_dim1];
	    }
L60:
	    c__ += a[*n + k * a_dim1] * c1;
	}
	b[*n + j * b_dim1] = c__;
/* L70: */
	b[j + *n * b_dim1] = c__;
    }


    return 0;
} /* mnvar_ */



/* ******************************************************************** */

/* Subroutine */ int vneg_(x, n)
doublereal *x;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*   Subroutine to negate a vector */

/* ********************************************************************* 
*/


    /* Parameter adjustments */
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = -x[i__];
    }
    return 0;
} /* vneg_ */



/* ******************************************************************** */

/* Subroutine */ int mxma_(alpha, beta, np, nq, n, d__)
doublereal *alpha, *beta;
integer *np, *nq, *n;
doublereal *d__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal c1;


/*   Subroutine to find the first n coefficients of the MA(infinity) */
/*   representation of an ARMA(alpha,beta,np,nq) process. */

/* ******************************************************************** */


    /* Parameter adjustments */
    --alpha;
    --beta;
    --d__;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	d__[i__] = (float)0.;
    }
    if (*nq > 0) {
	i__1 = *nq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	    d__[i__] = beta[i__];
	}
    }


    if (*np > 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c__ = d__[i__];
	    if (*np == 0) {
		goto L40;
	    }
	    i__2 = min(i__,*np);
	    for (j = 1; j <= i__2; ++j) {
		c1 = 1.;
		if (i__ - j > 0) {
		    c1 = d__[i__ - j];
		}
/* L30: */
		c__ -= alpha[j] * c1;
	    }
L40:
	    d__[i__] = c__;
	}
    }


    return 0;
} /* mxma_ */


/* *********************************************************************** */

/* Subroutine */ int movxyr_(x, y, n)
doublereal *x, *y;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/************************************************************************
**/


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = y[i__];
    }

    return 0;
} /* movxyr_ */


/**************************************************************************/

/* Subroutine */ int mtpoly_(alpha, beta, p, q, gamma)
doublereal *alpha, *beta;
integer *p, *q;
doublereal *gamma;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;


/************************************************************************
**/



    /* Parameter adjustments */
    --gamma;
    --beta;
    --alpha;

    /* Function Body */
    i__1 = *p + *q;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	gamma[i__] = (float)0.;
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	gamma[i__] += alpha[i__];
    }
    i__1 = *q;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	gamma[i__] += beta[i__];
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *q;
	for (j = 1; j <= i__2; ++j) {
/* L40: */
	    gamma[i__ + j] += alpha[i__] * beta[j];
	}
    }


    return 0;
} /* mtpoly_ */


/* ****************************************************************** */

/* Subroutine */ int arpart_(alpha, p, ier)
doublereal *alpha;
integer *p, *ier;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal temp, c__;
    static integer j, k;
    static doublereal c2;


/*   Subroutine to find partials from AR coefficients. */

/* ****************************************************************** */


    /* Parameter adjustments */
    --alpha;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/* L10: */
	alpha[j] = -alpha[j];
    }
    *ier = *p;
    if ((d__1 = alpha[*p], abs(d__1)) >= (float)1.) {
	return 0;
    }
    *ier = 0;
    if (*p == 1) {
	return 0;
    }


    for (j = *p - 1; j >= 1; --j) {
	*ier = j;

	c__ = alpha[j + 1];
	c2 = 1. - c__ * c__;

	i__1 = (j + 1) / 2;
	for (k = 1; k <= i__1; ++k) {
	    temp = (alpha[k] + c__ * alpha[j + 1 - k]) / c2;
	    alpha[j + 1 - k] = (alpha[j + 1 - k] + c__ * alpha[k]) / c2;
/* L20: */
	    alpha[k] = temp;
	}

	if ((d__1 = alpha[j], abs(d__1)) >= (float)1.) {
	    return 0;
	}
/* L30: */
    }

    *ier = 0;


    return 0;
} /* arpart_ */


/***************************************************************************/

/* Subroutine */ int arsppk_(alpha, p, rvar, start, n, f, covri, ri, a, b, 
	c__, ier, omega, se)
doublereal *alpha;
integer *p;
doublereal *rvar, *start;
integer *n;
doublereal (*f) ();
doublereal *covri, *ri, *a, *b, *c__;
integer *ier;
doublereal *omega, *se;
{
    /* System generated locals */
    integer covri_dim1, covri_offset, a_dim1, a_offset, c_dim1, c_offset, 
	    i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double atan(), sin(), cos(), sqrt();

    /* Local variables */
    extern /* Subroutine */ int macv_();
    static integer i__;
    static doublereal twopi, omega1, fm, fp;
    static integer oi, it;
    extern /* Subroutine */ int freqcv_();
    static doublereal ri0, del, fpp, var;


/************************************************************************
***/



/*   if start.eq.0 find largest relative max in f: */

    /* Parameter adjustments */
    c_dim1 = *p;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    --b;
    a_dim1 = *p;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ri;
    covri_dim1 = *p;
    covri_offset = covri_dim1 + 1;
    covri -= covri_offset;
    --alpha;

    /* Function Body */
    if (*start == (float)0.) {
	fm = 0.;
	for (i__ = 2; i__ <= 128; ++i__) {
	    i__1 = i__ - 1;
	    i__2 = i__ + 1;
	    if ((*f)(&i__1) < (*f)(&i__) && (*f)(&i__) > (*f)(&i__2)) {
		if ((*f)(&i__) > fm) {
		    *start = (doublereal) (i__ - 1) / (float)256.;
		    fm = (*f)(&i__);
		}
	    }
/* L20: */
	}
    }

/*   error return if no relative max's */

    if (*start == (float)0.) {
	*ier = 1;
	goto L99;
    }

/*   Find Peak Frequency: */

    macv_(&alpha[1], p, &c_b54, &ri[1], &ri0);
    twopi = atan((float)1.) * (float)8.;
    it = 1;
    *omega = *start;

/*   start iterations: */

L35:
    fp = 0.;
    fpp = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	oi = i__;
	fp += oi * ri[i__] * sin(twopi * oi * *omega);
/* L40: */
	fpp += oi * oi * ri[i__] * cos(twopi * oi * *omega);
    }

/*   check for 0 second derivative: */

    if (abs(fpp) < (float)1e-20) {
	*ier = 2;
	goto L99;
    }

/*   check for convergence: */

    del = fp / (fpp * twopi);
    omega1 = *omega - del;
    if ((d__1 = del / *omega, abs(d__1)) < (float)1e-5) {

/*   yes: */


/*   check for convergence to 0 or .5: */

	if (omega1 < (float)1e-5 || omega1 > (float).49999) {
	    *ier = 3;
	    goto L99;
	}
	*start = omega1;
	goto L60;
    }

/*   no: */

    ++it;
    *omega = omega1;

/*   check for max # of iterations: */

    if (it > 100) {
	*ier = 4;
	goto L99;
    }
    goto L35;

/*   Asymptotic Standard Error: */

L60:
    freqcv_(start, p, &alpha[1], n, p, &covri[covri_offset], &ri[1], &a[
	    a_offset], &c__[c_offset], &b[1], &var);
    *omega = *start;
    *se = sqrt(var);
    *ier = 0;

L99:
    return 0;
} /* arsppk_ */


/* ***************************************************************** */

/* Subroutine */ int freqcv_(freq, np, alpha, nobs, ndim, covri, ri, a, c__, 
	b, var)
doublereal *freq;
integer *np;
doublereal *alpha;
integer *nobs, *ndim;
doublereal *covri, *ri, *a, *c__, *b, *var;
{
    /* System generated locals */
    integer covri_dim1, covri_offset, a_dim1, a_offset, c_dim1, c_offset, 
	    i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double atan(), sin(), cos();

    /* Local variables */
    static integer i__, j, k, l, m, ijmiv, ivpij;
    static doublereal c1, c2;
    extern /* Subroutine */ int schur_();
    static integer npmiv;
    static doublereal twopi, cc;
    static integer ij, iv;


/*   subroutine to find the asymptotic variance of a frequency */
/*   estimator freq for an ar(np) process having coeffs alpha. */

/*   input : */
/*            freq,np,alpha(1),alpha(np) */
/*            nobs: number of observations */
/*            ndim : row dimension of various arrays in calling prog */

/*   output : */
/*            var : asymptotic variance */
/*            covri,ri */

/* ******************************************************************** */



    /* Parameter adjustments */
    --b;
    --ri;
    --alpha;
    c_dim1 = *ndim;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    covri_dim1 = *ndim;
    covri_offset = covri_dim1 + 1;
    covri -= covri_offset;

    /* Function Body */
    schur_(&alpha[1], ndim, np, &a[a_offset]);
    twopi = atan((float)1.) * (float)8.;
    i__1 = *np;
    for (iv = 1; iv <= i__1; ++iv) {
	i__2 = *np;
	for (ij = 1; ij <= i__2; ++ij) {
	    cc = 0.;
	    ivpij = iv + ij;
	    ijmiv = ij - iv;
	    if (ivpij <= *np) {
		cc += alpha[ivpij];
	    }
	    if (ijmiv == 0) {
		cc += 1.;
	    }
	    if (ijmiv > 0) {
		cc += alpha[ijmiv];
	    }
/* L10: */
	    c__[iv + ij * c_dim1] = cc;
	}
    }
    i__2 = *np;
    for (j = 1; j <= i__2; ++j) {
	i__1 = j;
	for (k = 1; k <= i__1; ++k) {
	    cc = (float)0.;
	    i__3 = *np;
	    for (l = 1; l <= i__3; ++l) {
		c1 = c__[j + l * c_dim1];
		c2 = (float)0.;
		i__4 = *np;
		for (m = 1; m <= i__4; ++m) {
/* L20: */
		    c2 += a[l + m * a_dim1] * c__[k + m * c_dim1];
		}
/* L30: */
		cc += c1 * c2;
	    }
	    covri[j + k * covri_dim1] = cc;
/* L40: */
	    covri[k + j * covri_dim1] = cc;
	}
    }
    i__1 = *np;
    for (iv = 1; iv <= i__1; ++iv) {
	cc = alpha[iv];
	if (iv == *np) {
	    goto L60;
	}
	npmiv = *np - iv;
	i__2 = npmiv;
	for (j = 1; j <= i__2; ++j) {
/* L50: */
	    cc += alpha[j] * alpha[j + iv];
	}
L60:
	ri[iv] = cc;
    }

/*   find vector b : */

    i__1 = *np;
    for (iv = 1; iv <= i__1; ++iv) {
/* L90: */
	b[iv] = (doublereal) iv * sin(iv * *freq * twopi);
    }

/*   find var : */


/*   find -h''/4pi */

    cc = 0.;
    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L95: */
	cc += (doublereal) (i__ * i__) * ri[i__] * cos(twopi * (doublereal) 
		i__ * *freq);
    }
    cc = twopi * cc;
    *var = (float)0.;
    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (k = 1; k <= i__2; ++k) {
/* L100: */
	    *var += b[j] * covri[j + k * covri_dim1] * b[k];
	}
    }
    *var /= cc * cc;
    *var /= (doublereal) (*nobs);


    return 0;
} /* freqcv_ */


/* ********************************************************** */

/* Subroutine */ int mxcsd_(alph, beta, np, nq, n, mdim, wk1, r__, ainf, sd, 
	ier)
doublereal *alph, *beta;
integer *np, *nq, *n, *mdim;
doublereal *wk1, *r__, *ainf, *sd;
integer *ier;
{
    /* System generated locals */
    integer wk1_dim1, wk1_offset, ainf_dim1, ainf_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__, j, m;
    extern /* Subroutine */ int schur_(), swpk12_();
    static integer nppnq, ii, mm;
    static doublereal on;
    extern /* Subroutine */ int mxcvin_();


/*   subroutine to find the standard deviation of the estimated */
/*   coefficients (based on n observations) of an arma(np,nq,alph,beta) */
/*   process. */

/*   input : */
/*           np,nq,alph,beta */
/*           mdim : dimension of various arrays in */
/*           calling program (see dimensions below) */
/*           mdim.ge.(2*max(np,nq))+1 */

/*   output : */
/*           ainf */
/*           ier : 0 is normal return, 1 means thereis */
/*           a singular matrix in the procedure */
/*           sd(1),...,sd(np+nq) */

/*   subroutines called : mxcvin,schur,swpk12,decomp,solv */

/* ********************************************************** */


/*   iab : */

    /* Parameter adjustments */
    --alph;
    --beta;
    ainf_dim1 = *mdim;
    ainf_offset = ainf_dim1 + 1;
    ainf -= ainf_offset;
    wk1_dim1 = *mdim;
    wk1_offset = wk1_dim1 + 1;
    wk1 -= wk1_offset;
    --r__;
    --sd;

    /* Function Body */
    if (*np * *nq == 0) {
	goto L99;
    }
    *ier = 0;
    on = (doublereal) (*n);
    m = max(*np,*nq);
    mm = (m << 1) + 1;
    nppnq = *np + *nq;
    mxcvin_(&alph[1], &beta[1], np, nq, &mm, mdim, &ainf[ainf_offset], &r__[1]
	    , ier);
    if (*ier == 1) {
	goto L99;
    }
    i__1 = nppnq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nppnq;
	for (j = 1; j <= i__2; ++j) {
/* L1: */
	    ainf[i__ + j * ainf_dim1] = (float)0.;
	}
    }
    i__2 = *np;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *nq;
	for (j = 1; j <= i__1; ++j) {
	    ainf[i__ + (*np + j) * ainf_dim1] = -r__[i__ - j + m + 1];
/* L2: */
	}
    }

/*   iaa : */

    schur_(&alph[1], mdim, np, &wk1[wk1_offset]);
    swpk12_(&wk1[wk1_offset], mdim, np, &c__1, np, ier);
    if (*ier == 1) {
	goto L99;
    }
    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *np;
	for (j = i__; j <= i__2; ++j) {
/* L3: */
	    ainf[i__ + j * ainf_dim1] = wk1[i__ + j * wk1_dim1];
	}
    }

/*   ibb : */

    schur_(&beta[1], mdim, nq, &wk1[wk1_offset]);
    swpk12_(&wk1[wk1_offset], mdim, nq, &c__1, nq, ier);
    if (*ier == 1) {
	goto L99;
    }
    i__2 = *nq;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = *np + i__;
	i__1 = *nq;
	for (j = i__; j <= i__1; ++j) {
/* L4: */
	    ainf[ii + (*np + j) * ainf_dim1] = wk1[i__ + j * wk1_dim1];
	}
    }

    i__1 = nppnq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nppnq;
	for (j = i__; j <= i__2; ++j) {
/* L5: */
	    ainf[j + i__ * ainf_dim1] = ainf[i__ + j * ainf_dim1];
	}
    }
    swpk12_(&ainf[ainf_offset], mdim, &nppnq, &c__1, &nppnq, ier);
    on = (doublereal) (*n);
    i__2 = nppnq;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L10: */
	sd[i__] = sqrt(ainf[i__ + i__ * ainf_dim1] / on);
    }
L99:

    return 0;
} /* mxcsd_ */


/* ********************************************************** */

/* Subroutine */ int mxcvin_(alph, beta, np, nq, mm, ndim, a, r__, ier)
doublereal *alph, *beta;
integer *np, *nq, *mm, *ndim;
doublereal *a, *r__;
integer *ier;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int solv_();
    static integer i__, j, k, m, m1, ip[100];
    extern /* Subroutine */ int decomp_();
    static integer jm1, mp1, mp2, mp3;


/*   subroutine to calculate the cross covariances for */
/*   lags (-max(np,nq),...,max(np,nq)), (stored in */
/*   r(1),...,r(mm),mm=(2*max(np,nq))+1) of the two */
/*   dimensional process (x1(t),x2(t)) where x1(.) and */
/*   x2(.) are defined by : */

/*      sum(j=0,np) alph(j)*x1(t-j)=e(t) */
/*      sum(j=0,nq)beta(j)*x2(t-j)=e(t) */

/*   input : */
/*           np,nq,mm=2*max(np,nq)+1 */
/*           alph,beta */
/*           ndim : dimension of a in calling program */
/*           (ndim.ge.mm) */

/*   output : */
/*           r(1),...,r(mm) (r(j) is r12(j-max(np,nq)-1) */
/*           ier : 0 is normal return, 1 means a is singular */

/*   subroutines called : decomp,solv */

/* ********************************************************** */

/*      data nout/0/ */



    /* Parameter adjustments */
    --alph;
    --beta;
    --r__;
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    m = max(*np,*nq);
    mp1 = m + 1;
    mp2 = m + 2;
    i__1 = *mm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *mm;
	for (j = 1; j <= i__2; ++j) {
/* L1: */
	    a[i__ + j * a_dim1] = (float)0.;
	}
    }


    if (*np < *nq) {
	goto L50;
    }

/*   for np.ge.nq : */

    a[a_dim1 + 1] = (float)1.;
    i__2 = *np;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L5: */
	a[(i__ + 1) * a_dim1 + 1] = alph[i__];
    }
    i__2 = mp1;
    for (j = 2; j <= i__2; ++j) {
	jm1 = j - 1;
	a[j + a_dim1] = a[jm1 + *mm * a_dim1];
	i__1 = *mm;
	for (k = 2; k <= i__1; ++k) {
/* L7: */
	    a[j + k * a_dim1] = a[jm1 + (k - 1) * a_dim1];
	}
/* L6: */
    }
    m1 = m - *nq + 1;
    i__2 = *nq;
    for (j = 1; j <= i__2; ++j) {
/* L8: */
	a[mp2 + (m1 + j) * a_dim1] = beta[*nq - j + 1];
    }
    a[mp2 + (m1 + *nq + 1) * a_dim1] = (float)1.;
    mp3 = m + 3;
    i__2 = *mm;
    for (j = mp3; j <= i__2; ++j) {
	jm1 = j - 1;
	a[j + a_dim1] = a[jm1 + *mm * a_dim1];
	i__1 = *mm;
	for (k = 2; k <= i__1; ++k) {
/* L10: */
	    a[j + k * a_dim1] = a[jm1 + (k - 1) * a_dim1];
	}
/* L9: */
    }
    goto L100;

/*   for np.lt.nq : */

L50:
    a[a_dim1 + 1] = (float)1.;
    i__2 = *np;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L51: */
	a[(i__ + 1) * a_dim1 + 1] = alph[i__];
    }
    if (m == 1) {
	goto L70;
    }
    i__2 = m;
    for (j = 2; j <= i__2; ++j) {
	jm1 = j - 1;
	a[j + a_dim1] = a[jm1 + *mm * a_dim1];
	i__1 = *mm;
	for (k = 2; k <= i__1; ++k) {
/* L53: */
	    a[j + k * a_dim1] = a[jm1 + (k - 1) * a_dim1];
	}
/* L52: */
    }
L70:
    i__2 = *nq;
    for (j = 1; j <= i__2; ++j) {
/* L71: */
	a[mp1 + j * a_dim1] = beta[*nq - j + 1];
    }
    a[mp1 + (*nq + 1) * a_dim1] = (float)1.;
    i__2 = *mm;
    for (j = mp2; j <= i__2; ++j) {
	jm1 = j - 1;
	a[j + a_dim1] = a[jm1 + *mm * a_dim1];
	i__1 = *mm;
	for (k = 2; k <= i__1; ++k) {
/* L73: */
	    a[j + k * a_dim1] = a[jm1 + (k - 1) * a_dim1];
	}
/* L72: */
    }

/*   solve system of equations : */

L100:
    i__2 = *mm;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L101: */
	r__[i__] = (float)0.;
    }
    r__[mp1] = (float)1.;


    decomp_(mm, ndim, &a[a_offset], ip);
    if (ip[*mm - 1] == 0) {
	goto L125;
    }
    *ier = 0;
    solv_(mm, ndim, &a[a_offset], &r__[1], ip);

/*   r contains r12(-m),...,r12(m) */

    return 0;

L125:
    *ier = 1;
    return 0;
} /* mxcvin_ */

/**************************************************************************/

/* Subroutine */ int corrar_(rho, r0, p, alpha, rvar)
doublereal *rho, *r0;
integer *p;
doublereal *alpha, *rvar;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal last, temp;
    static integer j, k;


/************************************************************************
**/



    /* Parameter adjustments */
    --alpha;
    --rho;

    /* Function Body */
    alpha[1] = -rho[1];
    *rvar = 1 - rho[1] * rho[1];
    if (*p == 1) {
	goto L99;
    }

    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {

	last = -rho[j];
	i__2 = j - 1;
	for (k = 1; k <= i__2; ++k) {
/* L10: */
	    last -= alpha[k] * rho[j - k];
	}
	last /= *rvar;

	alpha[j] = last;
	i__2 = j / 2;
	for (k = 1; k <= i__2; ++k) {
	    temp = alpha[k];
	    alpha[k] += last * alpha[j - k];
	    if (k == j - k) {
		goto L20;
	    }
	    alpha[j - k] += last * temp;
L20:
	    ;
	}

	*rvar *= 1 - last * last;

/* L30: */
    }

L99:
    *rvar = *r0 * *rvar;


    return 0;
} /* corrar_ */

/* ******************************************************************* */

/* Subroutine */ int cvmx1_(r__, r0, np, nq, ndim, maxit, del, ip, al, alpha, 
	beta, rvar, ier)
doublereal *r__, *r0;
integer *np, *nq, *ndim, *maxit;
doublereal *del;
integer *ip;
doublereal *al, *alpha, *beta, *rvar;
integer *ier;
{
    /* System generated locals */
    integer al_dim1, al_offset, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int macv_(), hiyw_();
    static doublereal c__, d__[100], f[100], g[100];
    static integer i__, j;
    static doublereal t[100], c1;
    static integer n1;
    static doublereal dn[100], ra[100], fn[100], fp[100], ry[100];
    extern /* Subroutine */ int cvmawl_();
    static doublereal ra0, ry0, gam[100];


/*        subroutine cvmx1(r,r0,np,nq,ndim,maxit,del,ip,al,ra,ry,g, */
/*     1  alpha,beta,rvar,ier,t,f,fn,fp,d,dn,gam) */
/*   Subroutine to find arma paramters from covariances. */

/*   ndim.ge.np+nq */

/* ****************************************************************** */

/*        dimension r(ndim),al(ndim,ndim),ra(ndim),ry(ndim),alpha(ndim), 
*/
/*     1  beta(ndim),ip(ndim),g(ndim),t(ndim),f(ndim),fn(ndim),d(ndim), */
/*     1	dn(ndim),gam(ndim),fp(ndim) */

/*   find alpha : */

    /* Parameter adjustments */
    --beta;
    --alpha;
    al_dim1 = *ndim;
    al_offset = al_dim1 + 1;
    al -= al_offset;
    --ip;
    --r__;

    /* Function Body */
    hiyw_(&r__[1], r0, np, nq, ndim, &ip[1], &al[al_offset], &alpha[1], ier);
    if (*ier == 1) {
	*ier = 3;
	goto L99;
    }

/*   find ra : */

    macv_(&alpha[1], np, &c_b54, ra, &ra0);

/*   find ry : */

    ry0 = ra0 * *r0;
    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L12: */
	ry0 += ra[i__ - 1] * (float)2. * r__[i__];
    }
    i__1 = *nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__ = ra0 * r__[i__];
	i__2 = *np;
	for (j = 1; j <= i__2; ++j) {
	    if (i__ - j == 0) {
		goto L16;
	    }
	    n1 = (i__3 = i__ - j, abs(i__3));
	    c1 = r__[n1];
	    goto L14;
L16:
	    c1 = *r0;
L14:
	    c__ += ra[j - 1] * (c1 + r__[i__ + j]);
	}
/* L13: */
	ry[i__ - 1] = c__;
    }

/*   find corresponding betas : */

    cvmawl_(ry, &ry0, nq, maxit, del, t, f, fn, fp, g, d__, dn, gam, &beta[1],
	     rvar, ier);
L99:
    return 0;
} /* cvmx1_ */


/* ****************************************************************** */

/* Subroutine */ int hiyw_(r__, r0, np, nq, ndim, ip, al, alpha, ier)
doublereal *r__, *r0;
integer *np, *nq, *ndim, *ip;
doublereal *al, *alpha;
integer *ier;
{
    /* System generated locals */
    integer al_dim1, al_offset, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int solv_();
    static integer i__, j, n1;
    extern /* Subroutine */ int decomp_();


/*   Subroutine to solve high order Yule-Walker equations. ndim.ge.np+nq 
*/

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --alpha;
    al_dim1 = *ndim;
    al_offset = al_dim1 + 1;
    al -= al_offset;
    --ip;
    --r__;

    /* Function Body */
    *ier = 0;
    if (*np == 1) {
	goto L10;
    }

    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	alpha[i__] = -r__[*nq + i__];
	i__2 = *np;
	for (j = i__; j <= i__2; ++j) {
	    if (*nq + i__ - j == 0) {
		goto L3;
	    }
	    n1 = (i__3 = *nq + i__ - j, abs(i__3));
	    al[i__ + j * al_dim1] = r__[n1];
	    goto L1;
L3:
	    al[i__ + j * al_dim1] = *r0;
L1:
	    al[j + i__ * al_dim1] = r__[*nq + j - i__];
	}
    }

    decomp_(np, ndim, &al[al_offset], &ip[1]);
    if (ip[*np] == 0) {
	goto L99;
    }
    solv_(np, ndim, &al[al_offset], &alpha[1], &ip[1]);
    return 0;
L10:
    alpha[1] = -r__[*nq + 1] / r__[*nq];
    return 0;
L99:
    *ier = 1;
    return 0;
} /* hiyw_ */


/* ****************************************************************** */

/* Subroutine */ int cvmawl_(r__, r0, nq, maxit, del, t, f, fn, fp, g, d__, 
	dn, gam, beta, rvar, ier)
doublereal *r__, *r0;
integer *nq, *maxit;
doublereal *del, *t, *f, *fn, *fp, *g, *d__, *dn, *gam, *beta, *rvar;
integer *ier;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublereal d0, e1, g0, t0;
    static integer it;
    static doublereal dn0, fpk, eps, fpk2;


/*   Subroutine to do Wilson's algorithm for MA processes. */

/* ****************************************************************** */


/*   starting values: */

    /* Parameter adjustments */
    --beta;
    --gam;
    --dn;
    --d__;
    --g;
    --fp;
    --fn;
    --f;
    --t;
    --r__;

    /* Function Body */
    t0 = sqrt(*r0);
    i__1 = *nq;
    for (j = 1; j <= i__1; ++j) {
/* L10: */
	t[j] = r__[j] / t0;
    }

/*   start iterations (t plays the role of tau, g plays role of r): */

    i__1 = *maxit;
    for (it = 1; it <= i__1; ++it) {


	i__2 = *nq;
	for (j = 1; j <= i__2; ++j) {
/* L20: */
	    f[j] = t[j] / t0;
	}
	g0 = t0 * t0;
	i__2 = *nq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L22: */
	    g0 += t[i__] * t[i__];
	}
	i__2 = *nq;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__ = t0 * t[i__];
	    if (i__ == *nq) {
		goto L26;
	    }
	    i__3 = *nq - i__;
	    for (j = 1; j <= i__3; ++j) {
/* L25: */
		c__ += t[j] * t[j + i__];
	    }
L26:
	    g[i__] = c__;
/* L24: */
	}
	g0 = (g0 + *r0) / t0;
	i__2 = *nq;
	for (j = 1; j <= i__2; ++j) {
/* L30: */
	    g[j] = (g[j] + r__[j]) / t0;
	}

/*   work backward: */

	for (k = *nq; k >= 1; --k) {
	    gam[k] = g[k];
	    fp[k] = f[k];
	    fpk = fp[k];
	    fpk2 = 1. - fpk * fpk;
	    if (fpk2 <= (float)0.) {
		goto L299;
	    }
	    if (k == 1) {
		goto L70;
	    }
	    i__2 = k - 1;
	    for (j = 1; j <= i__2; ++j) {
/* L40: */
		fn[j] = (f[j] - fpk * f[k - j]) / fpk2;
	    }
	    i__2 = k - 1;
	    for (j = 1; j <= i__2; ++j) {
/* L50: */
		g[j] -= gam[k] * fn[k - j];
	    }
	    i__2 = k - 1;
	    for (j = 1; j <= i__2; ++j) {
/* L60: */
		f[j] = fn[j];
	    }
L70:
	    ;
	}

/*   work forward: */

	d0 = g0 / (float)2.;
	i__2 = *nq - 1;
	for (k = 0; k <= i__2; ++k) {
	    d__[k + 1] = gam[k + 1];
	    fpk = fp[k + 1];
	    fpk2 = 1. - fpk * fpk;
	    dn0 = (d0 - fpk * d__[k + 1]) / fpk2;
	    dn[k + 1] = (d__[k + 1] - fpk * d0) / fpk2;
	    if (k == 0) {
		goto L90;
	    }
	    i__3 = k;
	    for (j = 1; j <= i__3; ++j) {
/* L80: */
		dn[j] = (d__[j] - fpk * d__[k + 1 - j]) / fpk2;
	    }
L90:
	    d0 = dn0;
	    i__3 = k + 1;
	    for (j = 1; j <= i__3; ++j) {
/* L100: */
		d__[j] = dn[j];
	    }
/* L110: */
	}

/*   check convergence: */

	eps = (d__1 = t0 - d0, abs(d__1)) / abs(d0);
	i__2 = *nq;
	for (j = 1; j <= i__2; ++j) {
	    e1 = (d__1 = t[j] - d__[j], abs(d__1)) / (d__2 = d__[j], abs(d__2)
		    );
	    if (e1 > eps) {
		eps = e1;
	    }
/* L120: */
	}
	t0 = d0;
	i__2 = *nq;
	for (j = 1; j <= i__2; ++j) {
/* L130: */
	    t[j] = d__[j];
	}

/*   yes: */

	if (eps < *del) {
	    *ier = 0;
	    *rvar = t0 * t0;
	    i__2 = *nq;
	    for (j = 1; j <= i__2; ++j) {
/* L140: */
		beta[j] = t[j] / t0;
	    }
	    goto L99;
	}

/*   no: */

/* L150: */
    }

/*   nonconvergence in maxit iterations: */

    *ier = 1;
    goto L99;

/*   partial outside (-1,1): */

L299:
    *ier = 2;


L99:
    return 0;
} /* cvmawl_ */


/* *********************************************************** */

/* Subroutine */ int decomp_(n, ndim, a, ip)
integer *n, *ndim;
doublereal *a;
integer *ip;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal t;
    static integer kp1;


/*   matrix triangularization by gaussian elimination. */

/*   input.... */
/*     n = order of matrix */
/*     ndim = declared dimension of array a */
/*     a = matrix to be triangularized */

/*   output.... */
/*     a(i,j), i.le.j = upper triangular factor, u. */
/*     a(i,j), i.gt.j = multipliers = lower triangular */
/*                      factor, i-l */
/*     ip(k),  k.lt.n = index of k-th pivot row */
/*     ip(n) = (-1)**(number of interchanges) or 0 */

/*   use subroutine solve to obtain solution of linear system. */
/*   determ(a) = ip(n)*a(1,1)*a(2,2)*...*a(n,n) */
/*   if ip(n)=0, a is singular,  solve  will divide by zero */
/*   interchanges finished in u, only partly in l. */

/*   reference...algorithm 423 'linear equation solver' */
/*               by cleve b. moler, cacm, april 1972 p. 274 */

/* *********************************************************** */


    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    ip[*n] = 1;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (k == *n) {
	    goto L5;
	}
	kp1 = k + 1;
	m = k;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/* L1: */
	}
	ip[k] = m;
	if (m != k) {
	    ip[*n] = -ip[*n];
	}
	t = a[m + k * a_dim1];
	a[m + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
	if (t == (float)0.) {
	    goto L5;
	}
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L2: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] / t;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[m + j * a_dim1];
	    a[m + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
	    if (t == (float)0.) {
		goto L4;
	    }
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/* L3: */
		a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
	    }
L4:
	    ;
	}
L5:
	if (a[k + k * a_dim1] == (float)0.) {
	    ip[*n] = 0;
	}
/* L6: */
    }
    return 0;
} /* decomp_ */


/* ****************************************************************** */

/* Subroutine */ int solv_(n, ndim, a, b, ip)
integer *n, *ndim;
doublereal *a, *b;
integer *ip;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, m;
    static doublereal t;
    static integer kb, km1, nm1, kp1;


/*     solution of linear system  a*x=b */

/*     do not use if  decomp  has set ip(n)=0 */

/*     input.... */
/*     n = order of matrix a */
/*     ndim = declared dimension of array a */
/*     a = triangularized matrix obtained from subroutine |decomp| */
/*     b = right hand side vector */
/*     ip = pivot vector obtained from |decomp| */

/*     output.... */
/*     b = solution vector, x. */

/*     reference...algorithm 423 'linear equation solver' */
/*     by cleve b. moler, cacm, april 1972 p. 274 */

/* ****************************************************************** */


    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*n == 1) {
	goto L9;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = ip[k];
	t = b[m];
	b[m] = b[k];
	b[k] = t;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L7: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
    }
    i__2 = nm1;
    for (kb = 1; kb <= i__2; ++kb) {
	km1 = *n - kb;
	k = km1 + 1;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__1 = km1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L8: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
    }
L9:
    b[1] /= a[a_dim1 + 1];
    return 0;
} /* solv_ */



/* ******************************************************* */

/* Subroutine */ int macv_(beta, nq, sig, r__, r0)
doublereal *beta;
integer *nq;
doublereal *sig, *r__, *r0;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer nqmi;
    static doublereal c__;
    static integer i__, j;


/*     subroutine to calculate the autocovariances r0,r(1), */
/*     ...,r(nq) for a moving average process of order nq */
/*     with parameters beta(1),...,beta(nq), and sig (res var) */

/*     input : */
/*     nq,beta(1),...,beta(nq),sig */

/*     output : */
/*     r0,r(1),...,r(nq) */

/*     subroutines called : none */

/* ******************************************************* */


    /* Parameter adjustments */
    --r__;
    --beta;

    /* Function Body */
    c__ = 1.;
    i__1 = *nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	c__ += beta[i__] * beta[i__];
    }
    *r0 = c__ * *sig;

    i__1 = *nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__ = beta[i__];
	if (i__ == *nq) {
	    goto L2;
	}
	nqmi = *nq - i__;
	i__2 = nqmi;
	for (j = 1; j <= i__2; ++j) {
/* L3: */
	    c__ += beta[j] * beta[j + i__];
	}
L2:
	r__[i__] = c__ * *sig;
    }

    return 0;
} /* macv_ */


/* ********************************************************************** */

/* Subroutine */ int diffeq_(alpha, np, n, e, x)
doublereal *alpha;
integer *np, *n;
doublereal *e, *x;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;


/* ********************************************************************** 
*/



    /* Parameter adjustments */
    --alpha;
    --x;
    --e;

    /* Function Body */
    if (*np == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L5: */
	    x[i__] = e[i__];
	}
	return 0;
    }

    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = e[i__];
    }
    i__1 = *n;
    for (i__ = *np + 1; i__ <= i__1; ++i__) {
	c__ = e[i__];
	i__2 = *np;
	for (j = 1; j <= i__2; ++j) {
/* L30: */
	    c__ -= alpha[j] * x[i__ - j];
	}
/* L20: */
	x[i__] = c__;
    }


    return 0;
} /* diffeq_ */


/* ********************************************************************* */

/* Subroutine */ int dtarma_(start, eps, maxit, dat, ny, e, np, nq, alpha, 
	beta, nppnq, p, rv, am2ll, ier, var)
doublereal *start, *eps;
integer *maxit;
doublereal *dat;
integer *ny;
doublereal *e;
integer *np, *nq;
doublereal *alpha, *beta;
integer *nppnq;
doublereal *p, *rv, *am2ll;
integer *ier;
doublereal *var;
{
    /* Initialized data */

    static doublereal rcoeff = 1.;
    static doublereal ecoeff = 2.;
    static doublereal ccoeff = .5;

    /* System generated locals */
    integer p_dim1, p_offset, i__1, i__2;

    /* Local variables */
    static doublereal pbar[20], cmin, step[20], summ;
    static integer i__, j, l, n;
    static doublereal y[20], z__, pstar[20], ystar, p2star[20], dn, y2star;
    static integer nn;
    static doublereal vk[20], vl[20];
    extern doublereal armalk_();
    static doublereal vw[20];
    static integer konvge;
    static doublereal curmin;
    static integer jcount, kcount;
    static doublereal ynewlo, del;
    static integer ihi;
    static doublereal dnn, min__[20];
    static integer ilo;
    static doublereal ylo, sum;



/*   Subroutine to find exact ARMA MLE's as described in TIMESLAB. */

/*   Input: */
/*      start: starting values (those for alpha, then for beta) */
/*      eps  : convergence criterion */
/*      maxit: maximum number of function evaluations */
/*      dat  : data */
/*      ny   : number of data points */
/*      np,nq: order of ARMA */
/*      nppnq: np+nq */

/*   Output: */
/*      alpha, beta, rv: values of parameters when algorithm terminates */
/*      am2ll : -2log(likelihhod) at last values */
/*      ier   : <=0 means convergence, 1 means nonconvergence, >1 means */
/*              error. */
/*      var   : variance of vertices upon termination */

/*   Auxilliary: */
/*      p(np+nq,np+nq+1) */
/*      e(ny) */

/*   Used routines: armalk,partar,lkhood,wilson */

/* ******************************************************************** */


    /* Parameter adjustments */
    --start;
    --e;
    --dat;
    --alpha;
    --beta;
    p_dim1 = *nppnq;
    p_offset = p_dim1 + 1;
    p -= p_offset;

    /* Function Body */


    n = *np + *nq;
    konvge = 10;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L310: */
	step[i__ - 1] = (float)1.;
    }
    *ier = 0;
    kcount = *maxit;
    *maxit = 0;
    if (*eps <= (float)0.) {
	--(*maxit);
    }
    if (n > 20) {
	*maxit += -10;
    }
    if (konvge <= 0) {
	*maxit += -100;
    }
    if (*maxit < 0) {
	return 0;
    }

    jcount = konvge;
    dn = (doublereal) n;
    nn = n + 1;
    dnn = (doublereal) nn;
    del = (float)1.;

/*   construction of initial simplex: */

L1001:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	p[i__ + nn * p_dim1] = start[i__];
    }
    z__ = armalk_(&start[1], np, nq, ny, &dat[1], &alpha[1], &beta[1], &e[1], 
	    vw, vl, vk, rv, ier);
    if (*ier > 0) {
	return 0;
    }
    y[nn - 1] = z__;
    sum = z__;
    summ = z__ * z__;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	start[j] += step[j - 1] * del;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L3: */
	    p[i__ + j * p_dim1] = start[i__];
	}
	z__ = armalk_(&start[1], np, nq, ny, &dat[1], &alpha[1], &beta[1], &e[
		1], vw, vl, vk, rv, ier);
	if (*ier > 0) {
	    return 0;
	}
	y[j - 1] = z__;
	sum += z__;
	summ += z__ * z__;
/* L2: */
	start[j] -= step[j - 1] * del;
    }

/*   simplex construction complete */

/*   find highest and lowest y values. ynewlo (=y(ihi)) indicates */
/*   the vertex of the simplex to be replaced. */

L1000:
    ylo = y[0];
    ynewlo = ylo;
    ilo = 1;
    ihi = 1;
    i__1 = nn;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (y[i__ - 1] >= ylo) {
	    goto L4;
	}
	ylo = y[i__ - 1];
	ilo = i__;
L4:
	if (y[i__ - 1] <= ynewlo) {
	    goto L5;
	}
	ynewlo = y[i__ - 1];
	ihi = i__;
L5:
	;
    }
    sum -= ynewlo;
    summ -= ynewlo * ynewlo;

/*   calculate pbar. The centroid of the simplex vertices */
/*   excepting that with y value ynewlo. */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = (float)0.;
	i__2 = nn;
	for (j = 1; j <= i__2; ++j) {
/* L6: */
	    z__ += p[i__ + j * p_dim1];
	}
	z__ -= p[i__ + ihi * p_dim1];
/* L7: */
	pbar[i__ - 1] = z__ / dn;
    }

/*   reflection through the centroid */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L8: */
	pstar[i__ - 1] = (rcoeff + (float)1.) * pbar[i__ - 1] - rcoeff * p[
		i__ + ihi * p_dim1];
    }
    ystar = armalk_(pstar, np, nq, ny, &dat[1], &alpha[1], &beta[1], &e[1], 
	    vw, vl, vk, rv, ier);
    if (*ier > 0) {
	return 0;
    }
    ++(*maxit);
    if (ystar >= ylo) {
	goto L12;
    }

/*   successful reflection, so extension */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L9: */
	p2star[i__ - 1] = ecoeff * pstar[i__ - 1] + ((float)1. - ecoeff) * 
		pbar[i__ - 1];
    }
    y2star = armalk_(p2star, np, nq, ny, &dat[1], &alpha[1], &beta[1], &e[1], 
	    vw, vl, vk, rv, ier);
    if (*ier > 0) {
	return 0;
    }
    ++(*maxit);

/*   retain extension or contraction */

    if (y2star >= ylo) {
	goto L19;
    }
L10:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L11: */
	p[i__ + ihi * p_dim1] = p2star[i__ - 1];
    }
    y[ihi - 1] = y2star;
    goto L900;
/*  no extension */
L12:
    l = 0;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (y[i__ - 1] > ystar) {
	    ++l;
	}
/* L13: */
    }
    if (l > 1) {
	goto L19;
    }
    if (l == 0) {
	goto L15;
    }

/*  contraction on the reflection side of the centroid */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L14: */
	p[i__ + ihi * p_dim1] = pstar[i__ - 1];
    }
    y[ihi - 1] = ystar;

/*   contraction on the y(ihi) side of the centroid */

L15:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L16: */
	p2star[i__ - 1] = ccoeff * p[i__ + ihi * p_dim1] + ((float)1. - 
		ccoeff) * pbar[i__ - 1];
    }
    y2star = armalk_(p2star, np, nq, ny, &dat[1], &alpha[1], &beta[1], &e[1], 
	    vw, vl, vk, rv, ier);
    if (*ier > 0) {
	return 0;
    }
    ++(*maxit);
    if (y2star <= y[ihi - 1]) {
	goto L10;
    }

/*   contract whole simplex */

    sum = 0.;
    summ = 0.;
    i__1 = nn;
    for (j = 1; j <= i__1; ++j) {
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p[i__ + j * p_dim1] = (p[i__ + j * p_dim1] + p[i__ + ilo * p_dim1]
		    ) * (float).5;
/* L17: */
	    min__[i__ - 1] = p[i__ + j * p_dim1];
	}
	y[j - 1] = armalk_(min__, np, nq, ny, &dat[1], &alpha[1], &beta[1], &
		e[1], vw, vl, vk, rv, ier);
	if (*ier > 0) {
	    return 0;
	}
	sum += y[j - 1];
/* L18: */
	summ += y[j - 1] * y[j - 1];
    }
    *maxit += nn;
    goto L901;
/*   retain reflection */
L19:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	p[i__ + ihi * p_dim1] = pstar[i__ - 1];
    }
    y[ihi - 1] = ystar;
L900:
    sum += y[ihi - 1];
    summ += y[ihi - 1] * y[ihi - 1];
L901:
    --jcount;
    if (jcount != 0) {
	goto L1000;
    }

/*   check to see if minimum reached */

    if (*maxit > kcount) {
	goto L22;
    }
    jcount = konvge;
    cmin = (summ - sum * sum / dnn) / dn;
    curmin = cmin;

/*   curmin is the variance of the n+1 fn values at the vertices */

    if (curmin >= *eps) {
	goto L1000;
    }

/*   factorial tests to check that ynewlo is a local minimum */

L22:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L23: */
	min__[i__ - 1] = p[i__ + ihi * p_dim1];
    }
    ynewlo = y[ihi - 1];
    if (*maxit > kcount) {
	*ier = 1;
	return 0;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	del = step[i__ - 1] * (float).001;
	min__[i__ - 1] += del;
	z__ = armalk_(min__, np, nq, ny, &dat[1], &alpha[1], &beta[1], &e[1], 
		vw, vl, vk, rv, ier);
	if (*ier > 0) {
	    return 0;
	}
	if (z__ < ynewlo) {
	    goto L25;
	}
	min__[i__ - 1] = min__[i__ - 1] - del - del;
	z__ = armalk_(min__, np, nq, ny, &dat[1], &alpha[1], &beta[1], &e[1], 
		vw, vl, vk, rv, ier);
	if (*ier > 0) {
	    return 0;
	}
	if (z__ < ynewlo) {
	    goto L25;
	}
/* L24: */
	min__[i__ - 1] += del;
    }
    *am2ll = z__;
    *var = curmin;
    return 0;

/*   restart procedure */

L25:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L26: */
	start[i__] = min__[i__ - 1];
    }
    del = (float).001;
    ++(*maxit);
    goto L1001;
} /* dtarma_ */


/* ******************************************************************** */

doublereal armalk_(theta, np, nq, n, y, alpha, beta, e, vw, vl, vk, rv, ier)
doublereal *theta;
integer *np, *nq, *n;
doublereal *y, *alpha, *beta, *e, *vw, *vl, *vk, *rv;
integer *ier;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double exp(), atan(), log();

    /* Local variables */
    static doublereal fact, al2pi, c__;
    static integer i__;
    static doublereal e1, sumsq;
    extern /* Subroutine */ int lkhood_(), partar_();
    static doublereal tol;


/*   Subroutine to calculate -2log likelihood for ARMA(np,nq) process */

/* ********************************************************************* 
*/


/*   find -alpha and -beta corresponding to input transformed partials: */

    /* Parameter adjustments */
    --theta;
    --e;
    --y;
    --alpha;
    --beta;
    --vw;
    --vl;
    --vk;

    /* Function Body */
    if (*np > 0) {
	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    e1 = exp(-theta[i__]);
/* L5: */
	    alpha[i__] = ((float)1. - e1) / (e1 + (float)1.);
	}
	partar_(&alpha[1], np);
	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	    alpha[i__] = -alpha[i__];
	}
    }
    if (*nq > 0) {
	i__1 = *nq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    e1 = exp(-theta[*np + i__]);
/* L20: */
	    beta[i__] = ((float)1. - e1) / (e1 + (float)1.);
	}
	partar_(&beta[1], nq);
	i__1 = *nq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	    beta[i__] = -beta[i__];
	}
    }

/*   call evaluator: */

    tol = (float)1e-4;
    lkhood_(&alpha[1], np, &beta[1], nq, &y[1], &e[1], n, &sumsq, &fact, &vw[
	    1], &vl[1], &vk[1], &tol, ier);


    if (*ier > 0) {
	return ret_val;
    }
    if (*np > 0) {
	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L35: */
	    alpha[i__] = -alpha[i__];
	}
    }
    if (*nq > 0) {
	i__1 = *nq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	    beta[i__] = -beta[i__];
	}
    }
    *rv = sumsq / (doublereal) (*n);
/* 	al2pi=alog(8.*atan(1.0)) */
/* 	c=float(n)*(1.+al2pi+alog(float(rv))+alog(float(fact))) */


    al2pi = log(atan(1.) * 8.);
    c__ = (doublereal) (*n) * (al2pi + (float)1. + log(*rv) + log(fact));


    ret_val = c__;
    return ret_val;
} /* armalk_ */


/* ******************************************************************** */

/* Subroutine */ int lkhood_(alpha, np, beta, nq, y, e, n, sumsq, fact, vw, 
	vl, vk, tol, ier)
doublereal *alpha;
integer *np;
doublereal *beta;
integer *nq;
doublereal *y, *e;
integer *n;
doublereal *sumsq, *fact, *vw, *vl, *vk, *tol;
integer *ier;
{
    /* Initialized data */

    static real epsil1 = (float)1e-10;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), pow_dd();

    /* Local variables */
    static integer last, loop, mxpq;
    static doublereal a;
    static integer i__, j, k;
    static doublereal r__;
    static integer jfrom, nexti;
    static doublereal fn;
    static integer mxpqp1, mr;
    static doublereal detcar, detman;
    extern /* Subroutine */ int wilson_();
    static doublereal vl1, vw1, alf, flj, aor;
    static integer mpp1, mqp1, mrp1;


/*   Subroutine to find the 2 terms in the exact ARMA likelihood for */
/*   -alpha and -beta as in Melard (Applied Statistics, 1984, 104) */

/*   vw and vl are max(p,q+1)+1 vk is max(p,q+1) */

/* ******************************************************************** */

    /* Parameter adjustments */
    --alpha;
    --beta;
    --e;
    --y;
    --vw;
    --vl;
    --vk;

    /* Function Body */


/* Computing MAX */
    i__1 = *np, i__2 = *nq + 1;
    mr = max(i__1,i__2);
    mrp1 = mr + 1;
    *fact = 0.;
    detman = 1.;
    detcar = 0.;
    *sumsq = 0.;
    mxpq = max(*np,*nq);
    mxpqp1 = mxpq + 1;
    mqp1 = *nq + 1;
    mpp1 = *np + 1;


    wilson_(&alpha[1], np, &beta[1], nq, &vw[1], &mxpqp1, &vl[1], &mxpqp1, &
	    vk[1], &mxpq, ier);
    if (*ier > 0) {
	return 0;
    }
    vk[1] = vw[1];
    if (mr == 1) {
	goto L15;
    }
    i__1 = mr;
    for (k = 2; k <= i__1; ++k) {
	vk[k] = (float)0.;
	if (k > *np) {
	    goto L12;
	}
	i__2 = *np;
	for (j = k; j <= i__2; ++j) {
/* L11: */
	    vk[k] += alpha[j] * vw[j + 2 - k];
	}
L12:
	if (k > mqp1) {
	    goto L14;
	}
	i__2 = mqp1;
	for (j = k; j <= i__2; ++j) {
/* L13: */
	    vk[k] -= beta[j - 1] * vl[j + 1 - k];
	}
L14:
	;
    }


L15:
    r__ = vk[1];
    vl[mr] = (float)0.;
    i__1 = mr;
    for (j = 1; j <= i__1; ++j) {
	vw[j] = (float)0.;
	if (j != mr) {
	    vl[j] = vk[j + 1];
	}
	if (j <= *np) {
	    vl[j] += alpha[j] * r__;
	}
/* L16: */
	vk[j] = vl[j];
    }


    last = mpp1 - *nq;
    loop = *np;
    jfrom = mpp1;
    vw[mpp1] = (float)0.;
    vl[mxpqp1] = (float)0.;


    if (*n <= 0) {
	goto L50;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {


	if (i__ != last) {
	    goto L17;
	}
	loop = min(*np,*nq);
	jfrom = loop + 1;


	if (*nq <= 0) {
	    goto L30;
	}
L17:
	if (r__ <= epsil1) {
	    goto L40;
	}
	if ((d__1 = r__ - (float)1., abs(d__1)) < *tol && i__ > mxpq) {
	    goto L30;
	}


	detman *= r__;
L19:
	if (abs(detman) < (float)1.) {
	    goto L20;
	}
	detman *= (float).0625;
	detcar += (float)4.;
	goto L19;
L20:
	if (abs(detman) >= (float).0625) {
	    goto L21;
	}
	detman *= (float)16.;
	detcar += (float)-4.;
	goto L20;
L21:
	vw1 = vw[1];
	a = y[i__] - vw1;
	e[i__] = a / sqrt(r__);
	aor = a / r__;
	*sumsq += a * aor;
	vl1 = vl[1];
	alf = vl1 / r__;
	r__ -= alf * vl1;
	if (loop == 0) {
	    goto L23;
	}


	i__2 = loop;
	for (j = 1; j <= i__2; ++j) {
	    flj = vl[j + 1] + alpha[j] * vl1;
	    vw[j] = vw[j + 1] + alpha[j] * vw1 + aor * vk[j];
	    vl[j] = flj - alf * vk[j];
/* L22: */
	    vk[j] -= alf * flj;
	}
L23:
	if (jfrom > *nq) {
	    goto L25;
	}
	i__2 = *nq;
	for (j = jfrom; j <= i__2; ++j) {
	    vw[j] = vw[j + 1] + aor * vk[j];
	    vl[j] = vl[j + 1] - alf * vk[j];
/* L24: */
	    vk[j] -= alf * vl[j + 1];
	}
L25:
	if (jfrom > *np) {
	    goto L27;
	}
	i__2 = *np;
	for (j = jfrom; j <= i__2; ++j) {
/* L26: */
	    vw[j] = vw[j + 1] + alpha[j] * y[i__];
	}
L27:
/* L29: */
	;
    }
    goto L39;


L30:
    nexti = i__;
    *ier = -nexti;
    i__1 = *n;
    for (i__ = nexti; i__ <= i__1; ++i__) {
/* L31: */
	e[i__] = y[i__];
    }
    if (*np == 0) {
	goto L34;
    }
    i__1 = *n;
    for (i__ = nexti; i__ <= i__1; ++i__) {
	i__2 = *np;
	for (j = 1; j <= i__2; ++j) {
/* L32: */
	    e[i__] -= alpha[j] * y[i__ - j];
	}
/* L33: */
    }
L34:
    if (*nq == 0) {
	goto L37;
    }
    i__1 = *n;
    for (i__ = nexti; i__ <= i__1; ++i__) {
	i__2 = *nq;
	for (j = 1; j <= i__2; ++j) {
/* L35: */
	    e[i__] += beta[j] * e[i__ - j];
	}
/* L36: */
    }


L37:
    i__1 = *n;
    for (i__ = nexti; i__ <= i__1; ++i__) {
/* L38: */
	*sumsq += e[i__] * e[i__];
    }
L39:
    fn = (doublereal) (*n);
    d__1 = (float)1. / fn;
    d__2 = detcar / fn;
    *fact = pow_dd(&detman, &d__1) * pow_dd(&c_b225, &d__2);
    return 0;
L40:
    *ier = 8;
    return 0;
L50:
    *ier = 9;
    return 0;
} /* lkhood_ */

/***************************************************************************/

/* Subroutine */ int filt_(beta, beta0, m, n, x, y)
doublereal *beta, *beta0;
integer *m, *n;
doublereal *x, *y;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer j, t, tpm;


/*   y(t) = beta0*x(t) + sum(j=1,m) beta(j)x(t-j),  t=m+1,...,n */

/************************************************************************
***/


    /* Parameter adjustments */
    --y;
    --x;
    --beta;

    /* Function Body */
    i__1 = *n - *m;
    for (t = 1; t <= i__1; ++t) {
	tpm = t + *m;
	c__ = *beta0 * x[tpm];
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
/* L10: */
	    c__ += beta[j] * x[tpm - j];
	}
/* L20: */
	y[t] = c__;
    }

    return 0;
} /* filt_ */

/* ********************************************************************* */

/* Subroutine */ int movave_(x, n, k2, xave)
doublereal *x;
integer *n, *k2;
doublereal *xave;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, ii, jj, oi, ok, kp1;


/*   Find moving averages of length k2 (except at edges where */
/*   as many as available are used). */

/* ********************************************************************** 
*/


/*   Get sums: */

    /* Parameter adjustments */
    --xave;
    --x;

    /* Function Body */
    k = (*k2 - 1) / 2;
    xave[1] = x[1];
    xave[*n] = x[*n];
    ii = 1;
    jj = *n;
    i__1 = k + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ii += 2;
	jj += -2;
	xave[*n - i__ + 1] = xave[*n - i__ + 2] + x[jj] + x[jj + 1];
/* L10: */
	xave[i__] = xave[i__ - 1] + x[ii] + x[ii - 1];
    }

    kp1 = k + 1;
    i__1 = *n - k - 1;
    for (i__ = k + 2; i__ <= i__1; ++i__) {
/* L20: */
	xave[i__] = xave[i__ - 1] + x[i__ + k] - x[i__ - kp1];
    }

    i__1 = k + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	oi = (integer) ((doublereal) ((i__ << 1) - 1));
	xave[i__] /= oi;
/* L30: */
	xave[*n - i__ + 1] /= oi;
    }

    ok = (k << 1) + 1;
    i__1 = *n - k - 1;
    for (i__ = k + 2; i__ <= i__1; ++i__) {
/* L40: */
	xave[i__] /= ok;
    }

    return 0;
} /* movave_ */


/* *********************************************************************** */

/* Subroutine */ int movbox_(x, n, k2, ii, jj, xsumm, nouts, indout, outval)
doublereal *x;
integer *n, *k2;
doublereal *ii, *jj, *xsumm;
integer *nouts, *indout;
doublereal *outval;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int summ5_();
    static integer k, m;
    extern /* Subroutine */ int medadd_();
    static integer nl;
    extern /* Subroutine */ int meddel_();


/*   Find moving five-pt summaries and extreme values (as in a box-plot) 
*/
/*  for a time series x of length n. The moving window has k2 (odd) points
*/
/*   (except for ends where as many observations as possible are used). */

/*   xsumm is (5 x n): (m is number of points in the window for column i) 
*/

/*       Row 1 is largest x in [u4, u4 + 1.5 (u4-l4)] or u4 if none */
/*       Row 2 is upper fourth (median of largest (m+1)/2) */
/*       Row 3 is median */
/*       Row 4 is lower fourth (median of smallest (m+1)/2) */
/*       Row 5 is smallest x in [l4, l4 - 1.5 (u4-l4)] or l4 if none */

/*   indout and outval are arrays (integer and real) of indices and */
/*   corresponding values outside of [l4 - 1.5 (u4-l4), u4 + 1.5 (u4-l4)].
 */

/*   indout and outval can be no longer than n elements */

/*   The arrays ii and jj are integer work arrays of length k2. */

/* ***********************************************************************
 */


    /* Parameter adjustments */
    xsumm -= 6;
    --x;
    --jj;
    --ii;
    --outval;

    /* Function Body */
    k = (*k2 - 1) / 2;
    ii[1] = 1.;
    jj[1] = 1.;
    nl = 1;
    summ5_(&ii[1], &jj[1], &nl, &c__1, &x[1], &xsumm[6], nouts, indout, &
	    outval[1]);

    i__1 = k + 1;
    for (m = 2; m <= i__1; ++m) {
	medadd_(&ii[1], &jj[1], &nl, &x[1]);
	++nl;
	medadd_(&ii[1], &jj[1], &nl, &x[1]);
	++nl;
/* L10: */
	summ5_(&ii[1], &jj[1], &nl, &m, &x[1], &xsumm[m * 5 + 1], nouts, 
		indout, &outval[1]);
    }

    i__1 = *n - k;
    for (m = k + 2; m <= i__1; ++m) {
	meddel_(&ii[1], &jj[1], k2);
	i__2 = *k2 - 1;
	medadd_(&ii[1], &jj[1], &i__2, &x[1]);
/* L20: */
	summ5_(&ii[1], &jj[1], k2, &m, &x[1], &xsumm[m * 5 + 1], nouts, 
		indout, &outval[1]);
    }

    nl = *k2;
    i__1 = *n;
    for (m = *n - k + 1; m <= i__1; ++m) {
	meddel_(&ii[1], &jj[1], &nl);
	--nl;
	meddel_(&ii[1], &jj[1], &nl);
	--nl;
/* L30: */
	summ5_(&ii[1], &jj[1], &nl, &m, &x[1], &xsumm[m * 5 + 1], nouts, 
		indout, &outval[1]);
    }

    return 0;
} /* movbox_ */

/* ********************************************************************* */

/* Subroutine */ int summ5_(ii, jj, n, npos, x, xsumm, nouts, indout, outval)
doublereal *ii, *jj;
integer *n, *npos;
doublereal *x, *xsumm;
integer *nouts, *indout;
doublereal *outval;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal xmid, flow, d__;
    static integer i__, m;
    extern /* Subroutine */ int medii_();
    static integer nn;
    static doublereal fup;


/*   Subroutine to */


/* ********************************************************************* 
*/



    /* Parameter adjustments */
    --jj;
    --ii;
    --x;
    --xsumm;
    --indout;
    --outval;

    /* Function Body */
    m = (*n + 1) / 2;
    medii_(&ii[1], n, &x[1], &xsumm[3]);
    medii_(&ii[*n - m + 1], &m, &x[1], &xsumm[2]);
    medii_(&ii[1], &m, &x[1], &xsumm[4]);

    d__ = (xsumm[2] - xsumm[4]) * 1e5;
    fup = xsumm[2] + d__;
    flow = xsumm[4] - d__;
    xsumm[1] = xsumm[2];
    xsumm[5] = xsumm[4];

    nn = 0;
    i__1 = *n;
    for (i__ = m; i__ <= i__1; ++i__) {
	if (x[(integer) ii[i__]] > fup) {
	    goto L20;
	}
/* L10: */
	if (x[(integer) ii[i__]] > xsumm[2]) {
	    nn = i__;
	}
    }
L20:
    if (nn != 0) {
	xsumm[1] = x[(integer) ii[nn]];
    }

    nn = 0;
    for (i__ = m; i__ >= 1; --i__) {
	if (x[(integer) ii[i__]] < flow) {
	    goto L40;
	}
/* L30: */
	if (x[(integer) ii[i__]] < xsumm[4]) {
	    nn = i__;
	}
    }
L40:
    if (nn != 0) {
	xsumm[5] = x[(integer) ii[nn]];
    }

    xmid = x[*npos];
    if (xmid < flow || xmid > fup) {
	++(*nouts);
	indout[*nouts] = *npos;
	outval[*nouts] = xmid;
    }

    return 0;
} /* summ5_ */


/* ********************************************************************* */

/* Subroutine */ int medii_(ii, n, x, xmed)
doublereal *ii;
integer *n;
doublereal *x, *xmed;
{

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --x;
    --ii;

    /* Function Body */
    *xmed = (x[(integer) ii[(*n + 1) / 2]] + x[(integer) ii[(*n + 2) / 2]]) / 
	    2;

    return 0;
} /* medii_ */

/* *********************************************************************** */

/* Subroutine */ int movmed_(x, n, k2, ii, jj, xmed)
doublereal *x;
integer *n, *k2;
doublereal *ii, *jj, *xmed;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, m;
    extern /* Subroutine */ int medadd_(), meddel_();
    static integer nl;


/*   Find moving medians (xmed(1),...xmed(n)) of x(1),...,x(n) where */
/*   (except at ends) medians are of k2 x's. At end points, medians */
/*   are of as many points as possible (for example, xmed(1) and xmed(n) 
*/
/*   are medians of 1 point). */

/*   The arrays ii and jj are integer work arrays of length k2. */

/* ***********************************************************************
 */


    /* Parameter adjustments */
    --xmed;
    --x;
    --jj;
    --ii;

    /* Function Body */
    k = (*k2 - 1) / 2;
    ii[1] = 1.;
    jj[1] = 1.;
    xmed[1] = x[1];

    i__1 = k + 1;
    for (m = 2; m <= i__1; ++m) {
	i__2 = (m << 1) - 3;
	medadd_(&ii[1], &jj[1], &i__2, &x[1]);
	i__2 = (m << 1) - 2;
	medadd_(&ii[1], &jj[1], &i__2, &x[1]);
/* L10: */
	xmed[m] = x[(integer) ii[m]];
    }

    i__1 = *n - k;
    for (m = k + 2; m <= i__1; ++m) {
	meddel_(&ii[1], &jj[1], k2);
	i__2 = *k2 - 1;
	medadd_(&ii[1], &jj[1], &i__2, &x[1]);
/* L20: */
	xmed[m] = x[(integer) ii[k + 1]];
    }

    nl = *k2;
    i__1 = *n;
    for (m = *n - k + 1; m <= i__1; ++m) {
	meddel_(&ii[1], &jj[1], &nl);
	--nl;
	meddel_(&ii[1], &jj[1], &nl);
	--nl;
/* L30: */
	xmed[m] = x[(integer) ii[(nl + 1) / 2]];
    }

    return 0;
} /* movmed_ */


/* *********************************************************************** */

/* Subroutine */ int movord_(x, n, k2, nord, ii, jj, xord)
doublereal *x;
integer *n, *k2, *nord;
doublereal *ii, *jj, *xord;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, m;
    extern /* Subroutine */ int medadd_(), meddel_();
    static integer nl;


/*   Find moving nord order statistic (xord(1),...xord(n)) of */
/*   x(1),...,x(n) where (except at ends) order statistics are of k2 x's. 
*/
/*   At end points, they are of as many points as possible (for example, 
*/
/*   xord(1) and xord(n) are for 1 point). */

/*   The arrays ii and jj are integer work arrays of length k2. */

/* ***********************************************************************
 */


    /* Parameter adjustments */
    --xord;
    --x;
    --jj;
    --ii;

    /* Function Body */
    k = (*k2 - 1) / 2;
    ii[1] = 1.;
    jj[1] = 1.;
    xord[1] = x[1];

    i__1 = k + 1;
    for (m = 2; m <= i__1; ++m) {
	i__2 = (m << 1) - 3;
	medadd_(&ii[1], &jj[1], &i__2, &x[1]);
	i__2 = (m << 1) - 2;
	medadd_(&ii[1], &jj[1], &i__2, &x[1]);
	if (m == k + 1) {
	    xord[m] = x[(integer) ii[*nord]];
	}
	if (m < k + 1) {
	    xord[m] = x[(integer) ii[m]];
	}
/* L10: */
    }

    i__1 = *n - k;
    for (m = k + 2; m <= i__1; ++m) {
	meddel_(&ii[1], &jj[1], k2);
	i__2 = *k2 - 1;
	medadd_(&ii[1], &jj[1], &i__2, &x[1]);
/* L20: */
	xord[m] = x[(integer) ii[*nord]];
    }

    nl = *k2;
    i__1 = *n;
    for (m = *n - k + 1; m <= i__1; ++m) {
	meddel_(&ii[1], &jj[1], &nl);
	--nl;
	meddel_(&ii[1], &jj[1], &nl);
	--nl;
/* L30: */
	xord[m] = x[(integer) ii[(nl + 1) / 2]];
    }

    return 0;
} /* movord_ */


/* *********************************************************************** */

/* Subroutine */ int meddel_(ii, jj, nl)
doublereal *ii, *jj;
integer *nl;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer noff, i__, i1, nr;

/*   Delete a point in moving median */

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --jj;
    --ii;

    /* Function Body */
    nr = (integer) jj[1];
    noff = (integer) (ii[nr] - 1);

    if (nr < *nl) {
	i__1 = *nl;
	for (i__ = nr + 1; i__ <= i__1; ++i__) {
	    i1 = (integer) (ii[i__] - noff);
	    --jj[i1];
/* L10: */
	    ii[i__ - 1] = ii[i__];
	}
    }

    i__1 = *nl;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L20: */
	jj[i__ - 1] = jj[i__];
    }

    return 0;
} /* meddel_ */


/* *********************************************************************** */

/* Subroutine */ int medadd_(ii, jj, nl, x)
doublereal *ii, *jj;
integer *nl;
doublereal *x;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer noff;
    static doublereal xnew;
    static integer i__;
    extern integer bsrch_();
    static integer i1, ns;


/*   Add a point in moving median */

/* ***********************************************************************
 */



    /* Parameter adjustments */
    --x;
    --jj;
    --ii;

    /* Function Body */
    xnew = x[(integer) (ii[(integer) jj[*nl]] + 1)];
    noff = (integer) (ii[(integer) jj[1]] - 1);

    ns = bsrch_(&x[1], &ii[1], nl, &xnew);

    if (ns <= *nl) {
	i__1 = *nl;
	for (i__ = ns; i__ <= i__1; ++i__) {
	    i1 = (integer) (ii[i__] - noff);
/* L30: */
	    ++jj[i1];
	}
	i__1 = ns;
	for (i__ = *nl; i__ >= i__1; --i__) {
/* L40: */
	    ii[i__ + 1] = ii[i__];
	}
    }

    jj[*nl + 1] = (doublereal) ns;
    ii[ns] = (doublereal) (*nl + noff + 1);

    return 0;
} /* medadd_ */


/***************************************************************************/

integer bsrch_(x, ii, n, xnew)
doublereal *x, *ii;
integer *n;
doublereal *xnew;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer nl, nr, mid;


/*   Find location of xnew in x(ii(1)) <= x(ii(2)) <= ... <= x(ii(nl)) */

/************************************************************************
***/


    /* Parameter adjustments */
    --ii;
    --x;

    /* Function Body */
    nl = 0;
    nr = *n + 1;

L5:
    if (nr - nl == 1) {
	goto L99;
    }

    mid = (nl + nr) / 2;
    if (*xnew > x[(integer) ii[mid]]) {
	nl = mid;
	goto L5;
    } else {
	nr = mid;
	goto L5;
    }

L99:
    ret_val = nr;
    return ret_val;
} /* bsrch_ */


/***************************************************************************/

/* Subroutine */ int pacf_(x, n, m, e, f, temp, theta)
doublereal *x;
integer *n, *m;
doublereal *e, *f, *temp, *theta;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal xbar, part;
    static integer i__;
    extern /* Subroutine */ int crlag_();
    static doublereal apart;
    extern /* Subroutine */ int accvec_();
    extern doublereal dtprod_();
    extern /* Subroutine */ int movvec_(), submns_();
    static integer npm;


/*   x of length n */
/*   e, f, and temp of length n+m */
/*   theta of length m */

/************************************************************************
***/



    /* Parameter adjustments */
    --theta;
    --temp;
    --f;
    --e;
    --x;

    /* Function Body */
    submns_(&x[1], n, &xbar, &e[1]);
    npm = *n + *m;
    i__1 = npm;
    for (i__ = *n + 1; i__ <= i__1; ++i__) {
/* L10: */
	e[i__] = (float)0.;
    }
    crlag_(&e[1], &npm, &f[1]);


    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	part = dtprod_(&f[1], &e[1], &npm) / dtprod_(&f[1], &f[1], &npm);
	apart = -part;
	movvec_(&e[1], &npm, &temp[1]);
	accvec_(&temp[1], &apart, &f[1], &npm, &e[1]);
	accvec_(&f[1], &apart, &temp[1], &npm, &f[1]);
	crlag_(&f[1], &npm, &f[1]);
	theta[i__] = part;
/* L30: */
    }


    return 0;
} /* pacf_ */


/* *********************************************************************** */

/* Subroutine */ int accvec_(x, a, y, n, z__)
doublereal *x, *a, *y;
integer *n;
doublereal *z__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*   Find z=x+a*y */

/************************************************************************
**/



    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	z__[i__] = x[i__] + *a * y[i__];
    }

    return 0;
} /* accvec_ */


/***************************************************************************/

doublereal dtprod_(x, y, n)
doublereal *x, *y;
integer *n;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__;


/************************************************************************
***/


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    ret_val = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	ret_val += x[i__] * y[i__];
    }

    return ret_val;
} /* dtprod_ */


/****************************************************************************/

/* Subroutine */ int movvec_(x, n, y)
doublereal *x;
integer *n;
doublereal *y;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/************************************************************************
*****/


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = x[i__];
    }

    return 0;
} /* movvec_ */


/****************************************************************************
*/

/* Subroutine */ int crlag_(x, n, y)
doublereal *x;
integer *n;
doublereal *y;
{
    static doublereal temp;
    static integer i__;



/************************************************************************
*****/



    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    temp = x[*n];
    for (i__ = *n; i__ >= 2; --i__) {
/* L10: */
	y[i__] = x[i__ - 1];
    }
    y[1] = temp;

    return 0;
} /* crlag_ */


/****************************************************************************
*/

/* Subroutine */ int submns_(x, n, xbar, y)
doublereal *x;
integer *n;
doublereal *xbar, *y;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/************************************************************************
*****/



    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    *xbar = (float)0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	*xbar += x[i__];
    }
    *xbar /= (doublereal) (*n);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	y[i__] = x[i__] - *xbar;
    }

    return 0;
} /* submns_ */


/* *********************************************************************** */

/* Subroutine */ int partar_(alpha, p)
doublereal *alpha;
integer *p;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal temp, c__;
    static integer j, k;


/* ***********************************************************************
 */


    /* Parameter adjustments */
    --alpha;

    /* Function Body */
    if (*p > 1) {
	i__1 = *p;
	for (j = 2; j <= i__1; ++j) {
	    c__ = alpha[j];
	    i__2 = j / 2;
	    for (k = 1; k <= i__2; ++k) {
		temp = alpha[k] - c__ * alpha[j - k];
		alpha[j - k] -= c__ * alpha[k];
/* L10: */
		alpha[k] = temp;
	    }
/* L20: */
	}
    }

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/* L30: */
	alpha[j] = -alpha[j];
    }


    return 0;
} /* partar_ */


/* ********************************************************************** */

/* Subroutine */ int poly_(coeffs, degree, x, n, fx)
doublereal *coeffs;
integer *degree;
doublereal *x;
integer *n;
doublereal *fx;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal dc, xi;


/*   Evaluate polynomial */

/*   f(x)=coeffs(1) + sum(j=1,degree) coeffs(j+1)*x^j */

/*   for x(1),...,x(n) */

/* ***********************************************************************
 */



    /* Parameter adjustments */
    --coeffs;
    --fx;
    --x;

    /* Function Body */
    dc = coeffs[*degree + 1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = x[i__];
	c__ = dc;
	for (j = *degree; j >= 1; --j) {
/* L10: */
	    c__ = c__ * xi + coeffs[j];
	}
/* L20: */
	fx[i__] = c__;
    }


    return 0;
} /* poly_ */


/* ******************************************************************** */

/* Subroutine */ int rtpoly_(r__, np, wk, a)
doublereal *r__;
integer *np;
doublereal *wk, *a;
{
    /* Initialized data */

    static doublereal cone[2] = { 1.,0. };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int ccacc_();
    static integer i__, j;
    extern /* Subroutine */ int ccneg_(), ccdiv_(), ccmul_();
    extern doublereal ccreal_();
    static doublereal zi[2];


/*   Subroutine to find the coefficients a(1),...,a(np) of the polynomial 
*/

/*      1 + Sum(j=1,np) a(j)z**j */

/*   given its roots r(1),...,r(np) */

/*   Input:  r (complex), np */

/*   Output: a (real) */

/*   Work:   wk (complex array of length np) */

/* ********************************************************************** 
*/

    /* Parameter adjustments */
    --a;
    wk -= 3;
    r__ -= 3;

    /* Function Body */

    ccdiv_(cone, &r__[3], &wk[3]);
    ccneg_(&wk[3], &wk[3]);
    if (*np == 1) {
	goto L30;
    }

    i__1 = *np;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ccdiv_(cone, &r__[(i__ << 1) + 1], zi);
	ccneg_(zi, zi);
	ccmul_(zi, &wk[(i__ - 1 << 1) + 1], &wk[(i__ << 1) + 1]);
	if (i__ == 2) {
	    goto L15;
	}
	for (j = i__ - 1; j >= 2; --j) {
/* L10: */
	    ccacc_(&wk[(j << 1) + 1], zi, &wk[(j - 1 << 1) + 1], &wk[(j << 1) 
		    + 1]);
	}
L15:
	ccacc_(&wk[3], zi, cone, &wk[3]);
/* L20: */
    }

L30:
    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	a[i__] = ccreal_(&wk[(i__ << 1) + 1]);
    }

    return 0;
} /* rtpoly_ */


/* ********************************************************************** */

/* Subroutine */ int ccadd_(x, y, z__)
doublereal *x, *y, *z__;
{

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    z__[1] = x[1] + y[1];
    z__[2] = x[2] + y[2];
    return 0;
} /* ccadd_ */


/* ********************************************************************** */

/* Subroutine */ int ccsub_(x, y, z__)
doublereal *x, *y, *z__;
{

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    z__[1] = x[1] - y[1];
    z__[2] = x[2] - y[2];
    return 0;
} /* ccsub_ */


/* ********************************************************************** */

/* Subroutine */ int ccneg_(x, y)
doublereal *x, *y;
{

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    y[1] = -x[1];
    y[2] = -x[2];
    return 0;
} /* ccneg_ */


/* ********************************************************************** */

/* Subroutine */ int ccmul_(x, y, z__)
doublereal *x, *y, *z__;
{
    static doublereal temp;


/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    temp = x[1] * y[1] - x[2] * y[2];
    z__[2] = x[1] * y[2] + x[2] * y[1];
    z__[1] = temp;
    return 0;
} /* ccmul_ */


/* ********************************************************************** */

/* Subroutine */ int ccdiv_(x, y, z__)
doublereal *x, *y, *z__;
{
    static doublereal temp, y2;


/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    y2 = y[1] * y[1] + y[2] * y[2];
    temp = (x[1] * y[1] + x[2] * y[2]) / y2;
    z__[2] = (x[2] * y[1] - x[1] * y[2]) / y2;
    z__[1] = temp;
    return 0;
} /* ccdiv_ */


/* ********************************************************************** */

/* Subroutine */ int ccacc_(x, a, y, z__)
doublereal *x, *a, *y, *z__;
{
    static doublereal temp;


/* ********************************************************************** 
*/


/*   z = x + a * y */


    /* Parameter adjustments */
    --z__;
    --y;
    --a;
    --x;

    /* Function Body */
    temp = x[1] + a[1] * y[1] - a[2] * y[2];
    z__[2] = x[2] + a[2] * y[1] + a[1] * y[2];
    z__[1] = temp;
    return 0;
} /* ccacc_ */


/* ********************************************************************** */

/* Subroutine */ int cccopy_(x, y)
doublereal *x, *y;
{

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    x[1] = y[1];
    x[2] = y[2];
    return 0;
} /* cccopy_ */


/* ********************************************************************** */

doublereal ccimag_(x)
doublereal *x;
{
    /* System generated locals */
    doublereal ret_val;


/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = x[2];
    return ret_val;
} /* ccimag_ */


/* ********************************************************************** */

doublereal ccreal_(x)
doublereal *x;
{
    /* System generated locals */
    doublereal ret_val;


/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = x[1];
    return ret_val;
} /* ccreal_ */


/* ********************************************************************** */

/* Subroutine */ int ccplx_(a, b, z__)
doublereal *a, *b, *z__;
{

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --z__;

    /* Function Body */
    z__[1] = *a;
    z__[2] = *b;
    return 0;
} /* ccplx_ */


/* ********************************************************************** */

doublereal ccabs_(x)
doublereal *x;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt();


/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = sqrt(x[1] * x[1] + x[2] * x[2]);
    return ret_val;
} /* ccabs_ */


/* ********************************************************************** */

/* Subroutine */ int crmul_(a, x, y)
doublereal *a, *x, *y;
{

/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    y[1] = *a * x[1];
    y[2] = *a * x[2];
    return 0;
} /* crmul_ */


/* ********************************************************************** */

/* Subroutine */ int ccsqrt_(x, y)
doublereal *x, *y;
{
    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal r__, w, xx, yy;


/* ********************************************************************** 
*/


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    if (x[1] == (float)0. && x[2] == (float)0.) {
	y[1] = 0.;
	y[2] = 0.;
	return 0;
    }
    xx = abs(x[1]);
    yy = abs(x[2]);
    if (xx >= yy) {
	r__ = yy / xx;
	w = sqrt(xx) * sqrt((sqrt(r__ * r__ + (float)1.) + (float)1.) * (
		float).5);
    } else {
	r__ = xx / yy;
	w = sqrt(yy) * sqrt((r__ + sqrt(r__ * r__ + (float)1.)) * (float).5);
    }
    if (x[1] >= (float)0.) {
	y[1] = w;
	y[2] = x[2] / (w * (float)2.);
    } else {
	y[2] = -w;
	if (x[2] >= (float)0.) {
	    y[2] = w;
	}
	y[1] = x[2] / (y[2] * (float)2.);
    }
    return 0;
} /* ccsqrt_ */


/* ***************************************************************** */

/* Subroutine */ int polyrt_(ra, np, maxit, eps, a, a1, roots, ier)
doublereal *ra;
integer *np, *maxit;
real *eps;
doublereal *a, *a1, *roots;
integer *ier;
{
    /* Initialized data */

    static doublereal cone[2] = { 1.,0. };
    static doublereal czero[2] = { 0.,0. };

    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int ccacc_();
    static doublereal b[2], c__[2];
    extern /* Subroutine */ int root1_();
    static integer i__, j;
    static doublereal x[2];
    extern /* Subroutine */ int ccplx_();
    extern doublereal ccimag_();
    static integer jj;
    extern doublereal ccreal_();
    extern /* Subroutine */ int cccopy_();


/*   calculate roots of a(1)+a(2)z+a(3)z**2+...+a(np+1)z**np */

/* ***************************************************************** */

    /* Parameter adjustments */
    roots -= 3;
    a1 -= 3;
    a -= 3;
    --ra;

    /* Function Body */

    cccopy_(&a[3], cone);
    cccopy_(&a1[3], cone);

    i__1 = *np + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	a[(i__ << 1) + 1] = ra[i__ - 1];
	a[(i__ << 1) + 2] = (float)0.;
/* L1: */
	cccopy_(&a1[(i__ << 1) + 1], &a[(i__ << 1) + 1]);
    }
    for (j = *np; j >= 1; --j) {
	cccopy_(x, czero);
	root1_(&a1[3], &j, x, eps, &c__0, maxit, ier);
	if (*ier == 1) {
	    *ier = 2;
	    goto L99;
	}
/* Computing 2nd power */
	r__1 = *eps;
	if ((d__1 = ccimag_(x), abs(d__1)) <= r__1 * r__1 * (float)2. * (d__2 
		= ccreal_(x), abs(d__2))) {
	    d__1 = ccreal_(x);
	    ccplx_(&d__1, &c_b307, x);
	}
	cccopy_(&roots[(j << 1) + 1], x);
	cccopy_(b, &a1[(j + 1 << 1) + 1]);
	for (jj = j; jj >= 1; --jj) {
	    cccopy_(c__, &a1[(jj << 1) + 1]);
	    cccopy_(&a1[(jj << 1) + 1], b);
/* L2: */
	    ccacc_(c__, b, x, b);
	}
/* L3: */
    }
    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	root1_(&a[3], np, &roots[(j << 1) + 1], eps, &c__1, maxit, ier);
/* L4: */
	if (*ier == 1) {
	    goto L99;
	}
    }
    i__1 = *np;
    for (j = 2; j <= i__1; ++j) {
	cccopy_(x, &roots[(j << 1) + 1]);
	for (i__ = j - 1; i__ >= 1; --i__) {
	    if (ccreal_(&roots[(i__ << 1) + 1]) <= ccreal_(x)) {
		goto L10;
	    }
	    cccopy_(&roots[(i__ + 1 << 1) + 1], &roots[(i__ << 1) + 1]);
/* L5: */
	}
	i__ = 0;
L10:
	cccopy_(&roots[(i__ + 1 << 1) + 1], x);
/* L6: */
    }
L99:
    return 0;
} /* polyrt_ */


/* ******************************************************************** */

/* Subroutine */ int root1_(a, np, x, eps, iopt, maxit, ier)
doublereal *a;
integer *np;
doublereal *x;
real *eps;
integer *iopt, *maxit, *ier;
{
    /* Initialized data */

    static doublereal czero[2] = { 0.,0. };

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd();

    /* Local variables */
    static integer npol;
    extern /* Subroutine */ int ccacc_(), ccadd_();
    static doublereal b[2], d__[2], f[2], g[2], h__[2];
    static integer j;
    extern doublereal ccabs_();
    extern /* Subroutine */ int ccdiv_(), ccsub_(), ccmul_(), ccplx_();
    static doublereal dxold;
    extern /* Subroutine */ int crmul_();
    static doublereal g2[2], t1[2], t2[2], x1[2], gm[2], gp[2], dx[2];
    static integer it;
    static doublereal sq[2];
    extern /* Subroutine */ int cccopy_(), ccsqrt_();
    static doublereal cdx;


/*   calculate one root. */

/* ******************************************************************** */

    /* Parameter adjustments */
    --x;
    a -= 3;

    /* Function Body */

    *ier = 0;
    if (*iopt == 1) {
	dxold = ccabs_(&x[1]);
	npol = 0;
    }
    i__1 = *maxit;
    for (it = 1; it <= i__1; ++it) {
	cccopy_(b, &a[(*np + 1 << 1) + 1]);
	cccopy_(d__, czero);
	cccopy_(f, czero);
	for (j = *np; j >= 1; --j) {
	    ccacc_(d__, &x[1], f, f);
	    ccacc_(b, &x[1], d__, d__);
/* L1: */
	    ccacc_(&a[(j << 1) + 1], &x[1], b, b);
	}
	if (ccabs_(b) <= (float)1e-15) {
	    cccopy_(dx, czero);
	} else if (ccabs_(d__) <= (float)1e-15 && ccabs_(f) <= (float)1e-15) {
	    ccdiv_(b, &a[(*np + 1 << 1) + 1], t1);
	    d__2 = ccabs_(t1);
	    d__3 = (doublereal) ((float)1. / *np);
	    d__1 = pow_dd(&d__2, &d__3);
	    ccplx_(&d__1, &c_b307, dx);
	} else {
	    ccdiv_(d__, b, g);
	    ccmul_(g, g, g2);
	    ccplx_(&c_b318, &c_b307, t1);
	    ccdiv_(f, b, t2);
	    ccacc_(g2, t1, t2, h__);
	    d__1 = (doublereal) (*np);
	    crmul_(&d__1, h__, t1);
	    ccsub_(t1, g2, t2);
	    d__1 = (doublereal) (*np - 1);
	    crmul_(&d__1, t2, t1);
	    ccsqrt_(t1, sq);
	    ccadd_(g, sq, gp);
	    ccsub_(g, sq, gm);
	    if (ccabs_(gp) < ccabs_(gm)) {
		cccopy_(gp, gm);
	    }
	    d__1 = (doublereal) (*np);
	    ccplx_(&d__1, &c_b307, t1);
	    ccdiv_(t1, gp, dx);
	}
	ccsub_(&x[1], dx, x1);
	if (x[1] == x1[0] && x[2] == x1[1]) {
	    return 0;
	}
	cccopy_(&x[1], x1);
	if (*iopt == 1) {
	    ++npol;
	    cdx = ccabs_(dx);
	    if (npol > 9 && cdx >= dxold) {
		return 0;
	    }
	    dxold = cdx;
	} else {
	    if (ccabs_(dx) <= *eps * ccabs_(&x[1])) {
		return 0;
	    }
	}
/* L2: */
    }
    *ier = 1;
    return 0;
} /* root1_ */


/* ******************************************************* */

/* Subroutine */ int schur_(alph, ndim, np, a)
doublereal *alph;
integer *ndim, *np;
doublereal *a;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer j, k, l;
    static doublereal d1, d2;


/*   subroutine to form the schur matrix of the coefficients */
/*   alph(1),...,alph(np) (see pagano : when is an */
/*   autoregressive scheme stationary ) */

/*   input : */
/*           ndim : dimension of matrix a in main program */
/*           np, alph(1),...,alph(np) */

/*   output : */
/*           a */

/*   subroutines called : none */

/* ******************************************************* */


    /* Parameter adjustments */
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --alph;

    /* Function Body */
    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (k = j; k <= i__2; ++k) {
	    c__ = 0.;
	    i__3 = j;
	    for (l = 1; l <= i__3; ++l) {
		if (l == j) {
		    goto L2;
		}
		d1 = alph[j - l];
		goto L3;
L2:
		d1 = 1.;
L3:
		if (l == k) {
		    goto L4;
		}
		d2 = alph[k - l];
		goto L6;
L4:
		d2 = 1.;
L6:
		c__ = c__ + d1 * d2 - alph[*np + l - j] * alph[*np + l - k];
	    }
	    a[j + k * a_dim1] = c__;
/* L1: */
	    a[k + j * a_dim1] = a[j + k * a_dim1];
	}
    }
    return 0;
} /* schur_ */


/* ****************************************************************** */

/* Subroutine */ int marq_(y, n, ndim, maxit, eps, nords, coeffs, nt, npnt, 
	ntot, e, e1, xx, x, a, as, lags, ier, rvar, d__)
doublereal *y;
integer *n, *ndim, *maxit;
doublereal *eps;
integer *nords;
doublereal *coeffs;
integer *nt, *npnt, *ntot;
doublereal *e, *e1, *xx, *x, *a, *as;
integer *lags, *ier;
doublereal *rvar, *d__;
{
    /* System generated locals */
    integer x_dim1, x_offset, a_dim1, a_offset, as_dim1, as_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal beta[100], beta1[100], g[100], h__[100];
    static integer i__, j;
    static doublereal alpha[100], delta, a1[100];
    extern /* Subroutine */ int swpk12_();
    static doublereal f2, cc, gg;
    static integer nl[100];
    static doublereal lp[100], lq[100];
    static integer it;
    static doublereal ss;
    extern doublereal seaslk_(), dtprod_();
    static doublereal gg1, gg2, ss1;
    extern /* Subroutine */ int movxyr_();
    static doublereal ppi;
    static integer ntp;


/*   subroutine to implement the algorithm on page 504 of box and jenkins 
*/
/*   except in the language of timeslab. */

/* ****************************************************************** */


/*   ppi, f2, and delta have the meaning in B+J */

    /* Parameter adjustments */
    --y;
    as_dim1 = *ndim;
    as_offset = as_dim1 + 1;
    as -= as_offset;
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --nords;
    --coeffs;
    x_dim1 = *npnt;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --e;
    --e1;
    --xx;
    --d__;

    /* Function Body */
    ppi = 0.;
    f2 = 2.;
    delta = 0.;
    *npnt = *n + *nt;
    ntp = *ntot + 1;


/*   start iterations: */


    i__1 = *maxit;
    for (it = 1; it <= i__1; ++it) {

	ss = seaslk_(&y[1], n, &nords[1], &coeffs[1], lags, lp, lq, alpha, 
		beta, &e[1], a1, nl, &xx[1], nt);

/*   Get x matrix: */

	i__2 = *ntot;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    coeffs[i__] += delta;
	    ss1 = seaslk_(&y[1], n, &nords[1], &coeffs[1], lags, lp, lq, 
		    alpha, beta, &e1[1], a1, nl, &xx[1], nt);
	    coeffs[i__] -= delta;
	    i__3 = *npnt;
	    for (j = 1; j <= i__3; ++j) {
/* L20: */
		x[j + i__ * x_dim1] = (e[j] - e1[j]) / delta;
	    }
/* L10: */
	}

/*   Get a, d, g: */

	i__2 = *ntot;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    for (j = 1; j <= i__3; ++j) {
		a[i__ + j * a_dim1] = dtprod_(&x[i__ * x_dim1 + 1], &x[j * 
			x_dim1 + 1], npnt);
/* L30: */
		a[j + i__ * a_dim1] = a[i__ + j * a_dim1];
	    }
	    g[i__ - 1] = dtprod_(&x[i__ * x_dim1 + 1], &e[1], npnt);
	    d__[i__] = sqrt(a[i__ + i__ * a_dim1]);
/* L15: */
	}
	if (*maxit == 1) {
	    *ier = 1;
	    goto L200;
	}


L25:
	i__2 = *ntot;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    as[ntp + i__ * as_dim1] = g[i__ - 1] / d__[i__];
	    as[i__ + ntp * as_dim1] = as[ntp + i__ * as_dim1];
	    i__3 = *ntot;
	    for (j = 1; j <= i__3; ++j) {
		as[i__ + j * as_dim1] = a[i__ + j * a_dim1] / (d__[i__] * d__[
			j]);
/* L40: */
	    }
	}
	i__3 = *ntot;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L50: */
	    as[i__ + i__ * as_dim1] += ppi;
	}
	as[ntp + ntp * as_dim1] = (float)1.;

/*   solve system: */

	swpk12_(&as[as_offset], &ntp, &ntp, &c__1, ntot, ier);

/*   check for errors in solving system: */

	if (*ier != 0) {
	    *ier = 2;
	    goto L99;
	}

/*   update and get new sum of squares: */

	i__3 = *ntot;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    h__[i__ - 1] = as[i__ + ntp * as_dim1] / d__[i__];
/* L60: */
	    beta1[i__ - 1] = coeffs[i__] + h__[i__ - 1];
	}
	ss1 = seaslk_(&y[1], n, &nords[1], beta1, lags, lp, lq, alpha, beta, &
		e1[1], a1, nl, &xx[1], nt);

/*   check for decrease in sum of squares: */

	if (ss1 >= ss) {

/*   no: */

	    if (ppi > (float)200.) {
		*ier = 3;
		goto L99;
	    }
	    ppi *= f2;
	    goto L25;
	}

/*   yes: */

/*   check for convergence: */

	i__3 = *ntot;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L70: */
	    if ((d__1 = h__[i__ - 1], abs(d__1)) > *eps * (d__2 = coeffs[i__],
		     abs(d__2))) {
		goto L75;
	    }
	}

/*   yes: */

	movxyr_(&coeffs[1], beta1, ntot);
	*ier = 0;
	goto L200;

/*   no: */

L75:
	movxyr_(&coeffs[1], beta1, ntot);
	if (ppi < (float)1e-20) {
	    *ier = 4;
	    goto L99;
	}
	ppi /= f2;

/*   end iteration: */

/* L100: */
    }

/*   ran out of iterations: */

    *ier = 1;
L200:
L99:

/*   get standard errors: */

    *rvar = ss / *n;
    swpk12_(&a[a_offset], ndim, ntot, &c__1, ntot, ier);
    i__1 = *ntot;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	d__[i__] = sqrt(*rvar * a[i__ + i__ * a_dim1]);
    }

/*   get constant: */

    if (nords[5] == 1) {

	cc = coeffs[*ntot];
	gg1 = (float)1.;
	if (nords[1] > 0) {
	    i__1 = nords[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
		gg1 += coeffs[i__];
	    }
	}
	gg2 = (float)1.;
	if (nords[2] > 0) {
	    i__1 = nords[2];
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
		gg2 += coeffs[nords[1] + i__];
	    }
	}
	gg = gg1 * gg2;
	cc *= gg;
	coeffs[*ntot] = cc;
	d__[*ntot] = sqrt(*rvar * d__[*ntot] * gg * gg);
    }

    return 0;
} /* marq_ */



/* ********************************************************************* */

doublereal seaslk_(x, n, nords, coeffs, lags, lp, lq, alpha, beta, e, a, nl, 
	xx, nt)
doublereal *x;
integer *n, *nords;
doublereal *coeffs;
integer *lags;
doublereal *lp, *lq, *alpha, *beta, *e, *a;
integer *nl;
doublereal *xx;
integer *nt;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer maxp, maxq;
    static doublereal c__;
    static integer i__, j, k;
    static doublereal c1;
    extern /* Subroutine */ int convt_();
    static integer ii, ij;
    static doublereal ss, amu;
    static integer nlp, nlq;



/*   Function to evaluate a seasonal ARMA sum of squares. */

/* ********************************************************************** 
*/



    /* Parameter adjustments */
    --xx;
    --a;
    --e;
    --beta;
    --alpha;
    --lq;
    --lp;
    --coeffs;
    --x;

    /* Function Body */
    convt_(nords, &coeffs[1], lags, &maxp, &maxq, &nlp, &nlq, &lp[1], &lq[1], 
	    &alpha[1], &beta[1], &a[1], nl, &amu);



/*   find eta: */

    if (*nt == 0) {
	goto L715;
    }
    if (maxq > 0) {
	i__1 = *n - maxp + maxq;
	for (i__ = *n - maxp + 1; i__ <= i__1; ++i__) {
/* L310: */
	    e[i__] = 0.;
	}
    }
    for (i__ = *n - maxp; i__ >= 1; --i__) {
	c__ = x[i__] - amu;
	if (nlp > 0) {
	    i__1 = nlp;
	    for (j = 1; j <= i__1; ++j) {
/* L320: */
		c__ += alpha[j] * (x[(integer) (i__ + lp[j])] - amu);
	    }
	}
	if (nlq > 0) {
	    i__1 = nlq;
	    for (k = 1; k <= i__1; ++k) {
/* L340: */
		c__ -= beta[k] * e[(integer) (i__ + lq[k])];
	    }
	}
/* L360: */
	e[i__] = c__;
    }

/*   now get x(0),x(-1),... (x(i)=xx(nt+i), i=0,...,1-nt) */

    i__1 = -(*nt) + 1;
    for (i__ = 0; i__ >= i__1; --i__) {
	c__ = (float)0.;
	if (nlq > 0) {
	    i__2 = nlq;
	    for (k = 1; k <= i__2; ++k) {
		ii = (integer) (i__ + lq[k]);
		if (ii >= 1) {
		    c__ += beta[k] * e[ii];
		}
/* L420: */
	    }
	}
	if (nlp > 0) {
	    i__2 = nlp;
	    for (j = 1; j <= i__2; ++j) {
		ii = (integer) (i__ + lp[j]);
		if (ii > 0) {
		    c1 = x[ii] - amu;
		}
		if (ii <= 0) {
		    c1 = xx[*nt + ii] - amu;
		}
/* L440: */
		c__ -= alpha[j] * c1;
	    }
	}
/* L400: */
	xx[*nt + i__] = c__ + amu;
    }

/*   now get e(-nt+1),...,e(n) and put them in indices 1,...,n+nt of e: */


L715:
    ss = (float)0.;
    ij = 1 - *nt;
    i__1 = *n;
    for (i__ = -(*nt) + 1; i__ <= i__1; ++i__) {
	if (i__ > 0) {
	    c__ = x[i__] - amu;
	}
	if (i__ <= 0) {
	    c__ = xx[*nt + i__] - amu;
	}
	if (nlp > 0) {
	    i__2 = nlp;
	    for (j = 1; j <= i__2; ++j) {
		ii = (integer) (i__ - lp[j]);
		c1 = (float)0.;
		if (ii > 0) {
		    c1 = x[ii] - amu;
		}
		if (ii <= 0 && ii >= ij) {
		    c1 = xx[*nt + ii] - amu;
		}
/* L20: */
		c__ += alpha[j] * c1;
	    }
	}
	if (nlq > 0) {
	    i__2 = nlq;
	    for (k = 1; k <= i__2; ++k) {
		ii = (integer) (i__ - lq[k]);
		if (ii >= ij) {
		    c__ -= beta[k] * e[*nt + ii];
		}
/* L40: */
	    }
	}
	e[*nt + i__] = c__;
/* L60: */
	ss += c__ * c__;
    }


    ret_val = ss;
    return ret_val;
} /* seaslk_ */


/* ******************************************************************* */

/* Subroutine */ int convt_(nords, coeffs, lags, maxp, maxq, nlp, nlq, lp, lq,
	 alpha, beta, a, nl, amu)
integer *nords;
doublereal *coeffs;
integer *lags, *maxp, *maxq, *nlp, *nlq;
doublereal *lp, *lq, *alpha, *beta, *a;
integer *nl;
doublereal *amu;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal c1;
    static integer n1, n2, n3, n4, jj, np, nq, npl, nql;


/*   Subroutine to convert multiplicative stuff to subset arma. */

/*   Input: */

/*      nords  : p,P,q,Q,M */
/*      coeffs : full AR, subset AR, full MA, subset MA coefficients */
/*               and (if M=1) mean */
/*      lags   : subset AR and MA lags */

/*   Output: */

/*      maxp, maxq  : maximum AR and MA lags in subset form */
/*      nlp, nlq    : number of AR and MA lags in subset form */
/*      lp, lq      : arrays of length nlp, nlq having AR and MA lags */
/*      alpha, beta : arrays of length nlp, nlq having AR and MA coeffs */
/*      amu         : mean */

/*   Auxilliary: */

/*      a, nl : arrays of length max(maxp,maxq) */


/* ******************************************************************* */



    /* Parameter adjustments */
    --nl;
    --a;
    --beta;
    --alpha;
    --lq;
    --lp;
    --lags;
    --coeffs;
    --nords;

    /* Function Body */
    np = nords[1];
    npl = nords[2];
    nq = nords[3];
    nql = nords[4];
    if (nords[5] == 0) {
	*amu = 0.;
    }
    if (nords[5] == 1) {
	*amu = coeffs[np + npl + nq + nql + 1];
    }
    n1 = 0;
    n2 = n1 + np;
    n3 = n2 + npl;
    n4 = n3 + nq;


    *maxp = np;
    *maxq = nq;
    if (npl > 0) {
	*maxp += lags[npl];
    }
    if (nql > 0) {
	*maxq += lags[npl + nql];
    }


    *nlp = 0;
    *nlq = 0;


    if (*maxp == 0) {
	goto L50;
    }


    i__1 = *maxp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nl[i__] = 0;
/* L10: */
	a[i__] = (float)0.;
    }


    if (np > 0) {
	i__1 = np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    nl[i__] = 1;
/* L15: */
	    a[i__] = coeffs[i__];
	}
    }


    if (npl == 0) {
	goto L30;
    }

    i__1 = np;
    for (i__ = 0; i__ <= i__1; ++i__) {
	c1 = 1.;
	if (i__ > 0) {
	    c1 = coeffs[i__];
	}
	i__2 = npl;
	for (j = 1; j <= i__2; ++j) {
	    jj = i__ + lags[j];
	    nl[jj] = 1;
/* L20: */
	    a[jj] += c1 * coeffs[n2 + j];
	}
    }

L30:

    i__2 = *maxp;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (nl[i__] == 0) {
	    goto L40;
	}
	++(*nlp);
	lp[*nlp] = (doublereal) i__;
	alpha[*nlp] = a[i__];
L40:
	;
    }


L50:


    if (*maxq == 0) {
	goto L100;
    }

    i__2 = *maxq;
    for (i__ = 1; i__ <= i__2; ++i__) {
	nl[i__] = 0;
/* L60: */
	a[i__] = (float)0.;
    }

    if (nq > 0) {
	i__2 = nq;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    nl[i__] = 1;
/* L65: */
	    a[i__] = coeffs[n3 + i__];
	}
    }

    if (nql == 0) {
	goto L80;
    }

    i__2 = nq;
    for (i__ = 0; i__ <= i__2; ++i__) {
	c1 = 1.;
	if (i__ > 0) {
	    c1 = coeffs[n3 + i__];
	}
	i__1 = nql;
	for (j = 1; j <= i__1; ++j) {
	    jj = i__ + lags[npl + j];
	    nl[jj] = 1;
/* L70: */
	    a[jj] += c1 * coeffs[n4 + j];
	}
    }

L80:


    i__1 = *maxq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nl[i__] == 0) {
	    goto L90;
	}
	++(*nlq);
	lq[*nlq] = (doublereal) i__;
	beta[*nlq] = a[i__];
L90:
	;
    }


L100:
    return 0;
} /* convt_ */




/* ******************************************************************* */

/* Subroutine */ int sspr_(x, n, nords, coeffs, lags, conf, ntf, ntl, nhl, e, 
	xx, xp, xpl, xpu, ier, npds, rvar)
doublereal *x;
integer *n, *nords;
doublereal *coeffs;
integer *lags;
real *conf;
integer *ntf, *ntl, *nhl;
doublereal *e, *xx, *xp, *xpl, *xpu;
integer *ier, *npds;
doublereal *rvar;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double pow_di(), log(), pow_dd(), sqrt(), exp();

    /* Local variables */
    static doublereal alam, beta[100];
    static integer imin, ijks, nthh, maxp, maxq;
    extern /* Subroutine */ int vmin_();
    static doublereal xmin;
    static integer ntph, ntot;
    static doublereal alam1;
    static integer maxp2, maxq2;
    static doublereal a[100], c__;
    static integer i__, j, k;
    static doublereal alpha[100];
    extern doublereal fctr10_();
    static doublereal a1[100], c1, c2, c3;
    extern /* Subroutine */ int convt_();
    static doublereal am;
    static integer nd, ii, jj, nh, nl[100];
    static doublereal lp[100], lq[100];
    static integer ns, nt, iptlam;
    extern /* Subroutine */ int movxyr_();
    static integer ndd;
    static doublereal amu;
    static integer nlp, nlq, ntt;



/* ******************************************************************* */



    /* Parameter adjustments */
    --xpu;
    --xpl;
    --xp;
    --xx;
    --e;
    --coeffs;
    --nords;
    --x;

    /* Function Body */
    *ier = 0;

/*   express model (without differencing) as subset ARMA: */

    convt_(&nords[1], &coeffs[1], lags, &maxp, &maxq, &nlp, &nlq, lp, lq, 
	    alpha, beta, a, nl, &amu);

/*   adjust for differences: */

    nd = nords[6];
    ndd = nords[7];
    ns = nords[8];
    maxq2 = maxq;
    maxp2 = maxp + nd + ndd * ns;

    if (maxp2 == 0) {
	goto L26;
    }
    i__1 = maxp2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ - 1] = 0.;
/* L10: */
	nl[i__ - 1] = 0;
    }
    i__1 = nlp;
    for (i__ = 0; i__ <= i__1; ++i__) {
	c1 = 1.;
	ii = 0;
	if (i__ > 0) {
	    ii = (integer) lp[i__ - 1];
	    c1 = alpha[i__ - 1];
	}
	i__2 = nd;
	for (j = 0; j <= i__2; ++j) {
	    c2 = pow_di(&c_b375, &j) * fctr10_(&nd, &j);
	    i__3 = ndd;
	    for (k = 0; k <= i__3; ++k) {
		c3 = pow_di(&c_b375, &k) * fctr10_(&ndd, &k);
		ijks = ii + j + k * ns;
		if (ijks > 0) {
		    nl[ijks - 1] = 1;
		    a[ijks - 1] += c1 * c2 * c3;
		}
/* L20: */
	    }
	}
    }
    nlp = 0;
    i__3 = maxp2;
    for (i__ = 1; i__ <= i__3; ++i__) {
	if (nl[i__ - 1] == 1) {
	    ++nlp;
	    lp[nlp - 1] = (doublereal) i__;
	    alpha[nlp - 1] = a[i__ - 1];
	}
/* L25: */
    }
L26:

/*   now nlp,nlq,alpha,beta,lp,lq are lags and coefficients of */
/*   full model */

/*   transform the series: */

    ntot = nords[1] + nords[2] + nords[3] + nords[4] + nords[5];
    am = coeffs[ntot + 1];
    alam = coeffs[ntot + 2];

/*   add the constant to make data positive: */

    i__3 = *n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L200: */
	x[i__] += am;
    }

/*   power transform: */

    iptlam = 1;
    if (alam == (float)1.) {
	goto L231;
    }

/*   check data + constant all positive: */

    vmin_(&x[1], n, &xmin, &imin);
    if (xmin <= (float)0.) {
	*ier = 1;
	goto L99;
    }

/*   do log transform: */

    if (alam == (float)0.) {
	iptlam = 2;
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L210: */
	    x[i__] = log(x[i__]);
	}
    }

/*   other power transforms: */

    if (abs(alam) > (float)1e-4 && abs(alam) < (float)1.) {
	iptlam = 3;
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L220: */
	    x[i__] = pow_dd(&x[i__], &alam);
	}
    }

L231:

/*   forecasts: first one possible is x(maxp2+maxq2+1): */


/*   find one step ahead forecast errors for x(maxp2+1),...,x(n): */

    i__3 = *n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L250: */
	e[i__] = (float)0.;
    }
    i__3 = *n;
    for (i__ = maxp2 + 1; i__ <= i__3; ++i__) {
	c__ = x[i__] - amu;
	if (nlp > 0) {
	    i__2 = nlp;
	    for (j = 1; j <= i__2; ++j) {
		if (i__ - lp[j - 1] >= 1.) {
		    c__ += alpha[j - 1] * x[(integer) (i__ - lp[j - 1])];
		}
/* L270: */
	    }
	}
	if (nlq > 0) {
	    i__2 = nlq;
	    for (j = 1; j <= i__2; ++j) {
		if (i__ - lq[j - 1] >= 1.) {
		    c__ -= beta[j - 1] * e[(integer) (i__ - lq[j - 1])];
		}
/* L280: */
	    }
	}
/* L260: */
	e[i__] = c__;
    }

/*   find forecasts: */

    *npds = 1;

    i__3 = *ntl;
    for (nt = *ntf; nt <= i__3; ++nt) {
	movxyr_(&xx[1], &x[1], &nt);

	i__2 = *nhl;
	for (nh = 1; nh <= i__2; ++nh) {
	    ntph = nt + nh;
	    c__ = amu;
	    if (maxq2 > 0) {
		i__1 = nlq;
		for (j = 1; j <= i__1; ++j) {
		    jj = (integer) (ntph - lq[j - 1]);
		    if (jj < 1 || jj > nt) {
			goto L370;
		    }
		    c__ += beta[j - 1] * e[(integer) (ntph - lq[j - 1])];
L370:
		    ;
		}
	    }
	    if (maxp2 > 0) {
		i__1 = nlp;
		for (j = 1; j <= i__1; ++j) {
/* L380: */
		    c__ -= alpha[j - 1] * xx[(integer) (ntph - lp[j - 1])];
		}
	    }
/* L360: */
	    xx[ntph] = c__;
	}

	movxyr_(&xp[*npds], &xx[nt + 1], nhl);
	*npds += *nhl;
/* L350: */
    }
    --(*npds);

/*   find psi weights (put them in xx): */

    if (maxp2 > 0) {
	i__3 = maxp2;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L410: */
	    a[i__ - 1] = (float)0.;
	}
	i__3 = nlp;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L420: */
	    a[(integer) lp[i__ - 1] - 1] = alpha[i__ - 1];
	}
    }
    if (maxq2 > 0) {
	i__3 = maxq2;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L430: */
	    a1[i__ - 1] = (float)0.;
	}
	i__3 = nlq;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L440: */
	    a1[(integer) lq[i__ - 1] - 1] = beta[i__ - 1];
	}
    }

    i__3 = *nhl;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L450: */
	xx[i__] = (float)0.;
    }
    if (maxq2 > 0) {
	i__3 = maxq2;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L451: */
	    xx[i__] = a1[i__ - 1];
	}
    }
    if (maxp2 > 0) {
	i__3 = *nhl;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    c__ = xx[i__];
	    i__2 = min(i__,maxp2);
	    for (j = 1; j <= i__2; ++j) {
		c1 = (float)1.;
		if (i__ - j > 0) {
		    c1 = xx[i__ - j];
		}
/* L452: */
		c__ -= a[j - 1] * c1;
	    }
/* L453: */
	    xx[i__] = c__;
	}
    }

/*        do 450 i=1,nhl */
/*           c=0. */
/*           if(i.le.maxq2) c=a1(i) */
/*           if(maxp2.gt.0) then */
/*              do 460 j=1,min0(maxp2,i) */
/*                 c1=1. */
/*                 if(i-j.ne.0) c1=xx(i-j) */
/* 460          c=c-a(j)*c1 */
/*           endif */
/* 450    xx(i)=c */

    e[1] = *rvar;
    if (*nhl > 1) {
	i__3 = *nhl;
	for (i__ = 2; i__ <= i__3; ++i__) {
	    c__ = (float)1.;
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
/* L480: */
		c__ += xx[j] * xx[j];
	    }
/* L470: */
	    e[i__] = *rvar * c__;
	}
    }
    i__3 = *nhl;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L471: */
	xx[i__] = sqrt(e[i__]);
    }

    i__3 = *nhl;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L490: */
	e[i__] = *conf * sqrt(e[i__]);
    }

    i__3 = *ntl;
    for (nt = *ntf; nt <= i__3; ++nt) {
	ntt = (nt - *ntf) * *nhl;
	i__2 = *nhl;
	for (nh = 1; nh <= i__2; ++nh) {
	    nthh = ntt + nh;
	    xpl[nthh] = xp[nthh] - e[nh];
/* L500: */
	    xpu[nthh] = xp[nthh] + e[nh];
	}
    }

/*   transform back: */

    if (iptlam == 3) {
	vmin_(&xpl[1], npds, &xmin, &imin);
	if (xmin <= (float)0.) {
	    *ier = 1;
	    goto L99;
	}
	alam1 = (float)1. / alam;
	i__2 = *npds;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xp[i__] = pow_dd(&xp[i__], &alam1);
	    xpl[i__] = pow_dd(&xpl[i__], &alam1);
/* L510: */
	    xpu[i__] = pow_dd(&xpu[i__], &alam1);
	}
    }

    if (iptlam == 2) {
	i__2 = *npds;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xp[i__] = exp(xp[i__]);
	    xpl[i__] = exp(xpl[i__]);
/* L520: */
	    xpu[i__] = exp(xpu[i__]);
	}
    }

/* L600: */
    i__2 = *npds;
    for (i__ = 1; i__ <= i__2; ++i__) {
	xp[i__] -= am;
	xpl[i__] -= am;
/* L610: */
	xpu[i__] -= am;
    }


L99:
    return 0;
} /* sspr_ */


/***************************************************************************/

doublereal fctr10_(n, k)
integer *n, *k;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static doublereal c__;
    static integer i__;
    static doublereal c1, c2;
    static integer kk;


/*   function to find the binomial coefficient { n choose k }. */

/************************************************************************
**/


    if (*k == 0 || *k == *n) {
	ret_val = 1.;
	return ret_val;
    }
/* Computing MIN */
    i__1 = *k, i__2 = *n - *k;
    kk = min(i__1,i__2);
    c__ = 1.;
    i__1 = kk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c1 = (doublereal) (*n - i__ + 1);
	c2 = (doublereal) i__;
/* L10: */
	c__ = c__ * c1 / c2;
    }
    ret_val = c__;
    return ret_val;
} /* fctr10_ */


/* ********************************************************* */

/* Subroutine */ int vmin_(x, n, xmin, ind)
doublereal *x;
integer *n;
doublereal *xmin;
integer *ind;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* ********************************************************* */


    /* Parameter adjustments */
    --x;

    /* Function Body */
    *ind = 1;
    *xmin = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] < *xmin) {
	    *xmin = x[i__];
	    *ind = i__;
	}
/* L10: */
    }
    return 0;
} /* vmin_ */


/* ********************************************************************** */

/* Subroutine */ int swpk12_(a, ndim, n, k1, k2, ier)
doublereal *a;
integer *ndim, *n, *k1, *k2, *ier;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal diag;
    static integer i__, j, k;


/*   Subroutine to sweep the (n x n) matrix a on diagonals k1 through k2. 
*/

/*   The output integer ier is 1 if an error occurs and 0 otherwise. */

/************************************************************************
***/



    /* Parameter adjustments */
    a_dim1 = *ndim;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *ier = 1;

    i__1 = *k2;
    for (k = *k1; k <= i__1; ++k) {
	diag = a[k + k * a_dim1];
	if (abs(diag) < (float)1e-20) {
	    return 0;
	}

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/* L10: */
		if (i__ != k && j != k) {
		    a[i__ + j * a_dim1] -= a[i__ + k * a_dim1] * a[k + j * 
			    a_dim1] / diag;
		}
	    }
	}

	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] / diag;
/* L20: */
	    a[k + i__ * a_dim1] /= diag;
	}

	a[k + k * a_dim1] = 1 / diag;

/* L30: */
    }

    *ier = 0;


    return 0;
} /* swpk12_ */

/* ******************************************************************** */

/* Subroutine */ int wilson_(alpha, np, beta, nq, acf, ma, cvli, mxpqp1, alph,
	 mxpq, ier)
doublereal *alpha;
integer *np;
doublereal *beta;
integer *nq;
doublereal *acf;
integer *ma;
doublereal *cvli;
integer *mxpqp1;
doublereal *alph;
integer *mxpq, *ier;
{
    /* Initialized data */

    static real epsil2 = (float)1e-20;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer mikp, i__, j, k, j1, miim1p, kc;
    static doublereal div;



/*   subroutine to find mx covariances */

/* ******************************************************************** */


    /* Parameter adjustments */
    --alph;
    --cvli;
    --acf;
    --beta;
    --alpha;

    /* Function Body */


    *ier = 0;
    if (*np < 0 || *nq < 0) {
	*ier = 1;
    }
    if (*mxpq != max(*np,*nq)) {
	*ier = 2;
    }
    if (*mxpqp1 != *mxpq + 1) {
	*ier = 3;
    }
    if (*ma < *mxpqp1) {
	*ier = 4;
    }
    if (*ier > 0) {
	return 0;
    }


    acf[1] = 1.;
    cvli[1] = 1.;
    if (*ma == 1) {
	return 0;
    }
    i__1 = *ma;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L1: */
	acf[i__] = 0.;
    }
    if (*mxpqp1 == 1) {
	return 0;
    }
    i__1 = *mxpqp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L2: */
	cvli[i__] = 0.;
    }
    i__1 = *mxpq;
    for (k = 1; k <= i__1; ++k) {
/* L9: */
	alph[k] = 0.;
    }


    if (*nq == 0) {
	goto L18;
    }
    i__1 = *nq;
    for (k = 1; k <= i__1; ++k) {
	cvli[k + 1] = -beta[k];
	acf[k + 1] = -beta[k];
	kc = *nq - k;
	if (kc == 0) {
	    goto L12;
	}
	i__2 = kc;
	for (j = 1; j <= i__2; ++j) {
/* L11: */
	    acf[k + 1] += beta[j] * beta[j + k];
	}
L12:
	acf[1] += beta[k] * beta[k];
/* L13: */
    }


L18:
    if (*np == 0) {
	return 0;
    }
    i__1 = *np;
    for (k = 1; k <= i__1; ++k) {
	alph[k] = alpha[k];
/* L19: */
	cvli[k] = alpha[k];
    }


    i__1 = *mxpq;
    for (k = 1; k <= i__1; ++k) {
	kc = *mxpq - k;
	if (kc >= *np) {
	    goto L24;
	}
	div = (float)1. - alph[kc + 1] * alph[kc + 1];
	if (div <= epsil2) {
	    goto L70;
	}
	if (kc == 0) {
	    goto L29;
	}
	i__2 = kc;
	for (j = 1; j <= i__2; ++j) {
/* L23: */
	    alph[j] = (cvli[j] + alph[kc + 1] * cvli[kc + 1 - j]) / div;
	}
L24:
	if (kc >= *nq) {
	    goto L26;
	}
/* Computing MAX */
	i__2 = kc + 1 - *np;
	j1 = max(i__2,1);
	i__2 = kc;
	for (j = j1; j <= i__2; ++j) {
/* L25: */
	    acf[j + 1] += acf[kc + 2] * alph[kc + 1 - j];
	}
L26:
	if (kc >= *np) {
	    goto L29;
	}
	i__2 = kc;
	for (j = 1; j <= i__2; ++j) {
/* L27: */
	    cvli[j] = alph[j];
	}
L29:
	;
    }


    acf[1] *= (float).5;
    i__1 = *mxpq;
    for (k = 1; k <= i__1; ++k) {
	if (k > *np) {
	    goto L33;
	}
	div = (float)1. - alph[k] * alph[k];
	i__2 = k + 1;
	for (j = 1; j <= i__2; ++j) {
/* L31: */
	    cvli[j] = (acf[j] + alph[k] * acf[k + 2 - j]) / div;
	}
	i__2 = k + 1;
	for (j = 1; j <= i__2; ++j) {
/* L32: */
	    acf[j] = cvli[j];
	}
L33:
	;
    }


    i__1 = *ma;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = i__ - 1;
	miim1p = min(i__2,*np);
	if (miim1p == 0) {
	    goto L43;
	}
	i__2 = miim1p;
	for (j = 1; j <= i__2; ++j) {
/* L42: */
	    acf[i__] += alpha[j] * acf[i__ - j];
	}
L43:
	;
    }
    acf[1] *= (float)2.;


    cvli[1] = (float)1.;
    if (*nq <= 0) {
	goto L60;
    }
    i__1 = *nq;
    for (k = 1; k <= i__1; ++k) {
	cvli[k + 1] = -beta[k];
	if (*np == 0) {
	    goto L53;
	}
	mikp = min(k,*np);
	i__2 = mikp;
	for (j = 1; j <= i__2; ++j) {
/* L52: */
	    cvli[k + 1] += alpha[j] * cvli[k + 1 - j];
	}
L53:
	;
    }

L60:
    return 0;


L70:
    *ier = 5;
    return 0;
} /* wilson_ */

