/* fanny.f -- translated by f2c (version 19950110).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>
#include "f2c.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define pow_dd(a,b) (pow(*a,*b))

/* Subroutine */ int fanny_(
integer *nn, integer *jpp, integer *kk,
doublereal *x, doublereal *dss,
integer *jdyss,
doublereal *valmd,
integer *jtmd, integer *ndyst, integer *nsend, integer *nelem, integer *negbr,
doublereal *syl, doublereal *p, doublereal *dp, doublereal *pt,
integer *nfuzz,
doublereal *esp, doublereal *ef, doublereal *dvec, doublereal *ttsyl, doublereal *eda, doublereal *edb, doublereal *obj,
integer *ncluv,
doublereal *sylinf, doublereal *eps)
{
    /* System generated locals */
    integer x_dim1, x_offset, p_dim1, p_offset, dp_dim1, dp_offset, 
	    sylinf_dim1, sylinf_offset;

    /* Local variables */
    static integer l;
    static doublereal s;
    extern /* Subroutine */ int caddy_();
    static integer nhalf, jhalt, ktrue;
    extern /* Subroutine */ int fygur_(), fuzzy_(), dysta3_();

/* C */
/* C   PROGRAM FOR FUZZY CLUSTER ANALYSIS */
/* C */
/* C  dimension of NSEND,NEGBR,NELEM,NCLUV,DVEC,SYL is MAXNN: */
/* C  dim. X(MAXNN,MAXPP),P(MAXNN,MAXKK),DP(MAXNN,MAXKK),DSS(MAXHH): */
/* C  dim. VALMD,JTMD,ESP,EF,PT,NFUZZ(MAXKK): */
/* C */
/* C   WHERE: */
/* C         NN = NUMBER OF OBJECTS */
/* C         JPP = NUMBER OF VARIABLES FOR CLUSTERING */
/* C         KK = NUMBER OF CLUSTERS */
/* C         MAXHH = (MAXNN*(MAXNN-1))/2 + 1 */
/* C */
    /* Parameter adjustments */
    sylinf_dim1 = *nn;
    sylinf_offset = sylinf_dim1 + 1;
    sylinf -= sylinf_offset;
    --ncluv;
    --dvec;
    --syl;
    --negbr;
    --nelem;
    --nsend;
    --dss;
    --jtmd;
    --valmd;
    x_dim1 = *nn;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --ef;
    --esp;
    --nfuzz;
    --pt;
    dp_dim1 = *nn;
    dp_offset = dp_dim1 + 1;
    dp -= dp_offset;
    p_dim1 = *nn;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --obj;

    /* Function Body */
    if (*jdyss == 1) {
	goto L125;
    }
    jhalt = 0;
    dysta3_(nn, jpp, &x[x_offset], &dss[1], ndyst, &jtmd[1], &valmd[1], &
	    jhalt);
    if (jhalt == 0) {
	goto L125;
    }
    *jdyss = -1;
    return 0;
/* C */
L125:
    s = (float)0.;
    nhalf = *nn * (*nn - 1) / 2 + 1;
    l = 1;
L130:
    ++l;
    if (dss[l] > s) {
	s = dss[l];
    }
    if (l < nhalf) {
	goto L130;
    }
    fuzzy_(nn, &nhalf, &p[p_offset], &dp[dp_offset], &pt[1], &dss[1], &esp[1],
	     &ef[1], eda, edb, kk, &obj[1], eps);
    caddy_(nn, &p[p_offset], kk, &ktrue, &nfuzz[1], &ncluv[1], &pt[1], &nelem[
	    1]);
    if (ktrue <= 1) {
	goto L140;
    }
    if (ktrue >= *nn) {
	goto L140;
    }
    fygur_(&ktrue, nn, kk, &nhalf, &ncluv[1], &nsend[1], &nelem[1], &negbr[1],
	     &syl[1], &dvec[1], &pt[1], ttsyl, &dss[1], &s, &sylinf[
	    sylinf_offset]);
L140:
    ;
} /* fanny_ */

/* C */
/* C */
/* C */
/* Subroutine */ int dysta3_(
integer *nn, integer *jpp,
doublereal *x, doublereal *dss,
integer *ndyst, integer *jtmd,
doublereal *valmd,
integer *jhalt)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer j, k, l, nnsub, npres;
    static doublereal rpres;
    static integer lplus;
    static doublereal pp, clk;
    static integer nlk;

    /* Parameter adjustments */
    --dss;
    --valmd;
    --jtmd;
    x_dim1 = *nn;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    pp = (doublereal) (*jpp);
    nnsub = *nn - 1;
    nlk = 0;
    i__1 = nnsub;
    for (l = 1; l <= i__1; ++l) {
	lplus = l + 1;
	i__2 = *nn;
	for (k = lplus; k <= i__2; ++k) {
	    clk = (float)0.;
	    ++nlk;
	    npres = 0;
	    i__3 = *jpp;
	    for (j = 1; j <= i__3; ++j) {
		if (jtmd[j] >= 0) {
		    goto L40;
		}
		if (x[l + j * x_dim1] == valmd[j]) {
		    goto L30;
		}
		if (x[k + j * x_dim1] == valmd[j]) {
		    goto L30;
		}
L40:
		++npres;
		if (*ndyst != 1) {
		    goto L50;
		}
		clk += (x[l + j * x_dim1] - x[k + j * x_dim1]) * (x[l + j * 
			x_dim1] - x[k + j * x_dim1]);
		goto L30;
L50:
		clk += (d__1 = x[l + j * x_dim1] - x[k + j * x_dim1], abs(
			d__1));
L30:
		;
	    }
	    rpres = (doublereal) npres;
	    if (npres != 0) {
		goto L60;
	    }
	    *jhalt = 1;
	    dss[nlk] = (float)-1.;
	    goto L20;
L60:
	    if (*ndyst != 1) {
		goto L70;
	    }
	    dss[nlk] = sqrt(clk * (pp / rpres));
	    goto L20;
L70:
	    dss[nlk] = clk * (pp / rpres);
L20:
	    ;
	}
/* L100: */
    }
} /* dysta3_ */

/* C */
/* C */
/* Subroutine */ int fuzzy_(
integer *nn, integer *hh,
doublereal *p, doublereal *dp, doublereal *pt, doublereal *dss, doublereal *esp, doublereal *ef, doublereal *eda, doublereal *edb,
integer *k,
doublereal *obj, doublereal *eps)
{
    /* System generated locals */
    integer p_dim1, p_offset, dp_dim1, dp_offset, i__1, i__2, i__3;
    doublereal d__1;


    /* Local variables */
    static doublereal reen, rkme, cryt;
    static integer j, l, m;
    static doublereal r;
    static integer kaunt, nnsub, j1, j2;
    static doublereal rvers;
    static integer nd;
    static doublereal dt;
    static integer lx;
    static doublereal zk, xx, ddd;
    static integer ndk;
    static doublereal ann, crt;
    static integer nyt;

/* C */
/* C     R IS THE EXPONENT, STRICTLY LARGER THAN 1.0 */
/* C     EPS IS THE PRECISION FOR THE ITERATIONS */
/* C     NYT IS THE MAXIMAL NUMBER OF ITERATIONS */
/* C */
    /* Parameter adjustments */
    --dss;
    --ef;
    --esp;
    --pt;
    dp_dim1 = *nn;
    dp_offset = dp_dim1 + 1;
    dp -= dp_offset;
    p_dim1 = *nn;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --obj;

    /* Function Body */
    r = (float)2.;
    nyt = 500;
/* C */
/* C   INITIAL FUZZY CLUSTERING */
/* C */
    nnsub = *nn - 1;
    rvers = (float)1. / r;
    rkme = (doublereal) (*k - 1);
    i__1 = *nn;
    for (m = 1; m <= i__1; ++m) {
	i__2 = *k;
	for (l = 1; l <= i__2; ++l) {
	    dp[m + l * dp_dim1] = (float)0.;
	    p[m + l * p_dim1] = (float).1 / rkme;
/* L20: */
	}
/* L30: */
    }
    ndk = *nn / *k;
    nd = ndk;
    l = 1;
    i__1 = *nn;
    for (m = 1; m <= i__1; ++m) {
	p[m + l * p_dim1] = (float).9;
	if (m < nd) {
	    goto L35;
	}
	nd += ndk;
	++l;
	if (l == *k) {
	    nd = *nn;
	}
L35:
	i__2 = *k;
	for (lx = 1; lx <= i__2; ++lx) {
	    p[m + lx * p_dim1] = pow_dd(&p[m + lx * p_dim1], &r);
/* L40: */
	}
/* L50: */
    }
/* C */
/* C   INITIAL CRITERION VALUE */
/* C */
    cryt = (float)0.;
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	esp[l] = (float)0.;
	ef[l] = (float)0.;
	i__2 = *nn;
	for (m = 1; m <= i__2; ++m) {
	    esp[l] += p[m + l * p_dim1];
	    i__3 = *nn;
	    for (j = 1; j <= i__3; ++j) {
		if (j == m) {
		    goto L80;
		}
		j2 = min(m,j);
		j1 = (j2 - 1) * *nn - j2 * (j2 + 1) / 2 + max(m,j);
		dp[m + l * dp_dim1] += p[j + l * p_dim1] * dss[j1];
		ef[l] += p[j + l * p_dim1] * p[m + l * p_dim1] * dss[j1];
L80:
		;
	    }
/* L90: */
	}
	cryt += ef[l] / (esp[l] * (float)2.);
/* L100: */
    }
    crt = cryt;
    reen = (float)1. / (r - (float)1.);
/* C */
/* C   START OF ITERATIONS */
/* C */
    kaunt = 1;
    m = 0;
/* C */
/* C   THE NEW MEMBERSHIP COEFFICIENTS OF THE OBJECTS ARE CALCULATED, */
/* C   AND THE RESULTING VALUE OF THE CRITERION IS COMPUTED. */
/* C */
L200:
    ++m;
    dt = (float)0.;
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	d__1 = esp[l] * (float)2. * esp[l] / (esp[l] * (float)2. * dp[m + l * 
		dp_dim1] - ef[l]);
	pt[l] = pow_dd(&d__1, &reen);
	dt += pt[l];
/* L210: */
    }
    xx = (float)0.;
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	pt[l] /= dt;
	if (pt[l] <= (float)0.) {
	    xx += pt[l];
	}
/* L220: */
    }
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	if (pt[l] <= (float)0.) {
	    pt[l] = (float)0.;
	}
	d__1 = pt[l] / (1 - xx);
	pt[l] = pow_dd(&d__1, &r);
	esp[l] = esp[l] + pt[l] - p[m + l * p_dim1];
	i__2 = *nn;
	for (j = 1; j <= i__2; ++j) {
	    if (j == m) {
		goto L230;
	    }
	    j2 = min(m,j);
	    j1 = (j2 - 1) * *nn - j2 * (j2 + 1) / 2 + max(m,j);
	    ddd = (pt[l] - p[m + l * p_dim1]) * dss[j1];
	    dp[j + l * dp_dim1] += ddd;
	    ef[l] += p[j + l * p_dim1] * (float)2. * ddd;
L230:
	    ;
	}
	p[m + l * p_dim1] = pt[l];
/* L240: */
    }
    if (m < *nn) {
	goto L200;
    }
    cryt = (float)0.;
    *eda = (float)0.;
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	ann = (doublereal) (*nn);
	*eda += esp[l] / ann;
	cryt += ef[l] / (esp[l] * (float)2.);
/* L250: */
    }
/* C */
/* C   CRITERION IS PRINTED AND TESTED FOR CONVERGENCE */
/* C */
    if (crt / cryt - (float)1. <= *eps) {
	goto L500;
    }
    if (kaunt < nyt) {
	goto L300;
    }
    goto L500;
L300:
    m = 0;
    ++kaunt;
    crt = cryt;
    goto L200;
/* C */
/* C   NON-FUZZYNESS INDEX OF LIBERT IS COMPUTED */
/* C */
L500:
    obj[1] = (doublereal) kaunt;
    obj[2] = cryt;
    zk = (doublereal) (*k);
    *edb = (zk * *eda - (float)1.) / (zk - (float)1.);
    i__1 = *nn;
    for (m = 1; m <= i__1; ++m) {
	i__2 = *k;
	for (l = 1; l <= i__2; ++l) {
	    p[m + l * p_dim1] = pow_dd(&p[m + l * p_dim1], &rvers);
/* L510: */
	}
/* L520: */
    }
    return 0;
} /* fuzzy_ */

/* C */
/* C */
/* Subroutine */ int caddy_(
integer *nn,
doublereal *p,
integer *k, integer *ktrue, integer *nfuzz, integer *ncluv,
doublereal *rdraw,
integer *nelem)
{
    /* System generated locals */
    integer p_dim1, p_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer ksup, ktry, l, m, kleft, kwalk, nbest;
    static doublereal pbest;
    static integer knext, jstay, lfuzz;

    /* Parameter adjustments */
    --nelem;
    --ncluv;
    --rdraw;
    --nfuzz;
    p_dim1 = *nn;
    p_offset = p_dim1 + 1;
    p -= p_offset;

    /* Function Body */
    pbest = p[p_dim1 + 1];
    nbest = 1;
    i__1 = *k;
    for (l = 2; l <= i__1; ++l) {
	if (p[l * p_dim1 + 1] <= pbest) {
	    goto L10;
	}
	pbest = p[l * p_dim1 + 1];
	nbest = l;
L10:
	;
    }
    nfuzz[1] = nbest;
    ncluv[1] = 1;
    *ktrue = 1;
    i__1 = *nn;
    for (m = 2; m <= i__1; ++m) {
	pbest = p[m + p_dim1];
	nbest = 1;
	i__2 = *k;
	for (l = 2; l <= i__2; ++l) {
	    if (p[m + l * p_dim1] <= pbest) {
		goto L30;
	    }
	    pbest = p[m + l * p_dim1];
	    nbest = l;
L30:
	    ;
	}
	jstay = 0;
	i__2 = *ktrue;
	for (ktry = 1; ktry <= i__2; ++ktry) {
	    if (nfuzz[ktry] != nbest) {
		goto L40;
	    }
	    ncluv[m] = ktry;
	    jstay = 1;
L40:
	    ;
	}
	if (jstay == 1) {
	    goto L20;
	}
	++(*ktrue);
	nfuzz[*ktrue] = nbest;
	ncluv[m] = *ktrue;
L20:
	;
    }
    if (*ktrue >= *k) {
	goto L100;
    }
    knext = *ktrue + 1;
    i__1 = *k;
    for (kwalk = knext; kwalk <= i__1; ++kwalk) {
	i__2 = *k;
	for (kleft = 1; kleft <= i__2; ++kleft) {
	    jstay = 0;
	    ksup = kwalk - 1;
	    i__3 = ksup;
	    for (ktry = 1; ktry <= i__3; ++ktry) {
		if (nfuzz[ktry] != kleft) {
		    goto L80;
		}
		jstay = 1;
L80:
		;
	    }
	    if (jstay == 1) {
		goto L70;
	    }
	    nfuzz[kwalk] = kleft;
	    goto L60;
L70:
	    ;
	}
L60:
	;
    }
L100:
    i__1 = *nn;
    for (m = 1; m <= i__1; ++m) {
	i__2 = *k;
	for (l = 1; l <= i__2; ++l) {
	    lfuzz = nfuzz[l];
	    rdraw[l] = p[m + lfuzz * p_dim1];
/* L120: */
	}
	i__2 = *k;
	for (l = 1; l <= i__2; ++l) {
	    p[m + l * p_dim1] = rdraw[l];
/* L130: */
	}
/* L110: */
    }
} /* caddy_ */

/* C */
/* C */
/* Subroutine */ int fygur_(integer *ktrue, integer *nn, integer *kk, integer *hh, integer *ncluv, integer *nsend, integer *nelem, integer *negbr, doublereal *syl, doublereal *srank, doublereal *avsyl, doublereal *ttsyl, doublereal *dss, doublereal *s, doublereal *sylinf)
{
    /* System generated locals */
    integer sylinf_dim1, sylinf_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer lang;
    static doublereal dysa, dysb;
    static integer nclu, j, l, lplac, numcl;
    static doublereal symax;
    static integer nsylr;
    static doublereal db;
    static integer nj, nl, nbb, mjl, njl;
    static doublereal att, btt, rnn;
    static integer ntt;
    static doublereal rtt;

    /* Parameter adjustments */
    sylinf_dim1 = *nn;
    sylinf_offset = sylinf_dim1 + 1;
    sylinf -= sylinf_offset;
    --srank;
    --syl;
    --negbr;
    --nelem;
    --nsend;
    --ncluv;
    --avsyl;
    --dss;

    /* Function Body */
    nsylr = 0;
    *ttsyl = (float)0.;
    i__1 = *ktrue;
    for (numcl = 1; numcl <= i__1; ++numcl) {
	ntt = 0;
	i__2 = *nn;
	for (j = 1; j <= i__2; ++j) {
	    if (ncluv[j] != numcl) {
		goto L30;
	    }
	    ++ntt;
	    nelem[ntt] = j;
L30:
	    ;
	}
	i__2 = ntt;
	for (j = 1; j <= i__2; ++j) {
	    nj = nelem[j];
	    dysb = *s * (float)1.1 + (float)1.;
	    negbr[j] = -1;
	    i__3 = *ktrue;
	    for (nclu = 1; nclu <= i__3; ++nclu) {
		if (nclu == numcl) {
		    goto L41;
		}
		nbb = 0;
		db = (float)0.;
		i__4 = *nn;
		for (l = 1; l <= i__4; ++l) {
		    if (ncluv[l] != nclu) {
			goto L43;
		    }
		    ++nbb;
		    if (l < nj) {
			goto L42;
		    }
		    if (l > nj) {
			goto L44;
		    }
		    goto L43;
L42:
		    mjl = *nn * (l - 1) + nj - l * (l + 1) / 2;
		    db += dss[mjl];
		    goto L43;
L44:
		    mjl = *nn * (nj - 1) + l - nj * (nj + 1) / 2;
		    db += dss[mjl];
L43:
		    ;
		}
		btt = (doublereal) nbb;
		db /= btt;
		if (db >= dysb) {
		    goto L41;
		}
		dysb = db;
		negbr[j] = nclu;
L41:
		;
	    }
	    if (ntt == 1) {
		goto L50;
	    }
	    dysa = (float)0.;
	    i__3 = ntt;
	    for (l = 1; l <= i__3; ++l) {
		nl = nelem[l];
		if (nj < nl) {
		    goto L46;
		}
		if (nj > nl) {
		    goto L47;
		}
		goto L45;
L46:
		njl = *nn * (nj - 1) + nl - nj * (nj + 1) / 2;
		dysa += dss[njl];
		goto L45;
L47:
		njl = *nn * (nl - 1) + nj - nl * (nl + 1) / 2;
		dysa += dss[njl];
L45:
		;
	    }
	    att = (doublereal) (ntt - 1);
	    dysa /= att;
	    if (dysa > (float)0.) {
		goto L51;
	    }
	    if (dysb > (float)0.) {
		goto L52;
	    }
L50:
	    syl[j] = (float)0.;
	    goto L40;
L52:
	    syl[j] = (float)1.;
	    goto L40;
L51:
	    if (dysb <= (float)0.) {
		goto L53;
	    }
	    if (dysb > dysa) {
		syl[j] = (float)1. - dysa / dysb;
	    }
	    if (dysb < dysa) {
		syl[j] = dysb / dysa - (float)1.;
	    }
	    if (dysb == dysa) {
		syl[j] = (float)0.;
	    }
	    goto L54;
L53:
	    syl[j] = (float)-1.;
L54:
	    if (syl[j] <= (float)-1.) {
		syl[j] = (float)-1.;
	    }
	    if (syl[j] >= (float)1.) {
		syl[j] = (float)1.;
	    }
L40:
	    ;
	}
	avsyl[numcl] = (float)0.;
	i__2 = ntt;
	for (j = 1; j <= i__2; ++j) {
	    symax = (float)-2.;
	    i__3 = ntt;
	    for (l = 1; l <= i__3; ++l) {
		if (syl[l] <= symax) {
		    goto L70;
		}
		symax = syl[l];
		lang = l;
L70:
		;
	    }
	    nsend[j] = lang;
	    srank[j] = syl[lang];
	    avsyl[numcl] += srank[j];
	    syl[lang] = (float)-3.;
/* L60: */
	}
	*ttsyl += avsyl[numcl];
	rtt = (doublereal) ntt;
	avsyl[numcl] /= rtt;
	if (ntt >= 2) {
	    goto L75;
	}
	++nsylr;
	sylinf[nsylr + sylinf_dim1] = (doublereal) numcl;
	sylinf[nsylr + (sylinf_dim1 << 1)] = (doublereal) negbr[1];
	sylinf[nsylr + sylinf_dim1 * 3] = (float)0.;
	sylinf[nsylr + (sylinf_dim1 << 2)] = (doublereal) nelem[1];
	goto L100;
L75:
	i__2 = ntt;
	for (l = 1; l <= i__2; ++l) {
	    ++nsylr;
	    lplac = nsend[l];
	    sylinf[nsylr + sylinf_dim1] = (doublereal) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (doublereal) negbr[lplac];
	    sylinf[nsylr + sylinf_dim1 * 3] = srank[l];
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (doublereal) nelem[lplac];
/* L80: */
	}
L100:
	;
    }
    rnn = (doublereal) (*nn);
    *ttsyl /= rnn;
} /* fygur_ */

