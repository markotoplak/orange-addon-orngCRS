/* twins.f -- translated by f2c (version 19950110).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
#ifdef __cplusplus
extern "C" {
#endif


#include <math.h>
#include "f2c.h"

int twins_(integer *nn,integer  *jpp, doublereal *x, doublereal *dys, doublereal *dys2,
integer *jdyss, doublereal *valmd, integer *jtmd, integer *ndyst, integer *jalg, 
integer *method, integer *kwan, integer *ner, doublereal *ban,doublereal  *coef,
integer *merge);
/* Subroutine */ int supcl_(
doublereal *dys,
integer *kka,integer  *kkb,
doublereal *arest,
integer *nn, integer *ner);


/* Subroutine */ int bandy_(
integer *nn,
doublereal *ban,
integer *ner,
doublereal *dc);

/* Subroutine */ int splyt_(
integer *nn, integer *kwan, integer *ner,
doublereal *ban, doublereal *dys,
integer *merge);


/* Subroutine */ int banag_(
integer *nn,
doublereal *ban,
integer *ner,
doublereal *ac);

/* Subroutine */ int averl_(
integer *nn, integer *kwan, integer *ner,
doublereal *ban,doublereal *dys,
integer *method, integer *merge);

int dysta4_(
	integer *nn, integer *jpp,
	doublereal *x, doublereal *dys,
	integer *ndyst, integer *jtmd,
	doublereal *valmd,
	integer *jhalt);

integer meet_(
integer *l, integer *j);



/* Table of constant values */

static integer c__1 = 1;

int twins_(integer *nn,integer  *jpp, doublereal *x, doublereal *dys, doublereal *dys2,
		   integer *jdyss, doublereal *valmd, integer *jtmd, integer *ndyst, integer *jalg, 
		   integer *method, integer *kwan, integer *ner, doublereal *ban,doublereal  *coef,
		   integer *merge) 
{
    /* System generated locals */
    integer x_dim1, x_offset, merge_dim1, merge_offset, i__1;
	
    /* Local variables */
    static integer i;
    static integer jhalt;
	
	/* C */
	/* C   THIS PROGRAM PERFORMS AGGLOMERATIVE NESTING (AGNES) USING THE */
	/* C   GROUP AVERAGE METHOD OF SOKAL AND MICHENER (1958), AS WELL AS */
	/* C   DIVISIVE ANALYSIS (DIANA) USING THE METHOD OF MCNAUGHTON-SMITH, */
	/* C   WILLIAMS, DALE, AND MOCKETT (1964). */
	/* C */
	/* C   LIST OF FUNCTIONS AND SUBROUTINES: */
	/* C       MAIN UNIT */
	/* C       FUNCTION MEET */
	/* C       SUBROUTINE DYSTA4 */
	/* C       SUBROUTINE AVERL */
	/* C       SUBROUTINE SPLYT */
	/* C       SUBROUTINE SUPCL */
	/* C */
	/* C   THE FOLLOWING VECTORS AND MATRICES MUST BE DIMENSIONED IN THE */
	/* C   MAIN PROGRAM ONLY: */
	/* C       KWAN(NN),NER(NN),BAN(NN) */
	/* C       X(NN,JPP),JTMD(JPP),VALMD(JPP),DYS((NN*(NN-1))/2 + 1) */
	/* C   WHERE: */
	/* C       NN = MAXIMAL NUMBER OF OBJECTS */
	/* C       JPP = MAXIMAL NUMBER OF VARIABLES USED IN THE ANALYSIS */
	/* C */
	/* C */
    /* Parameter adjustments */
    merge_dim1 = *nn - 1;
    merge_offset = merge_dim1 + 1;
    merge -= merge_offset;
    --ban;
    --ner;
    --kwan;
    --dys2;
    --dys;
    --jtmd;
    --valmd;
    x_dim1 = *nn;
    x_offset = x_dim1 + 1;
    x -= x_offset;
	
    /* Function Body */
    if (*jdyss == 0) {
		goto L70;
    }
    *jpp = 1;
    goto L100;
L70:
    jhalt = 0;
    dysta4_(nn, jpp, &x[x_offset], &dys[1], ndyst, &jtmd[1], &valmd[1], &
		jhalt);
    if (jhalt == 0) {
		goto L100;
    }
    *jdyss = -1;
    return 0;
L100:
    i__1 = *nn * (*nn - 1) / 2 + 1;
    for (i = 1; i <= i__1; ++i) {
		dys2[i] = dys[i];
		/* L110: */
    }
    if (*jalg == 2) {
		goto L200;
    }
    averl_(nn, &kwan[1], &ner[1], &ban[1], &dys[1], method, &merge[
		merge_offset]);
    banag_(nn, &ban[1], &ner[1], coef);
    goto L300;
L200:
    splyt_(nn, &kwan[1], &ner[1], &ban[1], &dys[1], &merge[merge_offset]);
    bandy_(nn, &ban[1], &ner[1], coef);
L300:
    ;
} /* twins_ */

/* C */
/* C */
/* Subroutine */ 
int dysta4_(
			integer *nn, integer *jpp,
			doublereal *x, doublereal *dys,
			integer *ndyst, integer *jtmd,
			doublereal *valmd,
			integer *jhalt)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer j, k, l, npres, lsubt;
    static doublereal rpres, pp, clk;
    static integer nlk;

    /* Parameter adjustments */
    --dys;
    --valmd;
    --jtmd;
    x_dim1 = *nn;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    pp = (doublereal) (*jpp);
    nlk = 1;
    dys[1] = (float)0.;
    i__1 = *nn;
    for (l = 2; l <= i__1; ++l) {
	lsubt = l - 1;
	i__2 = lsubt;
	for (k = 1; k <= i__2; ++k) {
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
	    dys[nlk] = (float)-1.;
	    goto L20;
L60:
	    if (*ndyst != 1) {
		goto L70;
	    }
	    dys[nlk] = sqrt(clk * (pp / rpres));
	    goto L20;
L70:
	    dys[nlk] = clk * (pp / rpres);
L20:
	    ;
	}
/* L100: */
    }
} /* dysta4_ */

/* C */
/* C */
/* Subroutine */ int averl_(
integer *nn, integer *kwan, integer *ner,
doublereal *ban,doublereal *dys,
integer *method, integer *merge)
{
    /* System generated locals */
    integer merge_dim1, merge_offset, i__1, i__2;
	
    /* Local variables */
    static doublereal dnew;
    static integer nclu, lnum, lput;
    static doublereal d;
    static integer j, l, lenda, lendb;
    static doublereal smald;
    static integer lmuch, llast, lnext, l1, l2, lfyrs;
    static doublereal fa, fb, fc;
    static integer la, lb;
    static doublereal ta, tb;
    static integer lq;
    static doublereal tq;
    static integer nmerge;
    static doublereal akb;
    static integer nab, lka, nej, naq, nbq, nlj, nns;
	
    /* Parameter adjustments */
    merge_dim1 = *nn - 1;
    merge_offset = merge_dim1 + 1;
    merge -= merge_offset;
    --dys;
    --ban;
    --ner;
    --kwan;
	
    /* Function Body */
    nclu = *nn - 1;
	/* C      INITIALIZATION */
    i__1 = *nn;
    for (l = 1; l <= i__1; ++l) {
		kwan[l] = 1;
		ner[l] = l;
		/* L10: */
    }
	/* C */
	/* C    FIND CLOSEST CLUSTERS */
	/* C */
    nmerge = 1;
L100:
    j = 1;
L80:
    ++j;
    if (kwan[j] == 0) {
		goto L80;
    }
    nej = meet_(&c__1, &j);
    smald = dys[nej] * (float)1.1 + (float)1.;
    nns = *nn - 1;
    i__1 = nns;
    for (l = 1; l <= i__1; ++l) {
		if (kwan[l] == 0) {
			goto L120;
		}
		lmuch = l + 1;
		i__2 = *nn;
		for (j = lmuch; j <= i__2; ++j) {
			if (kwan[j] == 0) {
				goto L110;
			}
			nlj = meet_(&l, &j);
			if (dys[nlj] > smald) {
				goto L110;
			}
			smald = dys[nlj];
			la = l;
			lb = j;
L110:
			;
		}
L120:
		;
    }
	/* C */
	/* C    MERGE-STRUCTURE FOR PLOTTING TREE IN S-PLUS */
	/* C */
    l1 = -la;
    l2 = -lb;
    if (nmerge == 1) {
		goto L121;
    }
    i__1 = nmerge - 1;
    for (j = 1; j <= i__1; ++j) {
		if (merge[j + merge_dim1] == l1 || merge[j + (merge_dim1 << 1)] == l1)
		{
			l1 = j;
		}
		if (merge[j + merge_dim1] == l2 || merge[j + (merge_dim1 << 1)] == l2)
		{
			l2 = j;
		}
		/* L122: */
    }
L121:
    merge[nmerge + merge_dim1] = l1;
    merge[nmerge + (merge_dim1 << 1)] = l2;
    ++nmerge;
	/* C */
	/* C    DETERMINE LFYRS AND LLAST */
	/* C */
    i__1 = *nn;
    for (l = 1; l <= i__1; ++l) {
		if (ner[l] == la) {
			lfyrs = l;
		}
		if (ner[l] == lb) {
			llast = l;
		}
		/* L200: */
    }
    ban[llast] = smald;
	/* C */
	/* C    IF THE TWO CLUSTERS ARE NEXT TO EACH OTHER, */
	/* C    NER MUST NOT BE CHANGED */
	/* C */
    lnext = lfyrs + kwan[la];
    if (lnext == llast) {
		goto L230;
    }
	/* C */
	/* C    UPDATING NER AND BAN */
	/* C */
    lput = lfyrs + kwan[la];
    lnum = llast - lput;
    i__1 = lnum;
    for (l = 1; l <= i__1; ++l) {
		lka = ner[lput];
		akb = ban[lput];
		lenda = llast + kwan[lb] - 2;
		lendb = lenda + 1;
		i__2 = lenda;
		for (j = lput; j <= i__2; ++j) {
			ner[j] = ner[j + 1];
			ban[j] = ban[j + 1];
			/* L210: */
		}
		ner[lendb] = lka;
		ban[lendb] = akb;
		/* L220: */
    }
	/* C */
	/* C    CALCULATE NEW DISSIMILARITIES */
	/* C */
L230:
    i__1 = *nn;
    for (lq = 1; lq <= i__1; ++lq) {
		if (lq == la || lq == lb) {
			goto L240;
		}
		if (kwan[lq] == 0) {
			goto L240;
		}
		naq = meet_(&la, &lq);
		nbq = meet_(&lb, &lq);
		if (*method == 2) {
			goto L300;
		}
		if (*method == 3) {
			goto L310;
		}
		if (*method == 4) {
			goto L320;
		}
		if (*method == 5) {
			goto L330;
		}
		/* C   GROUP AVERAGE METHOD */
		ta = (doublereal) kwan[la];
		tb = (doublereal) kwan[lb];
		fa = ta / (ta + tb);
		fb = tb / (ta + tb);
		dys[naq] = fa * dys[naq] + fb * dys[nbq];
		goto L240;
		/* C   SINGLE LINKAGE */
L300:
		dnew = dys[naq];
		if (dys[nbq] < dnew) {
			dnew = dys[nbq];
		}
		dys[naq] = dnew;
		goto L240;
		/* C   COMPLETE LINKAGE */
L310:
		dnew = dys[naq];
		if (dnew < dys[nbq]) {
			dnew = dys[nbq];
		}
		dys[naq] = dnew;
		goto L240;
		/* C   WARD'S METHOD */
L320:
		ta = (doublereal) kwan[la];
		tb = (doublereal) kwan[lb];
		tq = (doublereal) kwan[lq];
		fa = (ta + tq) / (ta + tb + tq);
		fb = (tb + tq) / (ta + tb + tq);
		fc = -tq / (ta + tb + tq);
		nab = meet_(&la, &lb);
		d = fa * dys[naq] * dys[naq] + fb * dys[nbq] * dys[nbq];
		d += fc * dys[nab] * dys[nab];
		dys[naq] = sqrt(d);
		goto L240;
		/* C   WEIGHTED AVERAGE LINKAGE */
L330:
		dys[naq] = (dys[naq] + dys[nbq]) / 2.;
L240:
		;
    }
	/* L250: */
    kwan[la] += kwan[lb];
    kwan[lb] = 0;
    --nclu;
    if (nclu > 0) {
		goto L100;
    }
} /* averl_ */

/* C */
/* C */
/* Subroutine */ int banag_(
integer *nn,
doublereal *ban,
integer *ner,
doublereal *ac)
{
    /* System generated locals */
    integer i__1;
	
    /* Local variables */
    static doublereal syze;
    static integer k, kafte, kearl;
    static doublereal rnn, sup;
	
    /* Parameter adjustments */
    --ner;
    --ban;
	
    /* Function Body */
    sup = (float)0.;
    i__1 = *nn;
    for (k = 2; k <= i__1; ++k) {
		if (ban[k] > sup) {
			sup = ban[k];
		}
		/* L70: */
    }
    *ac = (float)0.;
    i__1 = *nn;
    for (k = 1; k <= i__1; ++k) {
		kearl = k;
		if (k == 1) {
			kearl = 2;
		}
		kafte = k + 1;
		if (k == *nn) {
			kafte = *nn;
		}
		syze = ban[kearl];
		if (ban[kafte] < syze) {
			syze = ban[kafte];
		}
		*ac = *ac + (float)1. - syze / sup;
		/* L80: */
    }
    rnn = (doublereal) (*nn);
    *ac /= rnn;
	return 0;
} /* banag_ */

/* C */
/* C */
/* Subroutine */ int splyt_(
integer *nn, integer *kwan, integer *ner,
doublereal *ban, doublereal *dys,
integer *merge)
{
    /* System generated locals */
    integer merge_dim1, merge_offset, i__1, i__2;
	
    /* Local variables */
    static integer lmma, lmmb;
    static doublereal dyff;
    static integer lgrb;
    static doublereal dmin__;
    static integer jner, nclu, lner, lxxa;
    static doublereal rest;
    static integer lxxp, j, k, l, lchan, nhalf;
    static doublereal bdyff;
    static integer lndsd;
    static doublereal bygsd;
    static integer jaway;
    static doublereal arest;
    static integer l1, l2;
    static doublereal splyn, da, db;
    static integer ja, jb;
    static doublereal cs, sd;
    static integer nj, nmerge, jab, jma, jan, jbn, jmb, nlj, lmm, llq, lxf, 
		lxg, lmz, lxx, lxy;
	
	/* C */
	/* C    INITIALIZATION */
	/* C */
    /* Parameter adjustments */
    merge_dim1 = *nn - 1;
    merge_offset = merge_dim1 + 1;
    merge -= merge_offset;
    --dys;
    --ban;
    --ner;
    --kwan;
	
    /* Function Body */
    nclu = 1;
    nhalf = *nn * (*nn - 1) / 2 + 1;
    i__1 = *nn;
    for (l = 1; l <= i__1; ++l) {
		kwan[l] = 0;
		ban[l] = (float)0.;
		ner[l] = l;
		/* L10: */
    }
    kwan[1] = *nn;
    ja = 1;
	/* C */
	/* C    COMPUTATION OF DIAMETER OF DATA SET */
	/* C */
    cs = (float)0.;
    k = 0;
L20:
    ++k;
    if (dys[k] > cs) {
		cs = dys[k];
    }
    if (k < nhalf) {
		goto L20;
    }
	/* C */
	/* C    PREPARE FOR SPLITTING */
	/* C */
L30:
    jb = ja + kwan[ja] - 1;
    jma = jb;
	/* C */
	/* C    SPECIAL CASE OF A PAIR OF OBJECTS */
	/* C */
    if (kwan[ja] != 2) {
		goto L50;
    }
    kwan[ja] = 1;
    kwan[jb] = 1;
    jan = ner[ja];
    jbn = ner[jb];
    jab = meet_(&jan, &jbn);
    ban[jb] = dys[jab];
    goto L400;
	/* C */
	/* C    FINDING FIRST OBJECT TO BE SHIFTED */
	/* C */
L50:
    bygsd = (float)-1.;
    i__1 = jb;
    for (l = ja; l <= i__1; ++l) {
		lner = ner[l];
		sd = (float)0.;
		i__2 = jb;
		for (j = ja; j <= i__2; ++j) {
			jner = ner[j];
			nlj = meet_(&lner, &jner);
			sd += dys[nlj];
			/* L100: */
		}
		if (sd <= bygsd) {
			goto L110;
		}
		bygsd = sd;
		lndsd = l;
L110:
		;
    }
	/* C */
	/* C    SHIFTING THE FIRST OBJECT */
	/* C */
    --kwan[ja];
    kwan[jb] = 1;
    if (jb == lndsd) {
		goto L115;
    }
    lchan = ner[lndsd];
    lmm = jb - 1;
    i__1 = lmm;
    for (lmma = lndsd; lmma <= i__1; ++lmma) {
		lmmb = lmma + 1;
		ner[lmma] = ner[lmmb];
		/* L112: */
    }
    ner[jb] = lchan;
L115:
    splyn = (float)0.;
    jma = jb - 1;
	/* C */
	/* C    FINDING THE NEXT OBJECT TO BE SHIFTED */
	/* C */
L120:
    splyn += (float)1.;
    rest = (doublereal) (jma - ja);
    bdyff = (float)-1.;
    i__1 = jma;
    for (l = ja; l <= i__1; ++l) {
		lner = ner[l];
		da = (float)0.;
		i__2 = jma;
		for (j = ja; j <= i__2; ++j) {
			jner = ner[j];
			nlj = meet_(&lner, &jner);
			da += dys[nlj];
			/* L130: */
		}
		da /= rest;
		db = (float)0.;
		jmb = jma + 1;
		i__2 = jb;
		for (j = jmb; j <= i__2; ++j) {
			jner = ner[j];
			nlj = meet_(&lner, &jner);
			db += dys[nlj];
			/* L140: */
		}
		db /= splyn;
		dyff = da - db;
		if (dyff <= bdyff) {
			goto L150;
		}
		bdyff = dyff;
		jaway = l;
L150:
		;
    }
    jmb = jma + 1;
	/* C */
	/* C    SHIFTING THE NEXT OBJECT WHEN NECESSARY */
	/* C */
    if (bdyff <= (float)0.) {
		goto L200;
    }
    if (jma == jaway) {
		goto L165;
    }
    lchan = ner[jaway];
    lmz = jma - 1;
    i__1 = lmz;
    for (lxx = jaway; lxx <= i__1; ++lxx) {
		lxxp = lxx + 1;
		ner[lxx] = ner[lxxp];
		/* L160: */
    }
    ner[jma] = lchan;
L165:
    i__1 = jb;
    for (lxx = jmb; lxx <= i__1; ++lxx) {
		lxy = lxx - 1;
		if (ner[lxy] < ner[lxx]) {
			goto L180;
		}
		lchan = ner[lxy];
		ner[lxy] = ner[lxx];
		ner[lxx] = lchan;
		/* L170: */
    }
L180:
    --kwan[ja];
    kwan[jma] = kwan[jmb] + 1;
    kwan[jmb] = 0;
    --jma;
    jmb = jma + 1;
    if (jma != ja) {
		goto L120;
    }
	/* C */
	/* C    SWITCH THE TWO PARTS WHEN NECESSARY */
	/* C */
L200:
    if (ner[ja] < ner[jmb]) {
		goto L300;
    }
    lxxa = ja;
    i__1 = jb;
    for (lgrb = jmb; lgrb <= i__1; ++lgrb) {
		++lxxa;
		lchan = ner[lgrb];
		i__2 = lgrb;
		for (lxy = lxxa; lxy <= i__2; ++lxy) {
			lxf = lgrb - lxy + lxxa;
			lxg = lxf - 1;
			ner[lxf] = ner[lxg];
			/* L210: */
		}
		ner[lxg] = lchan;
		/* L220: */
    }
    llq = kwan[jmb];
    kwan[jmb] = 0;
    jma = ja + jb - jma - 1;
    jmb = jma + 1;
    kwan[jmb] = kwan[ja];
    kwan[ja] = llq;
	/* C */
	/* C    COMPUTE LEVEL FOR BANNER */
	/* C */
L300:
    if (nclu == 1) {
		ban[jmb] = cs;
    }
    if (nclu == 1) {
		goto L400;
    }
    supcl_(&dys[1], &ja, &jb, &arest, nn, &ner[1]);
    ban[jmb] = arest;
L400:
    ++nclu;
    if (nclu == *nn) {
		goto L500;
    }
	/* C */
	/* C    CONTINUE SPLITTING UNTIL ALL OBJECTS ARE SEPARATED */
	/* C */
    if (jb == *nn) {
		goto L430;
    }
L420:
    ja += kwan[ja];
    if (ja > *nn) {
		goto L430;
    }
    if (kwan[ja] <= 1) {
		goto L420;
    }
    goto L30;
L430:
    ja = 1;
    if (kwan[ja] == 1) {
		goto L420;
    }
    goto L30;
	/* C */
	/* C    MERGE-STRUCTURE FOR PLOTTING TREE IN S-PLUS */
	/* C */
L500:
    i__1 = *nn - 1;
    for (nmerge = 1; nmerge <= i__1; ++nmerge) {
		dmin__ = cs;
		i__2 = *nn;
		for (j = 2; j <= i__2; ++j) {
			if (kwan[j] < 0 || ban[j] > dmin__) {
				goto L560;
			}
			dmin__ = ban[j];
			nj = j;
L560:
			;
		}
		kwan[nj] = -1;
		l1 = -ner[nj - 1];
		l2 = -ner[nj];
		if (nmerge == 1) {
			goto L570;
		}
		i__2 = nmerge - 1;
		for (j = 1; j <= i__2; ++j) {
			if (merge[j + merge_dim1] == l1 || merge[j + (merge_dim1 << 1)] ==
				l1) {
				l1 = j;
			}
			if (merge[j + merge_dim1] == l2 || merge[j + (merge_dim1 << 1)] ==
				l2) {
				l2 = j;
			}
			/* L580: */
		}
L570:
		merge[nmerge + merge_dim1] = l1;
		merge[nmerge + (merge_dim1 << 1)] = l2;
		/* L550: */
    }
} /* splyt_ */

/* C */
/* C */
/* Subroutine */ int supcl_(
doublereal *dys,
integer *kka,integer  *kkb,
doublereal *arest,
integer *nn, integer *ner)
{
    /* System generated locals */
    integer i__1, i__2;
	
    /* Local variables */
    static integer jner, lner, j, l, kkc, kkd, mlj;
	
    /* Parameter adjustments */
    --ner;
    --dys;
	
    /* Function Body */
    kkc = *kkb - 1;
    *arest = (float)0.;
    i__1 = kkc;
    for (l = *kka; l <= i__1; ++l) {
		lner = ner[l];
		kkd = l + 1;
		i__2 = *kkb;
		for (j = kkd; j <= i__2; ++j) {
			jner = ner[j];
			mlj = meet_(&lner, &jner);
			if (dys[mlj] > *arest) {
				*arest = dys[mlj];
			}
			/* L10: */
		}
		/* L20: */
    }
    return 0;
} /* supcl_ */

/* C */
/* C */
/* Subroutine */ int bandy_(
integer *nn,
doublereal *ban,
integer *ner,
doublereal *dc)
{
    /* System generated locals */
    integer i__1;
	
    /* Local variables */
    static doublereal syze;
    static integer k, kafte, kearl;
    static doublereal rnn, sup;
	
    /* Parameter adjustments */
    --ner;
    --ban;
	
    /* Function Body */
    sup = (float)0.;
    i__1 = *nn;
    for (k = 2; k <= i__1; ++k) {
		if (ban[k] > sup) {
			sup = ban[k];
		}
		/* L70: */
    }
    *dc = (float)0.;
    i__1 = *nn;
    for (k = 1; k <= i__1; ++k) {
		kearl = k;
		if (k == 1) {
			kearl = 2;
		}
		kafte = k + 1;
		if (k == *nn) {
			kafte = *nn;
		}
		syze = ban[kearl];
		if (ban[kafte] < syze) {
			syze = ban[kafte];
		}
		*dc = *dc + (float)1. - syze / sup;
		/* L80: */
    }
    rnn = (doublereal) (*nn);
    *dc /= rnn;
	return 0;
} /* bandy_ */

integer meet_(
			  integer *l, integer *j)
{
    /* System generated locals */
    integer ret_val;
	
    if (*l > *j) {
		goto L10;
    }
    if (*l == *j) {
		goto L20;
    }
	/* C */
	/* C   L LESS THAN J */
	/* C */
    ret_val = (*j - 2) * (*j - 1) / 2 + *l + 1;
    return ret_val;
	/* C */
	/* C   J LESS THAN L */
	/* C */
L10:
    ret_val = (*l - 2) * (*l - 1) / 2 + *j + 1;
    return ret_val;
	/* C */
	/* C   J EQUALS L */
	/* C */
L20:
    ret_val = 1;
    return ret_val;
} /* meet_ */

#ifdef __cplusplus
}
#endif
