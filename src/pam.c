/* pam.f -- translated by f2c (version 19950110).
You must link the resulting object file with the libraries:
-lf2c -lm   (in that order)
*/
#ifdef __cplusplus
extern "C" {
#endif
	
	
#include <math.h>
#include "f2c.h"
	
	/* Subroutine */ int dark_(
	integer *kk, integer *nn, integer *hh,integer *ncluv,integer *nsend,integer *nelem,integer *negbr,
		doublereal *syl,doublereal  *srank,doublereal  *avsyl,doublereal  *ttsyl,doublereal  *dys,doublereal  *s,doublereal  *sylinf);
	/* Subroutine */ int cstat_(
	integer *kk, integer *nn, integer *nsend, integer *nrepr,
		doublereal *radus, doublereal *damer, doublereal *ttd, doublereal *separ, doublereal *z, doublereal *s,
		integer *hh,
		doublereal *dys,
		integer *ncluv, integer *nelem, integer *med, integer *nisol);
	/* Subroutine */ int bswap_(
	integer *kk, integer *nn, integer *nrepr,
		doublereal *dysma, doublereal *dysmb, doublereal *beter,
		integer *hh,
		doublereal *dys, doublereal *sky, doublereal *s, doublereal *obj);
	/* Subroutine */ int dysta_(
	integer *nn, integer *jpp,
		doublereal *x, doublereal  *dys,
		integer *ndyst, integer *jtmd,
		doublereal *valmd,
		integer *jhalt);
	
	
	/* Subroutine */ int pam_(
	integer *nn, integer *jpp, integer *kk,
		doublereal *x, doublereal *dys,
		integer *jdyss,
		doublereal *valmd,
		integer *jtmd, integer *ndyst, integer *nsend, integer *nrepr, integer *nelem,
		doublereal *radus, doublereal *damer, doublereal *ttd, doublereal *separ, doublereal *ttsyl,
		integer *med,
		doublereal *obj,
		integer *ncluv,
		doublereal *clusinf, doublereal *sylinf,
		integer *nisol)
	{
		/* System generated locals */
		integer x_dim1, x_offset, clusinf_dim1, clusinf_offset, sylinf_dim1, 
			sylinf_offset, i__1;
		
		/* Local variables */
		extern /* Subroutine */ int dark_();
		static integer k, l;
		static doublereal s;
		static integer nhalf, jhalt;
		extern /* Subroutine */ int bswap_(), cstat_(), dysta_();
		static doublereal sky;
		
		/* C */
		/* C    PARTITIONING AROUND MEDOIDS */
		/* C */
		/* C    CARRIES OUT A CLUSTERING USING THE K-MEDOID APPROACH. */
		/* C */
		/* Parameter adjustments */
		sylinf_dim1 = *nn;
		sylinf_offset = sylinf_dim1 + 1;
		sylinf -= sylinf_offset;
		--ncluv;
		--separ;
		--ttd;
		--damer;
		--radus;
		--nelem;
		--nrepr;
		--nsend;
		--dys;
		--jtmd;
		--valmd;
		x_dim1 = *nn;
		x_offset = x_dim1 + 1;
		x -= x_offset;
		--nisol;
		clusinf_dim1 = *kk;
		clusinf_offset = clusinf_dim1 + 1;
		clusinf -= clusinf_offset;
		--med;
		--obj;
		
		/* Function Body */
		if (*jdyss == 1) {
			goto L125;
		}
		jhalt = 0;
		dysta_(nn, jpp, &x[x_offset], &dys[1], ndyst, &jtmd[1], &valmd[1], &jhalt)
			;
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
		if (dys[l] > s) {
			s = dys[l];
		}
		if (l < nhalf) {
			goto L130;
		}
		bswap_(kk, nn, &nrepr[1], &radus[1], &damer[1], &ttd[1], &nhalf, &dys[1], 
			&sky, &s, &obj[1]);
		cstat_(kk, nn, &nsend[1], &nrepr[1], &radus[1], &damer[1], &ttd[1], &
			separ[1], &sky, &s, &nhalf, &dys[1], &ncluv[1], &nelem[1], &med[1]
			, &nisol[1]);
		i__1 = *kk;
		for (k = 1; k <= i__1; ++k) {
			clusinf[k + clusinf_dim1] = (doublereal) nrepr[k];
			clusinf[k + (clusinf_dim1 << 1)] = radus[k];
			clusinf[k + clusinf_dim1 * 3] = ttd[k];
			clusinf[k + (clusinf_dim1 << 2)] = damer[k];
			clusinf[k + clusinf_dim1 * 5] = separ[k];
			/* L135: */
		}
		if (*kk <= 1) {
			goto L140;
		}
		if (*kk >= *nn) {
			goto L140;
		}
		dark_(kk, nn, &nhalf, &ncluv[1], &nsend[1], &nelem[1], &nrepr[1], &radus[
			1], &damer[1], &ttd[1], ttsyl, &dys[1], &s, &sylinf[sylinf_offset]
			);
L140:
		;
		} /* pam_ */
		
		/* C */
		/* C */
		/* Subroutine */ int dysta_(
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
		} /* dysta_ */
		
		/* C */
		/* C */
		/* Subroutine */ int bswap_(
		integer *kk, integer *nn, integer *nrepr,
			doublereal *dysma, doublereal *dysmb, doublereal *beter,
			integer *hh,
			doublereal *dys, doublereal *sky, doublereal *s, doublereal *obj)
		{
			/* System generated locals */
			integer i__1, i__2, i__3;
			
			/* Local variables */
			static integer njaj;
			extern integer meet_();
			static integer nmax, j, k;
			static doublereal ammax, small;
			static integer kbest, nbest;
			static doublereal dzsky;
			static integer ja;
			static doublereal dz, cmd;
			static integer nkj, njn;
			static doublereal rnn;
			static integer nny;
			
			/* C */
			/* C    FIRST ALGORITHM: BUILD. */
			/* C */
			/* Parameter adjustments */
			--beter;
			--dysmb;
			--dysma;
			--nrepr;
			--dys;
			--obj;
			
			/* Function Body */
			nny = 0;
			i__1 = *nn;
			for (j = 1; j <= i__1; ++j) {
				nrepr[j] = 0;
				dysma[j] = *s * (float)1.1 + (float)1.;
				/* L17: */
			}
L20:
			i__1 = *nn;
			for (ja = 1; ja <= i__1; ++ja) {
				if (nrepr[ja] != 0) {
					goto L22;
				}
				beter[ja] = (float)0.;
				i__2 = *nn;
				for (j = 1; j <= i__2; ++j) {
					njaj = meet_(&ja, &j);
					cmd = dysma[j] - dys[njaj];
					if (cmd > (float)0.) {
						beter[ja] += cmd;
					}
					/* L21: */
				}
L22:
				;
			}
			ammax = (float)0.;
			i__1 = *nn;
			for (ja = 1; ja <= i__1; ++ja) {
				if (nrepr[ja] != 0) {
					goto L31;
				}
				if (beter[ja] < ammax) {
					goto L31;
				}
				ammax = beter[ja];
				nmax = ja;
L31:
				;
			}
			nrepr[nmax] = 1;
			++nny;
			i__1 = *nn;
			for (j = 1; j <= i__1; ++j) {
				njn = meet_(&nmax, &j);
				if (dys[njn] < dysma[j]) {
					dysma[j] = dys[njn];
				}
				/* L41: */
			}
			if (nny != *kk) {
				goto L20;
			}
			*sky = (float)0.;
			i__1 = *nn;
			for (j = 1; j <= i__1; ++j) {
				*sky += dysma[j];
				/* L51: */
			}
			rnn = (doublereal) (*nn);
			obj[1] = *sky / rnn;
			if (*kk == 1) {
				goto L75;
			}
			/* C */
			/* C    SECOND ALGORITHM: SWAP. */
			/* C */
L60:
			i__1 = *nn;
			for (j = 1; j <= i__1; ++j) {
				dysma[j] = *s * (float)1.1 + (float)1.;
				dysmb[j] = *s * (float)1.1 + (float)1.;
				i__2 = *nn;
				for (ja = 1; ja <= i__2; ++ja) {
					if (nrepr[ja] == 0) {
						goto L62;
					}
					njaj = meet_(&ja, &j);
					if (dys[njaj] >= dysma[j]) {
						goto L61;
					}
					dysmb[j] = dysma[j];
					dysma[j] = dys[njaj];
					goto L62;
L61:
					if (dys[njaj] >= dysmb[j]) {
						goto L62;
					}
					dysmb[j] = dys[njaj];
L62:
					;
				}
				/* L63: */
			}
			dzsky = (float)1.;
			i__1 = *nn;
			for (k = 1; k <= i__1; ++k) {
				if (nrepr[k] == 1) {
					goto L73;
				}
				i__2 = *nn;
				for (ja = 1; ja <= i__2; ++ja) {
					if (nrepr[ja] == 0) {
						goto L72;
					}
					dz = (float)0.;
					i__3 = *nn;
					for (j = 1; j <= i__3; ++j) {
						njaj = meet_(&ja, &j);
						nkj = meet_(&k, &j);
						if (dys[njaj] != dysma[j]) {
							goto L70;
						}
						small = dysmb[j];
						if (dys[nkj] < small) {
							small = dys[nkj];
						}
						dz = dz - dysma[j] + small;
						goto L71;
L70:
						if (dys[nkj] < dysma[j]) {
							dz = dz - dysma[j] + dys[nkj];
						}
L71:
						;
					}
					if (dz >= dzsky) {
						goto L72;
					}
					dzsky = dz;
					kbest = k;
					nbest = ja;
L72:
					;
				}
L73:
				;
			}
			if (dzsky >= -(float)1e-10) {
				goto L75;
			}
			nrepr[kbest] = 1;
			nrepr[nbest] = 0;
			*sky += dzsky;
			goto L60;
L75:
			rnn = (doublereal) (*nn);
			obj[2] = *sky / rnn;
} /* bswap_ */

/* C */
/* C */
/* Subroutine */ int cstat_(
integer *kk, integer *nn, integer *nsend, integer *nrepr,
doublereal *radus, doublereal *damer, doublereal *ttd, doublereal *separ, doublereal *z, doublereal *s,
integer *hh,
doublereal *dys,
integer *ncluv, integer *nelem, integer *med, integer *nisol)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
	
    /* Local variables */
    static integer kand, njaj;
    extern integer meet_();
    static integer mevj, nvna, jndz, numl, j, k, m, nplac;
    static doublereal dsmal;
    static integer ksmal, numcl, ja, jb, jk;
    static doublereal aja, ajb, dam;
    static integer nel, njm;
    static doublereal sep, rnn;
    static integer nvn, ntt;
    static doublereal rtt, ttt;
	
    /* Parameter adjustments */
    --nisol;
    --med;
    --nelem;
    --ncluv;
    --separ;
    --ttd;
    --damer;
    --radus;
    --nrepr;
    --nsend;
    --dys;
	
    /* Function Body */
    i__1 = *nn;
    for (j = 1; j <= i__1; ++j) {
		if (nrepr[j] == 1) {
			goto L120;
		}
		dsmal = *s * (float)1.1 + (float)1.;
		i__2 = *nn;
		for (k = 1; k <= i__2; ++k) {
			if (nrepr[k] == 0) {
				goto L110;
			}
			njaj = meet_(&k, &j);
			if (dys[njaj] >= dsmal) {
				goto L110;
			}
			dsmal = dys[njaj];
			ksmal = k;
L110:
			;
		}
		nsend[j] = ksmal;
		goto L130;
L120:
		nsend[j] = j;
L130:
		;
    }
    jk = 1;
    nplac = nsend[1];
    i__1 = *nn;
    for (j = 1; j <= i__1; ++j) {
		ncluv[j] = 0;
		if (nsend[j] == nplac) {
			ncluv[j] = 1;
		}
		/* L135: */
    }
    i__1 = *nn;
    for (ja = 2; ja <= i__1; ++ja) {
		nplac = nsend[ja];
		if (ncluv[nplac] != 0) {
			goto L145;
		}
		++jk;
		i__2 = *nn;
		for (j = 2; j <= i__2; ++j) {
			if (nsend[j] == nplac) {
				ncluv[j] = jk;
			}
			/* L140: */
		}
		if (jk == *kk) {
			goto L148;
		}
L145:
		;
    }
	/* C */
	/* C    ANALYSIS OF THE CLUSTERING. */
	/* C */
L148:
    i__1 = *kk;
    for (numcl = 1; numcl <= i__1; ++numcl) {
		ntt = 0;
		radus[numcl] = (float)-1.;
		ttt = (float)0.;
		i__2 = *nn;
		for (j = 1; j <= i__2; ++j) {
			if (ncluv[j] != numcl) {
				goto L150;
			}
			++ntt;
			m = nsend[j];
			nelem[ntt] = j;
			njm = meet_(&j, &m);
			ttt += dys[njm];
			if (dys[njm] > radus[numcl]) {
				radus[numcl] = dys[njm];
			}
L150:
			;
		}
		rtt = (doublereal) ntt;
		ttd[numcl] = ttt / rtt;
		med[numcl] = m;
		/* L160: */
    }
	/* L230: */
    rnn = (doublereal) (*nn);
    if (*kk != 1) {
		goto L240;
    }
    damer[1] = *s;
    nrepr[1] = *nn;
    goto L300;
	/* C */
	/* C    NUML = NUMBER OF L-CLUSTERS. */
	/* C */
L240:
    numl = 0;
    i__1 = *kk;
    for (k = 1; k <= i__1; ++k) {
		/* C */
		/* C    IDENTIFICATION OF CLUSTER K: */
		/* C       NEL=NUMBER OF OBJECTS */
		/* C       NELEM=VECTOR OF OBJECTS */
		/* C */
		nel = 0;
		i__2 = *nn;
		for (j = 1; j <= i__2; ++j) {
			if (ncluv[j] != k) {
				goto L23;
			}
			++nel;
			nelem[nel] = j;
L23:
			;
		}
		nrepr[k] = nel;
		if (nel != 1) {
			goto L24;
		}
		nvn = nelem[1];
		damer[k] = (float)0.;
		separ[k] = *s * (float)1.1 + (float)1.;
		i__2 = *nn;
		for (j = 1; j <= i__2; ++j) {
			if (j == nvn) {
				goto L250;
			}
			mevj = meet_(&nvn, &j);
			if (separ[k] > dys[mevj]) {
				separ[k] = dys[mevj];
			}
L250:
			;
		}
		/* C */
		/* C    IS CLUSTER K     1) AN L-CLUSTER ? */
		/* C                     2) AN L*-CLUSTER ? */
		/* C */
		if (separ[k] == (float)0.) {
			goto L400;
		}
		++numl;
L400:
		goto L35;
L24:
		dam = (float)-1.;
		sep = *s * (float)1.1 + (float)1.;
		kand = 1;
		i__2 = nel;
		for (ja = 1; ja <= i__2; ++ja) {
			nvna = nelem[ja];
			aja = (float)-1.;
			ajb = *s * (float)1.1 + (float)1.;
			i__3 = *nn;
			for (jb = 1; jb <= i__3; ++jb) {
				jndz = meet_(&nvna, &jb);
				if (ncluv[jb] == k) {
					goto L30;
				}
				if (dys[jndz] < ajb) {
					ajb = dys[jndz];
				}
				goto L25;
L30:
				if (dys[jndz] > aja) {
					aja = dys[jndz];
				}
L25:
				;
			}
			if (aja >= ajb) {
				kand = 0;
			}
			if (dam < aja) {
				dam = aja;
			}
			if (sep > ajb) {
				sep = ajb;
			}
			/* L26: */
		}
		separ[k] = sep;
		damer[k] = dam;
		if (kand == 0) {
			goto L35;
		}
		++numl;
		if (dam < sep) {
			goto L27;
		}
		/* C L-CLUSTER */
		nisol[k] = 1;
		goto L40;
		/* C L*-CLUSTER */
L27:
		nisol[k] = 2;
		goto L40;
L35:
		nisol[k] = 0;
L40:
		;
    }
L300:
    ;
} /* cstat_ */

/* C */
/* C */
/* Subroutine */ int dark_(
integer *kk,  integer *nn,  integer *hh,  integer *ncluv,  integer *nsend,  integer *nelem,  integer *negbr,
doublereal *syl, doublereal *srank, doublereal *avsyl, doublereal *ttsyl, doublereal *dys, doublereal *s, doublereal *sylinf)
{
    /* System generated locals */
    integer sylinf_dim1, sylinf_offset, i__1, i__2, i__3, i__4;
	
    /* Local variables */
    static integer lang;
    extern integer meet_();
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
    --avsyl;
    --srank;
    --syl;
    --negbr;
    --nelem;
    --nsend;
    --ncluv;
    --dys;
	
    /* Function Body */
    nsylr = 0;
    *ttsyl = (float)0.;
    i__1 = *kk;
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
			i__3 = *kk;
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
					mjl = meet_(&nj, &l);
					db += dys[mjl];
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
				njl = meet_(&nj, &nl);
				dysa += dys[njl];
				/* L45: */
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
	/* L96: */
} /* dark_ */

