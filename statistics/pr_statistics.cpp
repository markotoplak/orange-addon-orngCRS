#include "pr_statistics.h"

// Uses a combination of 3 different algorithms for different N & D.
double ks2(int n, double d) {
	if (n <= 50)
		return 1.0 - pks2(n,d);
	else if (n > 500)
		return ks2_asympt(n,d);
	else 
		if (d <= 1.2 / sqrt(float(n)))
			if (n <= 100)
				return 1.0 - pks2(n,d);
			else
				return ks2_asympt(n,d);
		else
			return 2. * ks1(n,d);
}


//! Code converted using TO_F90 by Alan Miller
//! Date: 1999-12-28  Time: 18:34:23
//!     ALGORITHM APPEARED IN COMM. ACM, VOL. 17, NO. 12, P. 703.

double pks2(int n, double d) {
	/*
! N IS THE SAMPLE SIZE USED.

! D IS THE MAXIMUM MAGNITUDE (OF THE DISCREPANCY BETWEEN THE EMPIRICAL AND
! PROPOSED DISTRIBUTIONS) IN EITHER THE POSITIVE OR NEGATIVE DIRECTION.
! PKS2 IS THE EXACT PROBABILITY OF OBTAINING A DEVIATION NO LARGER THAN D.
! THESE FORMULAS APPEAR AS (23) AND (24) IN J. DURBIN.  THE PROBABILITY THAT
! THE SAMPLE DISTRIBUTION FUNCTION LIES BETWEEN TWO PARALLEL STRAIGHT LINES.
! ANNALS OF MATHEMATICAL STATISTICS 39, 2(APRIL 1968), 398-411.
	*/
	double *q;
	double *fact;
	double sum, ci, fn, fnd, ft, fu, fv, sign, prob;
	int i, j, jmax, k, nd, ndd, nddp, ndp, ndt;

	if (n == 1)
		return 2.*d - 1.0;
	fn = n;
	fnd = fn*d;
	ndt = int(2*fnd);
	if (ndt < 1)
		return 0.0;
	
	q = new double[n+2];
	fact = new double[n+2];
	nd = int(fnd);
	ndd = MIN(2*nd,n);
	ndp = nd + 1;
	nddp = ndd + 1;

	fact[1] = 1.0;
	ci = 1.0;
	for (i = 1; i <= n; ++i) {
		fact[i+1] = fact[i]*ci;
		ci += 1.0;
	}
	q[1] = 1.0;
	if (ndd != 0) {
		ci = 1.0;
		for(i=1; i <= ndd; ++i) {
			q[i+1] = pow(ci,i)/fact[i+1];
			ci += 1.0;
		}
		if (ndp > n) goto l80;

		fv = ndp - fnd;
		jmax = int(fv) + 1;
		for (i = ndp; i <= ndd; ++i) {
			sum = 0.0;
			ft = fnd;
			k = i;
			fu = fv;
			for (j = 1; j <= jmax; ++j) {
				sum += pow(ft,j-2)/fact[j]*pow(fu,k)/fact[k+1];
				ft += 1.0;
				fu -= 1.0;
				--k;
			}
			q[i+1] -= 2.*fnd*sum;
			++jmax;
			fv += 1.0;
		}
		if (ndd == n) goto l80;
	}
	for (i = nddp; i <= n; ++i) {
		sum = 0.0;
		sign = 1.0;
		ft = 2.*fnd;
		for (j = 1; j <= ndt; ++j) {
		    ft -= 1.0;
			k = i-j+1;
		    sum += sign*pow(ft,j)/fact[j+1]*q[k];
			sign = -sign;
		}
	    q[i+1] = sum;
	}

l80:
	prob = q[n+1]*fact[n+1] / pow(fn,n);
	delete[] q,fact;
	return prob;
}



double lngamma(double z) {
//  Uses Lanczos-type approximation to ln(gamma) for z > 0.
//  Reference:
//       Lanczos, C. 'A precision approximation of the gamma
//               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
//  Accuracy: About 14 significant digits except for small regions
//            in the vicinity of 1 and 2.

	double a[9] = {0.9999999999995183, 676.5203681218835,
                      -1259.139216722289, 771.3234287757674, 
                      -176.6150291498386, 12.50734324009056, 
                      -0.1385710331296526, 0.9934937113930748e-05,
					  0.1659470187408462e-06};
	double lnsqrt2pi = 0.9189385332046727;
	double lanczos,tmp;

	int j;

	if (z <= EPSILON)
		return 0;

	lanczos = 0;
	tmp = z + 7;
	for (j = 8; j >= 1; --j) {
	  lanczos = lanczos + a[j]/tmp;
	  tmp = tmp - 1;
	}
	lanczos = lanczos + a[0];
	lanczos = log(lanczos) + lnsqrt2pi - (z + 6.5) + (z - 0.5)*log(z + 6.5);

	return lanczos;
}


double log_nCr(int n, int r) {
//! Evaluate log of the number of combinations of r out of n.
	double fn_val;
	int i, j;

//! If MIN(r, n-r) < 10, use multiplications & divisions, otherwise
//! use the gamma function.

	j = MIN(r, n-r);
	if (j < 10) {
		fn_val = 1.0;
		for (i = 1; i <= j; ++i)
			fn_val *= (n+1-i) / double(i);
		return log(fn_val);
	} else {
		return lngamma(double(n+1)) - lngamma(double(r+1)) - lngamma(double(n+1-r));
	}
}


double ks1(int n, double d) {
//! Calculates the single-tailed probability for the on sample KS test.

	double eps, lncj, lncj0, term, prob;
	int j, j0, jupper;

	eps = 10.0 * EPSILON;
	jupper = int(n * (1.0 - d));

//! Evaluate one of the largest terms in the sum,
//! and then sum up & down from term j0 until terms are negligible

	j0 = MIN(jupper, n/2);
	lncj0 = log_nCr(n,j0);
	prob = exp( lncj0 + (n-j0)*log(1.0 - d - j0/double(n)) + (j0-1)*log(d + j0/double(n)) );

	lncj = lncj0;
	for (j = j0+1; j <= jupper; ++j) {
	  lncj += log(double(n+1-j) / double(j));
	  term = exp( lncj + (n-j)*log(1.0 - d - j/double(n)) + (j-1)*log(d + j/double(n)) );
	  prob += term;
	  if (term < eps*prob)
		  break;
	}
	
	lncj = lncj0;
	for (j = j0-1; j >= 0; --j) {
	  lncj += log(double(j+1) / double(n-j));
	  term = exp(lncj + (n-j)*log(1.0 - d - j/double(n)) + (j-1)*log(d + j/double(n)) );
	  prob += term;
	  if (term < eps*prob)
		  break;
	}

	prob *= d;
	return prob;
}

double ks2_asympt(int n, double d) {
// ! Using Smirnov's asymptotic formaula for large N

	double sgn, sum, term, z;
	int i;

	z = d * sqrt( double(n) );
	sum = 0.0;
	i = 1;
	sgn = 1.0;
	while(true) {
	  term = exp( -2.0 * pow(i*z,2) );
	  sum += sgn*term;
	  if (term < EPSILON)
		  break;
	  ++i;
	  sgn = -sgn;
	}
	return 2.0*sum;
}


//  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

// Evaluates the tail area of the standardised normal curve
// from x to infinity if upper is .true. or
// from minus infinity to x if upper is .false.
double alnorm(double x, bool upper) {
	double norm_prob, con=1.28, z,y, ltone = 7.0, utzero = 18.66;
	double p, r, a2, b1, c1, c3, c5, d1, d3, d5;
	double q, a1, a3, b2, c2, c4, c6, d2, d4;
	bool up;
	
	p = 0.398942280444;
	q = 0.39990348504;
	r = 0.398942280385; 
	a1 = 5.75885480458;
	a2 = 2.62433121679; 
	a3 = 5.92885724438; 
	b1 = -29.8213557807;
	b2 = 48.6959930692;
    c1 = -3.8052e-8;
	c2 = 3.98064794e-4;       
    c3 = -0.151679116635;
	c4 = 4.8385912808;
	c5 = 0.742380924027;
	c6 = 3.99019417011; 
    d1 = 1.00000615302; 
	d2 = 1.98615381364;  
	d3 = 5.29330324926;
	d4 = -15.1508972451;
    d5 = 30.789933034;

	up = upper;
	z = x;
	if(z < 0.0) {
		up = !up;
		z = -z;
	}
	if(z <= ltone || up && z <= utzero) {
		y = 0.5*z*z;
	} else {
		norm_prob = 0.0;
	}
	if(z > con) {
		norm_prob = r*exp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
	} else {
		norm_prob = 0.5 - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
	}
	if(!up) 
		norm_prob = 1.0 - norm_prob;

	return norm_prob;
}

double gammad(double x, double p) {
	//!  ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3
	//!  Computation of the Incomplete Gamma Integral
	//!  Auxiliary functions required: ALNORM = algorithm AS66 (included) & LNGAMMA
	//!  Converted to be compatible with ELF90 by Alan Miller
	//!  N.B. The return parameter IFAULT has been removed as ELF90 allows only
	//!  one output parameter from functions.   An error message is issued instead.
	
	double gamma_prob;
	double pn1, pn2, pn3, pn4, pn5, pn6, tol = 1.0e-14, oflo = 1.0e+37;
	double xbig = 1.0e8, arg, c, rn, a, b, one = 1.0, zero = 0.0, an;
	double two = 2.0, elimit = -88.0, plimit = 1000.0, three = 3.0;
    double nine = 9;
	
	gamma_prob = zero;
	
	if	(p <= zero || x < EPSILON) {
		return 0.0;
	}
	
	//      Use a normal approximation if P > PLIMIT
	if (p > plimit) {
		pn1 = three * sqrt(p) * (pow(x/p,one/three) + one / (nine * p) - one);
		return alnorm(pn1, false);
	}
	
	//      If X is extremely large compared to P then set gamma_prob = 1
	if (x > xbig) {
		return one;
	}
	
	if (x <= one || x < p) {
		//!      Use Pearson's series expansion.
		//!      (Note that P is not large enough to force overflow in LNGAMMA)
		
		arg = p * log(x) - x - lngamma(p + one);
		c = one;
		gamma_prob = one;
		a = p;
		do {
			a = a + one;
			c = c * x / a;
			gamma_prob = gamma_prob + c;
		} while (c >= tol);
		
		arg = arg + log(gamma_prob);
		gamma_prob = zero;
		if (arg >= elimit) {
			gamma_prob = exp(arg);
		} 
	} else {
		//!      Use a continued fraction expansion
		
		arg = p * log(x) - x - lngamma(p);
		a = one - p;
		b = a + x + one;
		c = zero;
		pn1 = one;
		pn2 = x;
		pn3 = x + one;
		pn4 = x * b;
		gamma_prob = pn3 / pn4;
		do {
			a = a + one;
			b = b + two;
			c = c + one;
			an = a * c;
			pn5 = b * pn3 - an * pn1;
			pn6 = b * pn4 - an * pn2;
			if (fabs(pn6) > zero) {
				rn = pn5 / pn6;
				if(fabs(gamma_prob - rn) <= MIN(tol, tol * rn))
					break;
				gamma_prob = rn;
			}
			
			pn1 = pn3;
			pn2 = pn4;
			pn3 = pn5;
			pn4 = pn6;
			if (fabs(pn5) >= oflo) {
				//  !      Re-scale terms in continued fraction if terms are large
				
				pn1 = pn1 / oflo;
				pn2 = pn2 / oflo;
				pn3 = pn3 / oflo;
				pn4 = pn4 / oflo;
			}
		} while (true);
		arg = arg + log(gamma_prob);
		gamma_prob = one;
		if (arg >= elimit) {
			gamma_prob = one - exp(arg);
		}
	}
	return gamma_prob;
}

double PPND(double P,int& IER) {
	//C
	//C ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
	//C
	//C PRODUCES NORMAL DEVIATE CORRESPONDING TO LOWER TAIL AREA = P.
	//C
	//C	See also AS 241 which contains alternative routines accurate to
	//C	about 7 and 16 decimal digits.
	//C
	double SPLIT = 0.42;
	double A[] = {2.50662823884,-18.61500062529,41.39119773534,-25.44106049637};
	double B[] = {-8.47351093090,23.08336743743,-21.06224101826,3.13082909833};
	double C[] = {-2.78718931138,-2.29796479134,4.85014127135,2.32121276858};
	double D[] = {3.54388924762,1.63706781897};
	double ZERO = 0.0, ONE = 1.0, HALF = 0.5;
	double temp, Q, R;
	
	IER = 0;
	Q = P-HALF;
	if (fabs(Q) <= SPLIT) {
		//C
		//C 0.08 < P < 0.92
		//C
		R = Q*Q;
		return Q*(((A[3]*R + A[2])*R + A[1])*R + A[0])/((((B[3]*R + B[2])*R + B[1])*R + B[0])*R + ONE);
	} else {
		//C
		//C P < 0.08 OR P > 0.92, SET R = MIN(P,1-P)
		//C
		R = P;
		if (Q > ZERO) 
			R = ONE-P;
		if (R <= ZERO) {
			IER = 1;
			return ZERO;
		}
		R = sqrt(-log(R));
		temp = (((C[3]*R + C[2])*R + C[1])*R + C[0])/((D[1]*R + D[0])*R + ONE);
		if (Q < ZERO)
			return -temp;
		else
			return temp;
	}
}

double POLY(double *c, int nord, double x) {
	//c
	//c
	//C Algorithm AS 181.2   Appl. Statist.  (1982) Vol. 31, No. 2
	//c
	//C Calculates the algebraic polynomial of order nored-1 with
	//C array of coefficients c.  Zero order coefficient is c(1)
	//c
	double temp, p;
	int j,i,n2;
	
    temp = c[0];
    if(nord == 1)
		return temp;
    p = x*c[nord-1];
    if(nord != 2) {
		n2 = nord-2;
		j = n2;
		for (i = 0; i < n2; ++i) {
			p = (p+c[j])*x;
			--j;
		}
    }
    return temp + p;
}


// X: SORTED data
// n: q. of data (3-5000)
// A: weights?

// does not work when all x identical, or when the range too small
// n2 = n % 2
// R call: init=false, X, n, n,n2,single[n2], w, pw, ifault)
double SWILK(bool INIT, double *X, int N, int N1, int N2, double *A, double& W, double &PW, int& IFAULT) {
	//
	//C
	//C ALGORITHM AS R94 APPL. STATIST. (1995) VOL.44, NO.4
	//C
	//C Calculates the Shapiro-Wilk W test and its significance level
	//C
	
	double C1[6] = {0.0E0, 0.221157E0, -0.147981E0, -0.207119E1, 0.4434685E1, -0.2706056E1};
	double C2[6] = {0.0E0, 0.42981E-1, -0.293762E0, -0.1752461E1, 0.5682633E1, -0.3582633E1};
	double C3[4] = {0.5440E0, -0.39978E0, 0.25054E-1, -0.6714E-3};
	double C4[4] = {0.13822E1, -0.77857E0, 0.62767E-1, -0.20322E-2};
	double C5[4] = {-0.15861E1, -0.31082E0, -0.83751E-1, 0.38915E-2};
	double C6[3] = {-0.4803E0, -0.82676E-1, 0.30302E-2};
	double C7[2] = {0.164E0, 0.533E0};
	double C8[2] = {0.1736E0, 0.315E0};
	double C9[2] = {0.256E0, -0.635E-2};
	double G[2] =  {-0.2273E1, 0.459E0};
	double Z90 = 0.12816E1, Z95 = 0.16449E1, Z99 = 0.23263E1;
	double ZM = 0.17509E1, ZSS = 0.56268E0;
	double ZERO = 0.0, ONE = 1.0, TWO = 2.0;
	double BF1 = 0.8378E0, XX90 = 0.556E0, XX95 = 0.622E0;
	double THREE = 3.0, SQRTH = 0.70711E0;
	double QTR = 0.25E0, TH = 0.375E0, SMALL = 1E-19;
	double PI6 = 0.1909859E1, STQR = 0.1047198E1;
	
	double SUMM2, SSUMM2, FAC, RSN, AN, AN25, A1, A2, DELTA, RANGE;
	double SA, SX, SSX, SSA, SAX, ASA, XSX, SSASSX, W1, Y, XX, XI;
	double GAMMA, M, S, LD, BF, Z90F, Z95F, Z99F, ZFM, ZSD, ZBAR;
	
	int NCENS, NN2, I, I1, J;
	bool UPPER = true;
	
	PW  =  ONE;
	if (W > ZERO) 
		W = ONE;
	AN = N;
	IFAULT = 3;
	NN2 = N/2;
	if (N2 < NN2)
		return PW;
	IFAULT = 1;
	if (N < 3)
		return PW;
	
	//C If INIT is false, calculates coefficients for the test
	if (!INIT) {
		if (N == 3)
			A[0] = SQRTH;
		else {
			AN25 = AN + QTR;
			SUMM2 = ZERO;
			for(I = 0; I < N2; ++I) {
				int itmp;
				A[I] = PPND((I + 1 - TH)/AN25,itmp);
				SUMM2 += A[I]*A[I];
			}
			SUMM2 *= TWO;
			SSUMM2 = sqrt(SUMM2);
			RSN = ONE / sqrt(AN);
			A1 = POLY(C1, 6, RSN) - A[0] / SSUMM2;
		}
		//C
		//C Normalize coefficients
		//C
		if (N > 5) {
			I1 = 3;
			A2 = -A[1]/SSUMM2 + POLY(C2,6,RSN);
			FAC = sqrt((SUMM2 - TWO * A[0]*A[0] - TWO * A[1]*A[1])/(ONE - TWO * A1*A1 - TWO * A2*A2));
			A[0] = A1;
			A[1] = A2;
		} else {
			I1 = 2;
			FAC = sqrt((SUMM2 - TWO * A[0]*A[0])/(ONE - TWO * A1*A1));
			A[0] = A1;
		}
		for (I = I1; I <= NN2; ++I)
			A[I-1] = -A[I-1]/FAC;
		INIT = true;
	}
	
	if (N1 < 3)
		return PW;
	NCENS = N - N1;
	IFAULT = 4;
	if (NCENS < 0 || (NCENS > 0 && N < 20))
		return PW;
	IFAULT = 5;
	DELTA = float(NCENS)/AN;
	if (DELTA > 0.8)
		return PW;
	//C
	//C If W input as negative, calculate significance level of -W
	//C
	
	if (W < ZERO) {
		W1 = ONE + W;
		IFAULT = 0;
	} else {
		//C
		//C Check for zero range
		//C
		IFAULT = 6;
		RANGE = X[N1-1] - X[0];
		if (RANGE < SMALL)
			return PW;
		//C
		//C Check for correct sort order on range - scaled X
		//C
		IFAULT = 7;
		XX = X[0]/RANGE;
		SX = XX;
		SA = -A[0];
		J = N;
		for (I = 2; I <= N1; ++I) {
			XI = X[I-1]/RANGE;
			if (XX-XI > SMALL)
				IFAULT=7;
			SX += XI;
			if (I != J) 
				SA += SIGN(I - J) * A[MIN(I, J)-1];
			XX = XI;
			--J;
		}
		IFAULT = 0;
		if (N > 5000) 
			IFAULT = 2;
		//C Calculate W statistic as squared correlation
		//C between data and coefficients
		SA = SA/N1;
		SX = SX/N1;
		SSA = ZERO;
		SSX = ZERO;
		SAX = ZERO;
		J = N;
		for (I = 1; I <= N1; ++I) {
			if (I != J) 
				ASA = SIGN(I - J) * A[MIN(I, J)-1] - SA;
			else
				ASA = -SA;
			
			XSX = X[I-1]/RANGE - SX;
			SSA += ASA * ASA;
			SSX += XSX * XSX;
			SAX +=  ASA * XSX;
			--J;
		}
		
		//C W1 equals (1-W) claculated to avoid excessive rounding error
		//C for W very near 1 (a potential problem in very large samples)
		SSASSX = sqrt(SSA * SSX);
		W1 = (SSASSX - SAX) * (SSASSX + SAX)/(SSA * SSX);
	}
	
	W = ONE - W1;
	//C
	//C Calculate significance level for W (exact for N=3)
	//C
	if (N == 3) {
		PW = PI6 * (asin(sqrt(W)) - STQR);
		return PW;
	}
	Y = log(W1);
	XX = log(AN);
	M = ZERO;
	S = ONE;
	if (N <= 11) {
		GAMMA = POLY(G, 2, AN);
		if (Y >= GAMMA) {
			PW = SMALL;
			return PW;
		}
		Y = -log(GAMMA - Y);
		M = POLY(C3, 4, AN);
		S = exp(POLY(C4, 4, AN));
	} else {
		M = POLY(C5, 4, XX);
		S = exp(POLY(C6, 3, XX));
	}
	if (NCENS > 0) {
		//C
		//C Censoring by proportion NCENS/N.  Calculate mean and sd
		//C of normal equivalent deviate of W.
		//C
		LD = -log(DELTA);
		BF = ONE + XX * BF1;
		Z90F = Z90 + BF * pow(POLY(C7, 2, pow(XX90,XX)), LD);
		Z95F = Z95 + BF * pow(POLY(C8, 2, pow(XX95,XX)), LD);
		Z99F = Z99 + BF * pow(POLY(C9, 2, XX),LD);
		//C
		//C Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
		//C pseudo-mean and pseudo-sd of z as the slope and intercept
		//C
		ZFM = (Z90F + Z95F + Z99F)/THREE;
		ZSD = (Z90*(Z90F-ZFM)+Z95*(Z95F-ZFM)+Z99*(Z99F-ZFM))/ZSS;
		ZBAR = ZFM - ZSD * ZM;
		M += ZBAR * S;
		S *= ZSD;
	}
	PW = alnorm((Y - M)/S, UPPER);
	return PW;
}

double chi_squared(int ndf, double chi2) {
// Calculate the chi-squared distribution function
// ndf  = number of degrees of freedom
// chi2 = chi-squared value
// prob = probability of a chi-squared value <= chi2 (i.e. the left-hand
//        tail area)
	return gammad(0.5*chi2, 0.5*ndf);
}

static double poly(const float *cc, int nord, float x)
{
/* Algorithm AS 181.2	Appl. Statist.	(1982) Vol. 31, No. 2

	Calculates the algebraic polynomial of order nord-1 with
	array of coefficients cc.  Zero order coefficient is cc(1) = cc[0]
*/
    /* Local variables */
    int j;
    double p, ret_val;/* preserve precision! */

    ret_val = cc[0];
    if (nord > 1) {
	p = x * cc[nord-1];
	for (j = nord - 2; j > 0; j--)
	    p = (p + cc[j]) * x;

	ret_val += p;
    }
    return ret_val;
} /* poly */


double ppnd(double PROB) {
	
	//c
	//c     Modified 30-Jul-91 so if PROB is zero, returns BIG as 99.9999 not 1E+38.
	//c     Some logic modified to F77 standard.
	//c
	double A0,A1,A2,A3,B1,B2,B3,B4,C0,C1,C2,C3,D1,D2,Q,R,ONE,HALF,SPLIT,ZERO, TMP;
	//C
	//C      BASED ON ALGORITHM AS 111
	//C
	double big = 0.999999e+02;
	ONE=1.0;
	HALF=0.5;
	SPLIT=0.42;
	ZERO=0.0;
	A0 =  2.50662823884;
	A1 =-18.61500062529;
	A2 = 41.39119773534;
	A3 =-25.44106049637;
	B1 = -8.47351093090;
	B2 = 23.08336743743;
	B3 =-21.06224101826;
	B4 =  3.13082909833;
	C0 = -2.78718931138;
	C1 = -2.29796479134;
	C2 =  4.85014127135;
	C3 =  2.32121276858;
	D1 =  3.54388924762;
	D2 =  1.63706781897;
	
	Q=PROB-HALF;
	if(fabs(Q) <= SPLIT) {
		R=Q*Q;
		return Q*(((A3*R+A2)*R+A1)*R+A0)/((((B4*R+B3)*R+B2)*R+B1)*R+ONE);
	} else {
		R=PROB;
		if(Q > ZERO)
			R=ONE-PROB;
		if (R > ZERO) {
            R=sqrt(-log(R));
            TMP=(((C3*R+C2)*R+C1)*R+C0)/((D2*R+D1)*R+ONE);
		} else 
            TMP=big;
		if (Q < ZERO)
			return -TMP;
		else
			return TMP;
	}
}


void
swilk(int *init,/* logical: is a[] already initialized ? */
      float *x, int *n, int *n1, int *n2,
      float *a,/* coefficients a[] */
      double *w, double *pw, int *ifault)
{

/*	ALGORITHM AS R94 APPL. STATIST. (1995) vol.44, no.4, 547-551.

	Calculates the Shapiro-Wilk W test and its significance level
*/

    /* Initialized data */
    const static float zero = 0.f;
    const static float one = 1.f;
    const static float two = 2.f;
    const static float three = 3.f;

    const static float z90 = 1.2816f;
    const static float z95 = 1.6449f;
    const static float z99 = 2.3263f;
    const static float zm = 1.7509f;
    const static float zss = .56268f;
    const static float bf1 = .8378f;
    const static double xx90 = .556;
    const static double xx95 = .622;
    const static float sqrth = .70711f;/* sqrt(1/2) = .7071068 */
    const static float small = 1e-19f;
    const static float pi6 = 1.909859f;
    const static float stqr = 1.047198f;

    /* polynomial coefficients */
    const static float g[2] = { -2.273f,.459f };
    const static float
      c1[6] = { 0.f,.221157f,-.147981f,-2.07119f, 4.434685f, -2.706056f },
      c2[6] = { 0.f,.042981f,-.293762f,-1.752461f,5.682633f, -3.582633f };
    const static float c3[4] = { .544f,-.39978f,.025054f,-6.714e-4f };
    const static float c4[4] = { 1.3822f,-.77857f,.062767f,-.0020322f };
    const static float c5[4] = { -1.5861f,-.31082f,-.083751f,.0038915f };
    const static float c6[3] = { -.4803f,-.082676f,.0030302f };
    const static float c7[2] = { .164f,.533f };
    const static float c8[2] = { .1736f,.315f };
    const static float c9[2] = { .256f,-.00635f };

    /* System generated locals */
    float r__1;

/*
	Auxiliary routines : poly()  {below}
*/
    /* Local variables */
    int i, j, ncens, i1, nn2;

    float zbar, ssassx, summ2, ssumm2, gamma, delta, range;
    float a1, a2, an, bf, ld, m, s, sa, xi, sx, xx, y, w1;
    float fac, asa, an25, ssa, z90f, sax, zfm, z95f, zsd, z99f, rsn, ssx, xsx;

    /* Parameter adjustments */
    --a;

    *pw = 1.;
    if (*w >= 0.) {
	*w = 1.;
    }
    an = (float) (*n);
    nn2 = *n / 2;
    if (*n2 < nn2) {
	*ifault = 3; return;
    }
    if (*n < 3) {
	*ifault = 1; return;
    }

/*	If INIT is false, calculate coefficients a[] for the test */
    if (! (*init)) {
	if (*n == 3) {
	    a[1] = sqrth;
	} else {
	    an25 = an + .25;
	    summ2 = zero;
	    for (i = 1; i <= *n2; ++i) {
		a[i] = (float) ppnd((i - .375f) / an25);
		r__1 = a[i];
		summ2 += r__1 * r__1;
	    }
	    summ2 *= two;
	    ssumm2 = sqrt(summ2);
	    rsn = one / sqrt(an);
	    a1 = poly(c1, 6, rsn) - a[1] / ssumm2;

	    /* Normalize a[] */
	    if (*n > 5) {
		i1 = 3;
		a2 = -a[2] / ssumm2 + poly(c2, 6, rsn);
		fac = sqrt((summ2 - two * (a[1] * a[1]) - two * (a[2] * a[2]))
			 / (one - two * (a1 * a1) - two * (a2 * a2)));
		a[2] = a2;
	    } else {
		i1 = 2;
		fac = sqrt((summ2 - two * (a[1] * a[1])) /
			   ( one  - two * (a1 * a1)));
	    }
	    a[1] = a1;
	    for (i = i1; i <= nn2; ++i)
		a[i] /= - fac;
	}
	*init = (1);
    }
    if (*n1 < 3) {
	*ifault = 1;	return;
    }
    ncens = *n - *n1;
    if (ncens < 0 || (ncens > 0 && *n < 20)) {
	*ifault = 4;	return;
    }
    delta = (float) ncens / an;
    if (delta > .8f) {
	*ifault = 5;	return;
    }

/*	If W input as negative, calculate significance level of -W */

    if (*w < zero) {
	w1 = 1. + *w;
	*ifault = 0;
	goto L70;
    }

/*	Check for zero range */

    range = x[*n1 - 1] - x[0];
    if (range < small) {
	*ifault = 6;	return;
    }

/*	Check for correct sort order on range - scaled X */

    /* *ifault = 7; <-- a no-op, since it is set 0, below, in ANY CASE! */
    *ifault = 0;
    xx = x[0] / range;
    sx = xx;
    sa = -a[1];
    j = *n - 1;
    for (i = 1; i < *n1; --j) {
	xi = x[i] / range;
	if (xx - xi > small) {
	    /* Fortran had:	 print *, "ANYTHING"
	     * but do NOT; it *does* happen with sorted x (on Intel GNU/linux):
	     *  shapiro.test(c(-1.7, -1,-1,-.73,-.61,-.5,-.24, .45,.62,.81,1))
	     */
	    *ifault = 7;
	}
	sx += xi;
	++i;
	if (i != j)
	    sa += SIGN(i - j) * a[MIN(i,j)];
	xx = xi;
    }
    if (*n > 5000) {
	*ifault = 2;
    }

/*	Calculate W statistic as squared correlation
	between data and coefficients */

    sa /= *n1;
    sx /= *n1;
    ssa = ssx = sax = zero;
    j = *n - 1;
    for (i = 0; i < *n1; ++i, --j) {
	if (i != j)
	    asa = SIGN(i - j) * a[1+MIN(i,j)] - sa;
	else
	    asa = -sa;
	xsx = x[i] / range - sx;
	ssa += asa * asa;
	ssx += xsx * xsx;
	sax += asa * xsx;
    }

/*	W1 equals (1-W) claculated to avoid excessive rounding error
	for W very near 1 (a potential problem in very large samples) */

    ssassx = sqrt(ssa * ssx);
    w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
L70:
    *w = 1. - w1;

/*	Calculate significance level for W */

    if (*n == 3) {/* exact P value : */
	*pw = pi6 * (asin(sqrt(*w)) - stqr);
	return;
    }
    y = log(w1);
    xx = log(an);
    if (*n <= 11) {
	gamma = poly(g, 2, an);
	if (y >= gamma) {
	    *pw = small;/* FIXME: rather use an even smaller value, or NA ? */
	    return;
	}
	y = -log(gamma - y);
	m = poly(c3, 4, an);
	s = exp(poly(c4, 4, an));
    } else {/* n >= 12 */
	m = poly(c5, 4, xx);
	s = exp(poly(c6, 3, xx));
    }
    /*DBG printf("c(w1=%g, w=%g, y=%g, m=%g, s=%g)\n",w1,*w,y,m,s); */

    if (ncens > 0) {/* <==>  n > n1 */

/*	Censoring by proportion NCENS/N.
	Calculate mean and sd of normal equivalent deviate of W. */

	ld = -log(delta);
	bf = one + xx * bf1;
	r__1 = pow(xx90, (double) xx);
	z90f = z90 + bf * pow(poly(c7, 2, r__1), (double) ld);
	r__1 = pow(xx95, (double) xx);
	z95f = z95 + bf * pow(poly(c8, 2, r__1), (double) ld);
	z99f = z99 + bf * pow(poly(c9, 2, xx), (double)ld);

/*	Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
	pseudo-mean and pseudo-sd of z as the slope and intercept */

	zfm = (z90f + z95f + z99f) / three;
	zsd = (z90 * (z90f - zfm) +
	       z95 * (z95f - zfm) + z99 * (z99f - zfm)) / zss;
	zbar = zfm - zsd * zm;
	m += zbar * s;
	s *= zsd;
    }
    //*pw = pnorm((double) y, (double)m, (double)s, 0/* upper tail */, 0);
    *pw  = alnorm(double((y - m)/s), true);

    return;
} /* swilk */

