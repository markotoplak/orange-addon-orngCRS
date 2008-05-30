#ifndef _STAT_TESTS
#define _STAT_TESTS


#include <math.h>

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define SIGN(a) ((a) > 0.0 ? (1) : (-1))

const double EPSILON = .1192e-06;
const double INF = 1e32;

double ks2(int n, double d);
double pks2(int n, double d);
double lngamma(double z);
double log_nCr(int n, int r);
double ks1(int n, double d);
double ks2_asympt(int n, double d);
double alnorm(double x, bool upper);
double gammad(double x, double p);
double PPND(double P,int& IER);
double POLY(double *c, int nord, double x);
double chi_squared(int ndf, double chi2);


#endif