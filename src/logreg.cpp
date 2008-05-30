#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "lsq.h"
#include "pr_statistics.h"



void disaster() {
}


void logistic(// input
 			   int& ier,
			   int ngroups,			// # of examples
			   double **x,			// examples (real-valued)
			   int k,				// # attributes
			   double *s, double *n,// s-successes of n-trials
			   
			   // output
			   double& chisq,			// chi-squared
			   double& devnce,			// deviance
			   int& ndf,				// degrees of freedom
			   double *beta,			// fitted beta coefficients
			   double *se_beta,			// beta std.devs
			   double *fit,				// fitted probabilities for groups
			   double **cov_beta,		// approx covariance matrix
			   double *stdres,		 	// residuals
			   int *dependent,
			   double regularization
			   )
{
	int i, iter, j, ncov, pos;
	double *propn, *p, *wt, *xrow, TOL,
		*db, *bnew, dev_new, xb, *pnew, *covmat,
		*wnew, *range, var, *e, hii;
	bool *lindep;
	double penalization;
	lsq fitter;

	TOL = 0.0001; // tolerance
	
	ier = 0;
	ndf = ngroups - k - 1;
	if (ngroups < 2 || ndf < 0) {
		ier = 1; // not enough examples
		return;
	}
	for(i = 1; i <= ngroups; ++i) {
		if (s[i] < 0) {
			ier = 2;
			return;
		}
		if (n[i] < 0) {
			ier = 3;
			return;
		}
		if (n[i] < s[i]) {
			ier = 4;
			return;
		}
	}
	range = (double *)malloc(sizeof(double)*(k+1));
	for (i = 1; i <= k; ++i) {
		double min = INF;
		double max = -INF;
		
		for (j = 1; j <= ngroups; ++j) {
			if (x[j][i] < min)
				min = x[j][i];
			if (x[j][i] > max)
				max= x[j][i];
		}
		range[i] = max-min;
		if (range[i] < EPSILON*(fabs(min)+fabs(max))) {
			free(range);
			ier = 5; // constant variable
			//printf("variable %d is constant\n",i);
			dependent[i] = 1;
			return;
		}
	}
	
	++ngroups; ++k;
	propn = (double *)malloc(sizeof(double)*ngroups);
	p = (double *)malloc(sizeof(double)*ngroups);
	pnew = (double *)malloc(sizeof(double)*ngroups);
	wnew = (double *)malloc(sizeof(double)*ngroups);
	e = (double *)malloc(sizeof(double)*ngroups);
	wt = (double *)malloc(sizeof(double)*ngroups);
	
	xrow = (double *)malloc(sizeof(double)*k);
	db = (double *)malloc(sizeof(double)*k);	
	bnew = (double *)malloc(sizeof(double)*k);
	lindep = (bool *)malloc(sizeof(bool)*k);
	--k; --ngroups;
	
	for(i = 1; i <= ngroups; ++i) {
		/*
		printf("\n%2.2f %2.2f: ",s[i],n[i]);
		for (j = 1; j <= k; ++j) {
			printf("%5.2f ",x[i][j]);
		}*/
		propn[i] = s[i]/n[i];
		wt[i] = 1.0;
		p[i] = 0.5;
	}
	
	for(i = 0; i <= k; ++i) {
		beta[i] = 0.0;
	}
	
	iter = 1;
	
	do {
		fitter.startup(k,true);
		for (i = 1; i <= ngroups; ++i) {
			if (iter == 0) {
				xrow[0] = 0.25;
				for (j = 1; j <= k; ++j) {
					xrow[j] = 0.25*x[i][j];
				}
			} else {
				xrow[0] = p[i]*(1.0-p[i]);
				for (j = 1; j <= k; ++j) {
					xrow[j] = p[i]*(1.0-p[i])*x[i][j];
				}
			}
			fitter.includ(wt[i], xrow-1, propn[i]-p[i]);
		}
		
		//! Test for a singularity
		fitter.sing(lindep-1, ier);
		if (ier != 0) {
			for( i = 1; i<= k; ++i) {
				if (lindep[i]) {
					dependent[i] = 1;
					//printf("Variable number %d is linearly dependent upon earlier variables\n",i);
				}
			}
			ier = 6;
			return;
		}
		fitter.regcf(db-1, k+1, ier); // corrected
l10: 
		for (i = 0; i <= k; ++i) {
			bnew[i] = beta[i]+db[i];
		}
		
		//! Calculate new p(i)'s, weights & deviance
		
		dev_new = 0.0;
		for (i = 1; i <= ngroups; ++i) {
			xb = bnew[0];
			for (j = 1; j <= k; ++j) {
				xb += x[i][j] * bnew[j];
			}
			xb = exp(xb);
			pnew[i] = xb / (1.0 + xb);
			wnew[i] = (n[i])/(pnew[i]*(1.0 - pnew[i]));
			if (iter == 1) 
				wnew[i] = sqrt(wnew[i]);
			if (s[i] > 0) 
				dev_new = dev_new + s[i]*log(propn[i]/pnew[i]);
			if (s[i] < n[i]) 
				dev_new += (n[i]-s[i])*log((1.0-propn[i])/(1.0-pnew[i]));
		}
		dev_new = 2 * dev_new;
		//! If deviance has increased, reduce the step size.
		
		if (iter > 2 && dev_new > devnce*(1+TOL)) {
			for (i = 0; i <= k; ++i)
				db[i] = 0.5 * db[i];
			goto l10;
		}
		//! Replace betas, weights & p's with new values
		
		for (i = 0; i <= k; ++i) {
			beta[i] = bnew[i];
		}
		for (i = 1; i <= ngroups; ++i) {
			wt[i] = wnew[i];
			p[i] = pnew[i];
		}
		//! Test for convergence
		
//		printf("iter. %d, dev: %f\n",iter,devnce-dev_new);
		if (iter > 2 && devnce - dev_new < TOL)
			break;
		devnce = dev_new;
		iter = iter + 1;
		if (iter%20 == 0) {
			TOL *= 2;
		}
		if (iter > 2000) {
			ier = 8;
			return;
		}				
		//! Test for a very large beta
		
		for( i = 1; i<= k; ++i) {
			if(fabs(beta[i])*range[i] > 100) {
//				dependent[i] = 1;
//				printf("Coefficient for variable no. %d tending to infinity",i);
				
				ier = 7;
				return;
			}
		}
	} while(true);

	chisq = 0.0;
	for (i = 1; i <= ngroups; ++i) {
		double tt;

		e[i] = n[i]*p[i];
		tt = s[i]-e[i];
		chisq += tt*tt/e[i];
	}
	devnce = dev_new;

//! Calculate the approximate covariance matrix for the beta's, if ndf > 0.

	covmat = NULL;
	if (ndf > 0) {
		ncov = (k+1)*(k+2)/2;
		covmat = (double *)malloc(sizeof(double)*(ncov+1));
		fitter.cov(k+1, var, covmat, ncov, se_beta-1, ier);
		if (var < 1.0) {
			for (i = 1; i <= ncov; ++i) {
				covmat[i] /= var;
			}
			for (i = 0; i <= k; ++i) {
			    se_beta[i] /= sqrt(var);
			}
		}
		if (cov_beta != NULL) {
			pos = 1;
			for(i = 0;i<= k; ++i) {
				cov_beta[i][i] = covmat[pos];
				pos = pos + 1;
				for(j = i+1; j<= k; ++j) {
					cov_beta[i][j] = covmat[pos];
					cov_beta[j][i] = covmat[pos];
					pos = pos + 1;
				}
			}
		}
	}
	
	if (fit != NULL) {
		for (i = 1; i <= ngroups; ++i) {
			fit[i] = p[i];
		}
	}

	if (stdres != NULL) {
		for (i = 1; i <= ngroups; ++i) {
			xrow[0] = p[i]*(1.0 - p[i]);
			for (j = 1; j <= k; ++j) {
				xrow[j] = p[i]*(1.0 - p[i])*x[i][j];
			}
			fitter.hdiag(xrow-1, k+1, hii, ier);
			stdres[i] = (s[i]-n[i]*p[i]) / sqrt(n[i]*p[i]*(1.0-p[i])*(1.0-hii));
		}
	}

	if (covmat != NULL)
		free(covmat);
	free(propn); 
	free(p); 
	free(wnew); 
	free(pnew); 
	free(e); 
	free(wt);	
	free(xrow); 
	free(db); 
	free(bnew); 
	free(range); 
	free(lindep);
}
