#include <assert.h>
#include "svm.h"

struct CHInfo {
	long    n;			/* number of combinations */
	long	*merging;	/* merging pairs */
	long	*order;		/* merge ordering */
	double	*height;		/* merge height */
	double	ac;			/* agglomerative coefficient */
};

struct CMInfo {
	long	n;	 	 	/* number of combinations */
	long	k;	 	 	/* number of clusters */
	long	*mapping;	/* mapping */
	long	*medoids;	/* medoids */
	double *cdisp;		/* cluster dispersion */
	double disp;		/* dispersion */
};

struct CFInfo {
	long	n;	 	 	/* number of combinations */
	long	k;	 	 	/* number of clusters */
	long    iterations;
	double  value;
	long	*mapping;	/* mapping */
	double  *membership; 
	double  *cdisp;		/* cluster dispersion */
	double  disp;		/* dispersion */
};

struct LRInfo {
   int nn, k;
   double chisq;			// chi-squared
   double devnce;			// deviance
   int    ndf;				// degrees of freedom
   double *beta;			// fitted beta coefficients
   double *se_beta;			// beta std.devs
   double *fit;				// fitted probabilities for groups
   double **cov_beta;		// approx covariance matrix
   double *stdres;   		// residuals
   int    *dependent;		// dependent/redundant variables
   int	  error;
};


struct LRInput {
	long nn;
	long k;
	double **data;    //nn*k
	double *success;     //nn
	double *trials;
};


struct SVMInput {
	long nn;
	long total;
	long k;
	double **data;  //total*nn
	char   **masking;
	double *label;
};


struct SVMExample {
	long k;			 // attributes
	double *data;	 //total*nn
	char   *masking;
};

struct CInput {
	long nn;
	long jpp;
	double *data;
};

struct DInput {
	long nn;
	double *data;
};




#ifdef __cplusplus
extern "C" {
#endif

void LogReg(struct LRInput *in, struct LRInfo *OutValue);
void MCluster(struct CInput *in, long k,  int metric,  struct CMInfo *OutValue);
void HCluster(struct CInput *in,  int metric,  int method, struct CHInfo *OutValue);
void FCluster(struct CInput *in, long k,  int metric,  struct CFInfo *OutValue);
void DMCluster(struct DInput *in, long k, struct CMInfo *OutValue);
void DHCluster(struct DInput *in,  int method, struct CHInfo *OutValue);
void DFCluster(struct DInput *in, long k, struct CFInfo *OutValue);

void LRInfoCleanup(struct LRInfo *source);
void LRcleanup(struct LRInput *p);
void SVMcleanup(struct SVMInput *p);
void SVMexcleanup(struct SVMExample *p);

psvm_model SVMClassifier(struct svm_model *InValue);
struct svm_model *SVMLearn(struct SVMInput *input, int svm_type, int kernel_type, double degree,
		 double gamma, double coef0, double nu, double cache_size, double C, 
		 double eps, double p, int shrinking, int nr_weight, double *weight, 
		 double *weight_label);
double SVMClassify(psvm_model model, struct SVMExample *input);
double SVMClassifyM(psvm_model model, struct SVMExample *input);


struct XX {
	long nn;
	long k;
	double **data;    //nn*k
};
void Computer(struct XX *input, struct XX *OutValue);


#ifdef __cplusplus
}
#endif
