#include <assert.h>
#include "svm.h"

typedef struct svm_model *psvm_model;

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

struct SVMSparseInput {
	long nn;
	long elements;
	int *lengths;
	double **value;  
	int    **index;
	double *label;
};


struct SVMSparseExample {
	long nn;
	double *value;
	int	   *index;
};

struct SVMOut {
	long nn;
	double *v;
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

struct wsvm_model {
	struct svm_model *m;
	struct svm_node *x_space;
	struct svm_problem *prob;
};



#ifdef __cplusplus
extern "C" {
#endif

void LogReg(struct LRInput *in, double regularization, struct LRInfo *OutValue);
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
void SVMscleanup(struct SVMSparseInput *p);
void SVMsecleanup(struct SVMSparseExample *p);

double SVMClassify(psvm_model model, struct SVMExample *input);
void SVMClassifyP(psvm_model model, struct SVMExample *input, struct SVMOut *OutValue );
void SVMClassifyM(psvm_model model, struct SVMExample *input, struct SVMOut *OutValue );
double SVMClassifyS(psvm_model model, struct SVMSparseExample *input);
void SVMClassifyPS(psvm_model model, struct SVMSparseExample *input, struct SVMOut *OutValue );
void SVMClassifyMS(psvm_model model, struct SVMSparseExample *input, struct SVMOut *OutValue );
struct wsvm_model *SVMLearnS(struct SVMSparseInput *input, int svm_type, int kernel_type, double degree,
		 double gamma, double coef0, double nu, double cache_size, double C, 
		 double eps, double p, int shrinking, int probability, int nr_weight, double *weight, 
		 int *weight_label);
struct wsvm_model *SVMLearn(struct SVMInput *input, int svm_type, int kernel_type, double degree,
		 double gamma, double coef0, double nu, double cache_size, double C, 
		 double eps, double p, int shrinking, int probability, int nr_weight, double *weight, 
		 int *weight_label);

void svm_destroy_model(struct svm_model *model);
psvm_model SVMClassifier(struct svm_model *InValue);



struct XX {
	long nn;
	long k;
	double **data;    //nn*k
};
void Computer(struct XX *input, struct XX *OutValue);


#ifdef __cplusplus
}
#endif
