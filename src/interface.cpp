#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "interface.h"
#include "f2c.h"
#include "svm.h"

void logistic(// input
			   int& error,
			   int ngroups,			// # of examples
			   double **x,			// examples (real-valued)
			   int k,				// # attributes
			   double *s, double *n,		// s-successes of n-trials
			   
			   // output
			   double& chisq,			// chi-squared
			   double& devnce,			// deviance
			   int& ndf,				// degrees of freedom
			   double *beta,			// fitted beta coefficients
			   double *se_beta,			// beta std.devs
			   double *fit,				// fitted probabilities for groups
			   double **cov_beta,		// approx covariance matrix
			   double *stdres,			// residuals
			   int *dependent,
			   double regularization
			  );


extern "C" {
	int fanny_(integer *nn, integer *jpp, integer *kk, doublereal *x, 
		doublereal *dss, integer *jdyss, 
		doublereal *valmd, integer *jtmd, integer *ndyst, 
		integer *nsend, integer *nelem, integer *negbr, 
		doublereal *syl, doublereal *p, doublereal *dp, 
		doublereal *pt, integer *nfuzz, doublereal *esp, 
		doublereal *ef, doublereal *dvec, doublereal *ttsyl, 
		doublereal *eda, doublereal *edb, doublereal *obj, 
		integer *ncluv, doublereal *sylinf, doublereal *eps);
	
	int twins_(integer *nn,integer  *jpp, doublereal *x, doublereal *dys, doublereal *dys2,
		integer *jdyss, doublereal *valmd, integer *jtmd, integer *ndyst, integer *jalg, 
		integer *method, integer *kwan, integer *ner, doublereal *ban,doublereal  *coef,
		integer *merge);
	
	int pam_(
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
		integer *nisol);

struct svm_model *SVMLearn(struct SVMInput *input, int svm_type, int kernel_type, double degree,
		 double gamma, double coef0, double nu, double cache_size, double C, 
		 double eps, double p, int shrinking, int nr_weight=0, double *weight=NULL, 
		 double *weight_label=NULL) 
{
	struct svm_parameter param;		// set by parse_command_line
	struct svm_problem prob;		// set by read_problem
	struct svm_node *x_space;
	int elements, max_index, i, j, k;

	// default values
	param.svm_type = svm_type;
	param.kernel_type = kernel_type;
	param.degree = degree;
	param.gamma = gamma;	// 1/k
	param.coef0 = coef0;
	param.nu = nu;
	param.cache_size = cache_size;
	param.C = C;
	param.eps = eps;
	param.p = p;
	param.shrinking = shrinking;
	param.nr_weight = nr_weight;
	param.weight_label = weight_label;
	param.weight = weight;

	prob.l = input->nn;
	elements = input->total; //total is corrected + nn

	//printf("svm:1\n");

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node,elements);

	max_index = 0;
	j=0;
	for(i=0;i<prob.l;i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = input->label[i];
		for(k = 0; k < input->k; ++k)
		{
			if (input->masking[i][k] != 1) {
				x_space[j].index = k+1;
				x_space[j].value = input->data[i][k];
				++j;
			}
		}	
		if(j>=1 && x_space[j-1].index > max_index)
			max_index = x_space[j-1].index;
		x_space[j++].index = -1;
	}

	if(param.gamma == 0)
		param.gamma = 1.0/max_index;
	//printf("svm:2\n");
	return svm_train(&prob,&param);
}

psvm_model SVMClassifier(struct svm_model *InValue) {
	return InValue;
}

double SVMClassify(psvm_model model, struct SVMExample *input) {
	int i,j;
	struct svm_node *x;
	double r;

	x = (struct svm_node*)malloc(sizeof(struct svm_node)*(input->k+1));
	i = 0;
	for(j = 0; j < input->k; ++j) {
		if (input->masking[j] != 1) {
			x[i].index = j+1;
			x[i].value = input->data[j];
			++i;
		}
	}
	x[i].index = -1;

	r = svm_predict(model, x);

	free(x);

	return r;
}

double SVMClassifyM(psvm_model model, struct SVMExample *input) {
	int i,j;
	struct svm_node *x;
	double r;

	x = (struct svm_node*)malloc(sizeof(struct svm_node)*(input->k+1));
	i = 0;
	for(j = 0; j < input->k; ++j) {
		if (input->masking[j] != 1) {
			x[i].index = j+1;
			x[i].value = input->data[j];
			++i;
		}
	}
	x[i].index = -1;

	r = svm_predict_margin(model, x);

	free(x);

	return r;
}

void LogReg(struct LRInput *input, double regularization, struct LRInfo *O) {
	int i;
	
	O->nn = input->nn;
	O->k = input->k;

	//printf("#att:%d #ex:%d\n",O->k,O->nn);
	
/*
	for (i = 1; i <= input->nn; ++i) {
		printf("\nEx %d = %d :: ",i,input->success[i]);
		for (j = 1;  j<= input->k; ++j) {
			printf("%2.2f ",input->data[i][j]);
		}
	}
*/
	O->beta = (double*)malloc(sizeof(double)*(input->k+1));
	O->se_beta = (double*)malloc(sizeof(double)*(input->k+1));
	O->fit = (double*)malloc(sizeof(double)*(input->nn+1));
	O->stdres = (double*)malloc(sizeof(double)*(input->nn+1));
	O->cov_beta = (double**)malloc(sizeof(double*)*(input->k+1));
	O->dependent = (int *)malloc(sizeof(int)*(input->k+1));
	for(i = 0; i <= input->k; ++i) {
		O->cov_beta[i] = (double*)malloc(sizeof(double)*(input->k+1));
		O->dependent[i] = 0; // no dependence
	}
//	printf("+lr:2\n");
	logistic(O->error, input->nn,input->data,input->k,input->success,input->trials,
		O->chisq, O->devnce, O->ndf, O->beta, O->se_beta,
		O->fit, O->cov_beta, O->stdres, O->dependent, regularization
	);
//	printf("-lr:2\n");
	//printf("lr:3\n");
}

void LRInfoCleanup(struct LRInfo *source) {
	int i;
	
	for (i = 0; i <= source->k; ++i)
		free(source->cov_beta[i]);
	free(source->fit);
	free(source->beta);
	free(source->se_beta);
	free(source->stdres);
	free(source->cov_beta);
	free(source->dependent);
	free(source);
}

void LRcleanup(struct LRInput *p) {
	int i;
	if (p->data != NULL) {
		for (i=0; i <= p->nn; ++i)
			free(p->data[i]);
		free(p->data);
	}
	if (p->success != NULL) {
		free(p->success);
	}
	free(p);
}

void SVMexcleanup(struct SVMExample *p) {
	if (p->data != NULL) {
		free(p->data);
		free(p->masking);
	}
	free(p);
}

void SVMcleanup(struct SVMInput *p) {
	int i;
	if (p->data != NULL) {
		for (i=0; i <= p->nn; ++i) {
			free(p->data[i]);
			free(p->masking[i]);
		}
		free(p->data);
		free(p->masking);
	}
	if (p->label!= NULL) {
		free(p->label);
	}
	free(p);
}

void HCluster(struct CInput *input, int metric, int method, struct CHInfo *OutValue) {
	long i, res,nn,jpp;
	long *kwan, *ner, jdyss, nydist, *jmtd, meth, jalg, *merge;
	double *ban, *x, *dv2, *dv, *valmd, ac; // matrix

	nn = input->nn;
	jpp = input->jpp;
	x = input->data;

	dv = (double*)malloc((1 + (nn * (nn - 1))/2)*sizeof(double));
	dv2 = (double*)malloc((1 + (nn * (nn - 1))/2)*sizeof(double));
	//x = (double*)malloc(sizeof(double)*(nn*jpp));
	merge = OutValue->merging = (long*)malloc(sizeof(long)*(nn-1)*2);
	jmtd = (long*)malloc(sizeof(long)*jpp);
	valmd = (double*)malloc(sizeof(double)*jpp);
	kwan = (long*)malloc(sizeof(long)*nn);
	ner = OutValue->order = (long*)malloc(sizeof(long)*nn);
	ban = OutValue->height = (double*)malloc(sizeof(double)*nn);

	OutValue->n = nn;

	//nydist = 2; // 1-euclidean, 2-manhattan
	nydist = metric;
	//meth = 4;   1-average, 2-single, 3-complete, 4-ward, 5-weighted
	meth = method;


	jalg = 1;
	jdyss = 0;

	// valmisdat -- vrednost, ki oznacuje miss val
	for (i = 0; i < jpp; ++i) {
		jmtd[i] = 1; // no missing values
		valmd[i] = -0.5;
	}
	res = twins_(&nn, &jpp, x, dv, dv2, &jdyss, valmd, jmtd, &nydist, &jalg, &meth,
		         kwan, ner, ban, &ac, merge);

	OutValue->ac = ac;

//	cout << ac;
/*	cout << "\nMerging: \n";
	for (i = 0; i < nn-1; ++i) {
		cout << merge[i*2] << " " << merge[i*2+1] << "\n";
	}
	cout << "\nOrder: ";
	for (i = 0; i < nn; ++i) {
		cout << ner[i] << " ";
	}
	cout << "\nHeight: ";
	for (i = 0; i < nn; ++i) {
		cout << ban[i] << " ";
	}
	*/
	free(dv); free(dv2); free(jmtd); free(valmd); free(kwan);
}


void MCluster(struct CInput *input, long k, int metric, struct CMInfo *OutValue) {
	long i, res, nn, jpp;
	long jdyss, nydist, *jmtd, *a1, *a2, *a3, *med;
	double *x, *dv, *valmd, *b1, *b2, *b3, *avsil; // matrix
	double *obj, *clusinf, *silinf;
	long *clu, *isol;

	nn = OutValue->n = input->nn;
	assert(nn > 1);
	jpp = input->jpp;
	x = input->data;

	clu = OutValue->mapping = (long*)malloc(sizeof(long)*nn);
	med = OutValue->medoids = (long*)malloc(sizeof(long)*k);

	dv = (double*)malloc((1 + (nn * (nn - 1))/2)*sizeof(double));
	//x = (double*)malloc(sizeof(double)*(nn*jpp));
	jmtd = (long*)malloc(sizeof(long)*jpp);
	valmd = (double*)malloc(sizeof(double)*jpp);
	a1 = (long*)malloc(sizeof(long)*nn);
	a2 = (long*)malloc(sizeof(long)*nn);
	a3 = (long*)malloc(sizeof(long)*nn);
	b1 = (double*)malloc(sizeof(double)*nn);
	b2 = (double*)malloc(sizeof(double)*nn);
	b3 = (double*)malloc(sizeof(double)*nn);
	avsil = (double*)malloc(sizeof(double)*nn);
	obj = (double*)malloc(sizeof(double)*2);
	clusinf = (double*)malloc(sizeof(double)*5*k);
	silinf = (double*)malloc(sizeof(double)*4*nn);
	isol = (long*)malloc(sizeof(long)*k);

	/*
	for (i = 0; i < nn; ++i) {
		for (j = 0; j < jpp; ++j) {
			x[i+nn*j] = input[i][j];
		}
	}
	*/

	nydist = metric;
	//nydist = 2;

	jdyss = 0;
	// valmisdat -- vrednost, ki oznacuje miss val
	for (i = 0; i < jpp; ++i) {
		jmtd[i] = 1; // no missing values
		valmd[i] = -0.5;
	}

	res = pam_(&nn, &jpp, &k, 
		x, dv, &jdyss, valmd, jmtd, &nydist, 
		a1, a2, a3, b1, b2, avsil, b3, &(OutValue->disp),
		med, obj, clu, clusinf, silinf, isol
	);


/*	entropy = 0;
	for (i = 0; i < k; ++i) {
		sum = 0;
		for (j = 0; j < jpp; ++j) {
			if ((*cases[med[i]-1]->a)[j] > 1e-6)
				entropy -= (*cases[med[i]-1]->a)[j] * log((*cases[med[i]-1]->a)[j]);
			sum += (*cases[med[i]-1]->a)[j];
		}
		if ((1-sum) > 1e-6)
			entropy -= (1-sum) * log(1-sum);
	}
*/
	// copy the average silhouettes
	OutValue->cdisp = (double*)malloc(sizeof(double)*k);
	OutValue->k = k;
	for (i = 0; i < k; ++i) 
		OutValue->cdisp[i] = avsil[i];

	free(dv); 
	free(jmtd); 
	free(valmd); 
	free(a1); 
	free(a2);
	free(a3); 
	free(b1); 
	free(b2); 
	free(b3); 
	free(avsil);
	free(obj);
	free(clusinf);
	free(silinf);
	free(isol);
}

void FCluster(struct CInput *input, long k, int metric, struct CFInfo *OutValue) {
	long i, res, nn, jpp;
	long jdyss, nydist, *jmtd, *a1, *a2, *a3, *ik;
	double *x, *dv, *valmd, *b1, *dn, *dk1, *dk2, *avsil, *p, *m1; // matrix
	double *obj, *silinf;
	double eps, eda, edb, ttsil;
	long *clu;

	nn = OutValue->n = input->nn;
	OutValue->k = k;
	assert(nn > 1);
	jpp = input->jpp;
	x = input->data;

	clu = OutValue->mapping = (long*)malloc(sizeof(long)*nn);

	dv = (double*)malloc((1 + (nn * (nn - 1))/2)*sizeof(double));
	//x = (double*)malloc(sizeof(double)*(nn*jpp));
	jmtd = (long*)malloc(sizeof(long)*jpp);
	valmd = (double*)malloc(sizeof(double)*jpp);
	a1 = (long*)malloc(sizeof(long)*nn);
	a2 = (long*)malloc(sizeof(long)*nn);
	a3 = (long*)malloc(sizeof(long)*nn);
	b1 = (double*)malloc(sizeof(double)*nn);
	dn = (double*)malloc(sizeof(double)*nn);
	dk1 = (double*)malloc(sizeof(double)*k);
	dk2 = (double*)malloc(sizeof(double)*k);
	ik = (long*)malloc(sizeof(long)*k);
	OutValue->cdisp = avsil = (double*)malloc(sizeof(double)*k);
	obj = (double*)malloc(sizeof(double)*2);
	silinf = (double*)malloc(sizeof(double)*4*nn);

	OutValue->membership = p = (double*)malloc(sizeof(double)*nn*k);
	m1 = (double*)malloc(sizeof(double)*nn*k);

	eps = 1e-15;
	eda = edb = 0.0;

	ttsil = 0;
	nydist = metric;

	jdyss = 0;
	// valmisdat -- vrednost, ki oznacuje miss val
	for (i = 0; i < jpp; ++i) {
		jmtd[i] = 1; // no missing values
		valmd[i] = 1; // not missing, -1 for missing
	}

	res = fanny_(&nn, &jpp, &k, 
		x, dv /**/, &jdyss /**/, valmd, jmtd, 
		&nydist, 
		a1, a2, a3, b1, 
		p, m1,
		avsil, 
		ik,dk1,dk2,dn,
		&ttsil, 
		&eda /**/, &edb /**/,
		obj /**/, clu /**/, 
		silinf /**/, &eps
	);


/*	entropy = 0;
	for (i = 0; i < k; ++i) {
		sum = 0;
		for (j = 0; j < jpp; ++j) {
			if ((*cases[med[i]-1]->a)[j] > 1e-6)
				entropy -= (*cases[med[i]-1]->a)[j] * log((*cases[med[i]-1]->a)[j]);
			sum += (*cases[med[i]-1]->a)[j];
		}
		if ((1-sum) > 1e-6)
			entropy -= (1-sum) * log(1-sum);
	}
*/
	
	// copy the average silhouettes
	OutValue->iterations =  int(obj[0]);
	OutValue->value = obj[1];

	free(dv); 
	free(jmtd); 
	free(valmd); 
	free(a1); 
	free(a2);
	free(a3); 
	free(b1); 
	free(m1);
	
	free(dk1);
	free(dk2);
	free(ik);
	free(dn);

	free(obj);
	free(silinf);
}

void DHCluster(struct DInput *input, int method, struct CHInfo *OutValue) {
	long res,nn,jpp;
	long *kwan, *ner, jdyss, nydist, meth, jalg, *merge;
	double *ban, *x, *dv2, *dv, ac; // matrix
	double valmd;
	long jmtd;

	nn = input->nn;
	jpp = 1;
	dv = input->data;

	dv2 = (double*)malloc((1 + (nn * (nn - 1))/2)*sizeof(double));
	x = (double*)malloc(sizeof(double)*nn);
	merge = OutValue->merging = (long*)malloc(sizeof(long)*(nn-1)*2);
	kwan = (long*)malloc(sizeof(long)*nn);
	ner = OutValue->order = (long*)malloc(sizeof(long)*nn);
	ban = OutValue->height = (double*)malloc(sizeof(double)*nn);

	OutValue->n = nn;

	//nydist = 2; // 1-euclidean, 2-manhattan
	nydist = 0;
	//meth = 4;   1-average, 2-single, 3-complete, 4-ward, 5-weighted
	meth = method;


	jalg = 1;
	jdyss = 1;

	res = twins_(&nn, &jpp, x, dv, dv2, &jdyss, &valmd, &jmtd, &nydist, &jalg, &meth,
		         kwan, ner, ban, &ac, merge);

	OutValue->ac = ac;

//	cout << ac;
/*	cout << "\nMerging: \n";
	for (i = 0; i < nn-1; ++i) {
		cout << merge[i*2] << " " << merge[i*2+1] << "\n";
	}
	cout << "\nOrder: ";
	for (i = 0; i < nn; ++i) {
		cout << ner[i] << " ";
	}
	cout << "\nHeight: ";
	for (i = 0; i < nn; ++i) {
		cout << ban[i] << " ";
	}
	*/
	free(dv2); free(kwan); free(x);
}


void DMCluster(struct DInput *input, long k, struct CMInfo *OutValue) {
	long i, res, nn, jpp;
	long jdyss, nydist, jmtd, *a1, *a2, *a3, *med;
	double *x, *dv, valmd, *b1, *b2, *b3, *avsil; // matrix
	double *obj, *clusinf, *silinf;
	long *clu, *isol;

	nn = OutValue->n = input->nn;
	jpp = 1;

	clu = OutValue->mapping = (long*)malloc(sizeof(long)*nn);
	med = OutValue->medoids = (long*)malloc(sizeof(long)*k);

	dv = input->data;
	x = (double*)malloc(sizeof(double)*nn);
	a1 = (long*)malloc(sizeof(long)*nn);
	a2 = (long*)malloc(sizeof(long)*nn);
	a3 = (long*)malloc(sizeof(long)*nn);
	b1 = (double*)malloc(sizeof(double)*nn);
	b2 = (double*)malloc(sizeof(double)*nn);
	b3 = (double*)malloc(sizeof(double)*nn);
	avsil = (double*)malloc(sizeof(double)*nn);
	obj = (double*)malloc(sizeof(double)*2);
	clusinf = (double*)malloc(sizeof(double)*5*k);
	silinf = (double*)malloc(sizeof(double)*4*nn);
	isol = (long*)malloc(sizeof(long)*k);

	/*
	for (i = 0; i < nn; ++i) {
		for (j = 0; j < jpp; ++j) {
			x[i+nn*j] = input[i][j];
		}
	}
	*/

	nydist = 0;

	jdyss = 1;

	res = pam_(&nn, &jpp, &k, 
		x, dv, &jdyss, &valmd, &jmtd, &nydist, 
		a1, a2, a3, b1, b2, avsil, b3, &(OutValue->disp),
		med, obj, clu, clusinf, silinf, isol
	);


/*	entropy = 0;
	for (i = 0; i < k; ++i) {
		sum = 0;
		for (j = 0; j < jpp; ++j) {
			if ((*cases[med[i]-1]->a)[j] > 1e-6)
				entropy -= (*cases[med[i]-1]->a)[j] * log((*cases[med[i]-1]->a)[j]);
			sum += (*cases[med[i]-1]->a)[j];
		}
		if ((1-sum) > 1e-6)
			entropy -= (1-sum) * log(1-sum);
	}
*/
	// copy the average silhouettes
	OutValue->cdisp = (double*)malloc(sizeof(double)*k);
	OutValue->k = k;
	for (i = 0; i < k; ++i) 
		OutValue->cdisp[i] = avsil[i];

	free(a1); free(a2);
	free(a3); free(b1); free(b2); free(b3); free(avsil);free(obj);
	free(clusinf);free(silinf);free(isol);

	free(x);
}

void DFCluster(struct DInput *input, long k, struct CFInfo *OutValue) {
	long i, res, nn, jpp;
	long jdyss, nydist, *jmtd, *a1, *a2, *a3, *ik;
	double *x, *dv, *valmd, *b1, *dn, *dk1, *dk2, *avsil, *p, *m1; // matrix
	double *obj, *silinf;
	double eps, eda, edb, ttsil;
	long *clu;

	nn = OutValue->n = input->nn;
	OutValue->k = k;
	assert(nn > 1);
	jpp = 1;

	clu = OutValue->mapping = (long*)malloc(sizeof(long)*nn);

	dv = (double*)malloc((1 + (nn * (nn - 1))/2)*sizeof(double));
	x = (double*)malloc(sizeof(double)*(nn));
	jmtd = (long*)malloc(sizeof(long)*jpp);
	valmd = (double*)malloc(sizeof(double)*jpp);
	a1 = (long*)malloc(sizeof(long)*nn);
	a2 = (long*)malloc(sizeof(long)*nn);
	a3 = (long*)malloc(sizeof(long)*nn);
	b1 = (double*)malloc(sizeof(double)*nn);
	dn = (double*)malloc(sizeof(double)*nn);
	dk1 = (double*)malloc(sizeof(double)*k);
	dk2 = (double*)malloc(sizeof(double)*k);
	ik = (long*)malloc(sizeof(long)*k);
	OutValue->cdisp = avsil = (double*)malloc(sizeof(double)*k);
	obj = (double*)malloc(sizeof(double)*2);
	silinf = (double*)malloc(sizeof(double)*4*nn);

	OutValue->membership = p = (double*)malloc(sizeof(double)*nn*k);
	m1 = (double*)malloc(sizeof(double)*nn*k);

	eps = 1e-15;
	eda = edb = 0.0;

	nydist = 0;
	
	// CORRECT THE DISPLACEMENT: it's different from PAM and AGNES
	for (i = 0; i < (nn * (nn - 1))/2; ++i) {
		dv[i] = input->data[(nn * (nn - 1))/2 - i];
	}

	jdyss = 1;
	// valmisdat -- vrednost, ki oznacuje miss val
	for (i = 0; i < jpp; ++i) {
		jmtd[i] = 1; // no missing values
		valmd[i] = 1;
	}

	res = fanny_(&nn, &jpp, &k, 
		x, dv /**/, &jdyss /**/, valmd, jmtd, 
		&nydist, 
		a1, a2, a3, b1, 
		p, m1,
		avsil, 
		ik,dk1,dk2,dn,
		&ttsil, 
		&eda /**/, &edb /**/,
		obj /**/, clu /**/, 
		silinf /**/, &eps
	);


/*	entropy = 0;
	for (i = 0; i < k; ++i) {
		sum = 0;
		for (j = 0; j < jpp; ++j) {
			if ((*cases[med[i]-1]->a)[j] > 1e-6)
				entropy -= (*cases[med[i]-1]->a)[j] * log((*cases[med[i]-1]->a)[j]);
			sum += (*cases[med[i]-1]->a)[j];
		}
		if ((1-sum) > 1e-6)
			entropy -= (1-sum) * log(1-sum);
	}
*/
	
	// copy the average silhouettes
	OutValue->iterations =  int(obj[0]);
	OutValue->value = obj[1];

	free(dv); 
	free(x);
	free(jmtd); 
	free(valmd); 
	free(a1); 
	free(a2);
	free(a3); 
	free(b1); 
	free(m1);
	
	free(dk1);
	free(dk2);
	free(ik);
	free(dn);

	free(obj);
	free(silinf);
}

void xComputer(struct XX *input, struct XX *OutValue) {
	int i, j, k, mi,ma,ofs;
	double a,b,c,d, x,y;
	double *dd, *ee;


	ma = (input->k*(input->k-1)/2);
	OutValue->data = (double **)malloc(sizeof(double*)*ma);
	OutValue->nn = ma;
	OutValue->k = 4;
	for(i = 0; i < ma; ++i) {
		OutValue->data[i] = (double *)malloc(sizeof(double)*4);
		for(j = 0; j < 4; ++j)
			OutValue->data[i][j] = 0.0;
	}
	for(k = 0; k < input->nn; ++k) {
		ofs = 0;
		dd = input->data[k];
		for(i = 1; i < input->k; ++i) {
			x = dd[i];
			for (j = 0; j < i; ++j) {
				y = dd[j];
				ee = OutValue->data[ofs];
				if(x > y) {
					if(x < 1)
						ee[0] += 1-x;
					ee[2] += x-y;
					ee[3] += y;
				} else {
					if(y < 1)
						ee[0] += 1-y;
					ee[1] += y-x;
					ee[3] += x;
				}
				++ofs;
			}
		}
	}
}





void Computer(struct XX *input, struct XX *OutValue) {
	int i, j, k, mi,ma,ofs;
	double x,y,xx,yy, r;
	double *dd, *ee, *means, *devs, *data;


	ma = (input->k*(input->k-1)/2);
	OutValue->data = (double **)malloc(sizeof(double*));
	OutValue->nn = 1; // 0 - results, 1 - means, 2 - standard deviations
	OutValue->k = ma;
	for(i = 0; i < 1; ++i) {
		OutValue->data[i] = (double *)malloc(sizeof(double)*ma);
		for(j = 0; j < ma; ++j) {
			OutValue->data[i][j] = 0.0;
		}
	}
	data = OutValue->data[0];

	// means
	means = (double *)malloc(sizeof(double)*input->k);
	devs = (double *)malloc(sizeof(double)*input->k);
	for(i = 0; i < input->k; ++i) {
		means[i] = 0.0;
		devs[i] = 0.0;
	}
	for(k = 0; k < input->nn; ++k) {
		dd = input->data[k];
		for(i = 0; i < input->k; ++i) {
			means[i] += dd[i];
		}
	}
	for(i = 0; i < input->k; ++i) {
		means[i] /= input->nn;
	}
	// variances
	for(k = 0; k < input->nn; ++k) {
		dd = input->data[k];
		for(i = 0; i < input->k; ++i) {
			r = dd[i]-means[i];
			devs[i] += r*r;
		}
	}
	// correlations
	for(k = 0; k < input->nn; ++k) {
		ofs = 0;
		dd = input->data[k];
		for(i = 1; i < input->k; ++i) {
			x = dd[i];
			xx = x-means[i];
			for (j = 0; j < i; ++j) {
				y = dd[j];
				yy = y - means[j];
				data[ofs] += xx*yy;
				++ofs;
			}
		}
	}
	// normalization
	ofs = 0;
	for(i = 1; i < input->k; ++i) {
		for (j = 0; j < i; ++j) {
			data[ofs] /= sqrt(devs[i]*devs[j]);
			++ofs;
		}
	}

	free(means);
	free(devs);
}





}