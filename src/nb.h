#ifndef _NB_H
#define _NB_H

#pragma warning( disable : 4786 )

#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct NBGroup {
	int n;
	int *l;
};

struct NBList {
	int n;
	double *l;
};

struct NBModel {
	int nosingles;
	int *singles;
	int nogroups;
	struct NBGroup *groups;
};

struct NBResult {
	double kl_q;
	double kl_err;
	double b_q;
	double b_err;
	double er_q;
	double er_err;
};

struct NBInput {
	long nn;	   // number of examples
	long k;
	double a;
	long na;	   // number of attributes (k-1)
	int **data;    // nn*k
	int *l;        // nn
	int *card;	   // k
	int maxcvi;    // number of folds
	int *cvi;      // nn
};

struct NBInfo {
	struct NBInput *i;
	void *c;
};

void NBcleanup(struct NBInput *p);
void NBquality(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue);
void TANquality(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue);
struct NBInfo *NBprepare(struct NBInput *input);
void NBkill(struct NBInfo *p);

void NBsaveScores(struct NBInfo *in);
void NBrememberScores(struct NBInfo *in);
void NBcompareScores(struct NBInfo *in, struct NBResult *OutValue);
double NBcompareLists(int n, double *a, double *b);

void NBclassify(struct NBInfo *in, int *ex, struct NBList *OutValue);
void NBclassifyW(struct NBInfo *in, int *ex, struct NBList *OutValue);

void NBexportScores(struct NBInfo *in, int mode, struct NBList *OutValue);
void NBexportProbabilities(struct NBInfo *in, int mode, struct NBList *OutValue);
void NBstoreModel(struct NBInfo *in, double *w, struct NBModel *m);
void NBcompareScores(struct NBInfo *in, struct NBResult *OutValue);
void NBclassify(struct NBInfo *in, int *ex, struct NBList *OutValue);
void NBclassifyW(struct NBInfo *in, int *ex, struct NBList *OutValue);
void NBquality(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue);
void NBqualityW(struct NBInfo *in, double *w, struct NBModel *m, struct NBResult *OutValue);
void NBupdate(struct NBInfo *in, int attribute, int card, int *values);

void NBdivergence(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue);


#endif