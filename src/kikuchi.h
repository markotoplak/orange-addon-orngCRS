#ifndef _KIKUCHI_H
#define _KIKUCHI_H

#pragma warning( disable : 4786 )

#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct KGroup {
	int n;
	int times;
	int *l;
};

struct KList {
	int n;
	double *l;
};

struct KMatrix {
	int rows;
	int columns;
	double *l;
};

struct KModel {
	int nogroups;
	struct KGroup *groups;
};

struct KModels {
	int nomodels;
	struct KModel *models;
};

struct KInput {
	long nn;	   // number of examples
	long k;		   // number of information sources
	long na;	   // number of attributes (k-1)
	int **data;    // nn*k
	int *card;	   // k, cardinalities
};

struct KInfo {
	struct KInput *i;
	void *c;
};

void Kdie(struct KInfo *p);
struct KInfo *Kremember(struct KInput *input);
void Klearn(struct KInfo *in, int samples, int depth, struct KMatrix *OutValue);
void Kcleanup(struct KInput *p);

void Ktestmodels(struct KInfo *in, struct KModels *ms, int samples, struct KMatrix *OutValue);
void Kprepare(struct KInfo *in, double prior, struct KModel *m);
void Kuse(struct KInfo *in, int *ex, struct KList *OutValue);


#endif