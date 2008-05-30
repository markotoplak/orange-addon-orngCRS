#ifndef _KIKUCHI_H
#define _KIKUCHI_H

#pragma warning( disable : 4786 )

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <map>
#include "cMersenneTwister.h"

#define MAX_INTS (100000)


using namespace std;

struct KGroup {
	int n;
	int times;
	int *l;
};

struct KList {
	int n;
	double *l;
};

struct KArray {
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


struct KCModel {
	int n;
	int na;
	int *indices;
	int *magnitudes;
	int *stats;
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

#define safelog(prob) ((prob >= 1e-8)?(log(prob)):(-18.5))

class Stats {
public:
	int *v; // frequencies
	int total;
	int labelled;
	long allcomb;
	double divisor;
	//double xdivisor;
	//double depriorify;
	double prior;

	void model();
	void reset();
	int getDOF();
	int getDOFsup(int card);
	void set_unlabelled();
	void set_prior(double _p);
	virtual void update(int *ex) = 0;
	virtual double getprob(int *ex ) = 0;
	//virtual double getxprob(int *ex );
	Stats(long _allcomb);
	~Stats();
};


class Stats1 : public Stats {
public:
	int a;
	int nvaluesa;

	virtual void update(int *ex);
	virtual double getprob(int *ex);
	Stats1(int _a, int _nvalues);
};

class Stats2 : public Stats {
public:
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;

	Stats2(int _a, int _b, int _nvaluesa, int _nvaluesb);
	virtual void update(int *ex);
	virtual double getprob(int *ex);
};

class Stats3 : public Stats{
public:
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;
	int c,m;
	int nvaluesc;

	Stats3(int _a, int _b, int _c, int _nvaluesa, int _nvaluesb, int _nvaluesc);
	virtual void update(int *ex);
	virtual double getprob(int *ex);
};

class Stats4 : public Stats {
public:
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;
	int c,m;
	int nvaluesc;
	int d,n;
	int nvaluesd;

	Stats4(int _a, int _b, int _c, int _d, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd);
	virtual void update(int *ex);
	virtual double getprob(int *ex);
};


class Stats5 : public Stats{
public:
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;
	int c,m;
	int nvaluesc;
	int d,n;
	int nvaluesd;
	int e,o;
	int nvaluese;

	Stats5(int _a, int _b, int _c, int _d, int _e, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd, int _nvaluese);
	virtual void update(int *ex);
	virtual double getprob(int *ex);
};


struct KStat {
	int n;
	int idx;
	Stats *s;
	double *pred;
};

struct KRegion {
	int magnitude;
	KStat *ks;
};



struct KCRegion {
	int n;
	Stats *train, *test;
};

struct KModelCompact {
	int n;
	int *indices;
	int *magnitudes;
};

struct KEnsemble {
	int n;

	vector<KRegion *> *regions;
	double **predictions;
	vector<double> *weights;
	vector<KModelCompact> *indices;
};

class KCache {
	vector<KRegion *> *regions;
	KEnsemble *ensemble;
	vector<KStat *> *region_bank;
	double *predictions, *saved; // for all instances
	double *temp, *best, *results; // for one instance
	double prior;
	int pdof, padof;
	map<long,int> lookup;
	map<long,int>::iterator found;

	KStat *makeaStat(int *p, int n);
	KStat *findaStat(int *p, int n);
public:
	KInput *d;
	int dof, adof;
	double loss;

	void addModel(struct KModel *m);
	void testModel(struct KModel *m, double *result, double *dof);
	void CheckReversal(double *outresult);
	double getLoss(int *ex);
	void emptyModel();
	void deLogize(double *outresult);
	void Classify(int *ex, double *outresult);
	void ClassifyEnsemble(int *ex, double *outresult);
	void setEnsemble(struct KModels *m, struct KArray *weights);
	void emptyEnsemble();
	~KCache();
	KCache(struct KInput *input, double prior);
};

// stores the model and learns the statistics
void Ksetmodel(struct KInfo *in, struct KModel *m);
void Kaddmodel(struct KInfo *in, struct KModel *m);
void Ktestaddition(struct KInfo *in, struct KModel *m, struct KList *OutValue);

void Kdie(struct KInfo *p);
void Klearn(struct KInfo *in, int samples, int depth, struct KMatrix *OutValue);
void Kcleanup(struct KInput *p);

void Ktestmodels(struct KInfo *in, struct KModels *ms, int samples, struct KMatrix *OutValue);
void Kuse(struct KInfo *in, int *ex, struct KList *OutValue);
double Kvalidate(struct KInfo *in, int *ex);
void KgetDOF(struct KInfo *in, struct KList *OutValue);
void Kcheckreversal(struct KInfo *in, struct KList *OutValue);
struct KInfo *Kremember(struct KInput *input, double prior);
void Ksetensemble(struct KInfo *in, struct KModels *ms, struct KArray *weights);
void Kuseensemble(struct KInfo *in, int *ex, struct KList *OutValue);

#endif