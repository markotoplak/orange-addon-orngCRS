#ifndef _MODELS_H_
#define _MODELS_H_

#include <vector>
#include <map>
#include <math.h>
#include "kikuchi.h"

using namespace std;

#define safelog(prob) ((prob >= 1e-8)?(log(prob)):(-18.5))


class Stats {
public:
	vector< int > *v; // frequencies
	int total;
	long allcomb;
	double divisor;
	double prior;

	void model();
	void reset();
	void set_prior(double _p);
	virtual void update(int *ex) = 0;
	virtual double getprob(int *ex ) = 0;
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


struct KRegion {
	int magnitude;
	int n;
	Stats *s;
};

struct KCRegion {
	int n;
	Stats *train, *test;
};

struct KCModel {
	int n;
	int na;
	int *indices;
	int *magnitudes;
	int *stats;
};


double KirkwoodClass(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex);
double KirkwoodInt(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex);
double KirkwoodJoint(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex);

double KikuchiClass(const int *cards, KCModel *m, vector<KCRegion *> &regions, double *probs, double *tprobs, double *oddss, double *toddss, int *ex);
double KikuchiInt(const int *cards, KCModel *m, vector<KCRegion *> &regions, double *probs, double *tprobs, double *oddss, double *toddss, int *ex);

#endif