#ifndef _MODELS_H_
#define _MODELS_H_

#include <vector>
#include <map>
#include <math.h>
#include "kikuchi.h"

using namespace std;



double KirkwoodClass(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex);
double KirkwoodInt(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex);
double KirkwoodJoint(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex);

double KikuchiClass(const int *cards, KCModel *m, vector<KCRegion *> &regions, double *probs, double *tprobs, double *oddss, double *toddss, int *ex);
double KikuchiInt(const int *cards, KCModel *m, vector<KCRegion *> &regions, double *probs, double *tprobs, double *oddss, double *toddss, int *ex);

#endif