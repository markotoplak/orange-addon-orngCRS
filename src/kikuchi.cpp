#include "kikuchi.h"
#include "models.h"
#include "../statistics/pr_statistics.h"
#include "cMersenneTwister.h"
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#define MAX_INTS (100000)




#define DELETEVECTOR(vv) {for (int i = 0; i < vv->size(); ++i) delete (*vv)[i]; delete vv;}
#ifndef minx
template <class T> inline T minx(T x,T y) { return (x<y)?x:y; }
#endif
#ifndef maxx
template <class T> inline T maxx(T x,T y) { return (x>y)?x:y; }
#endif


Stats::Stats(long _allcomb) {
	allcomb = _allcomb;
	v = new vector< int >(allcomb,0);
	total = 0;
	prior = 0.0;
}

void Stats::set_prior(double _p) {
	prior = _p;
}

void Stats::model() {
	if(total+prior > 0) 
		divisor = 1.0/(total+allcomb*prior);
	else
		divisor = 0.0;
}

void Stats::reset() {
	for (int i = 0; i < allcomb; ++i) {
		(*v)[i] = 0;
	}
	total = 0;
}

Stats::~Stats() {
	delete v;
}


Stats1::Stats1(int _a, int _nvalues) : Stats(_nvalues) {
	nvaluesa = _nvalues; 
	a = _a; // attribute identifier
};

void Stats1::update(int *ex) {
	(*(v))[ex[a]]++;
	total++;
}

double Stats1::getprob(int *ex) {
	return (prior + (*(v))[ex[a]])*divisor; 
}

Stats2::Stats2(int _a, int _b, int _nvaluesa, int _nvaluesb) : Stats(_nvaluesa*_nvaluesb) {
	a = _a; // attribute identifier
	b = _b;
	nvaluesa = _nvaluesa; 
	nvaluesb = _nvaluesb; 
};

void Stats2::update(int *ex) {
	(*v)[ex[a]*nvaluesb + ex[b]]++;
	total++;
}

double Stats2::getprob(int *ex) {
	int x = ex[a]*nvaluesb + ex[b];
	return (prior + (*(v))[x])*divisor; 
}

Stats3::Stats3(int _a, int _b, int _c, int _nvaluesa, int _nvaluesb, int _nvaluesc) : Stats(_nvaluesa*_nvaluesb*_nvaluesc) {
	a = _a; // attribute identifier
	b = _b;
	c = _c;
	nvaluesa = _nvaluesa; 
	nvaluesb = _nvaluesb; 
	nvaluesc = _nvaluesc; 
	m = nvaluesb*nvaluesc;
};

void Stats3::update(int *ex) {
	(*v)[ex[a]*m + ex[b]*nvaluesc + ex[c]]++;
	total++;
}

double Stats3::getprob(int *ex) {
	int x = ex[a]*m + ex[b]*nvaluesc + ex[c];
	return (double)(prior+(*(v))[x])*divisor; 
}

Stats4::Stats4(int _a, int _b, int _c, int _d, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd) :
		Stats(_nvaluesa*_nvaluesb*_nvaluesc*_nvaluesd)
{
	a = _a; // attribute identifier
	b = _b;
	c = _c;
	d = _d;
	nvaluesa = _nvaluesa; 
	nvaluesb = _nvaluesb; 
	nvaluesc = _nvaluesc; 
	nvaluesd = _nvaluesd; 
	m = nvaluesd*nvaluesc;
	n = m*nvaluesb;
};

void Stats4::update(int *ex) {
	(*v)[ex[a]*n + ex[b]*m + ex[c]*nvaluesd + ex[d]]++;
	total++;
};

double Stats4::getprob(int *ex) {
	int x = ex[a]*n + ex[b]*m + ex[c]*nvaluesd + ex[d];
	return (double)(prior+(*(v))[x])*divisor; 
};

Stats5::Stats5(int _a, int _b, int _c, int _d, int _e, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd, int _nvaluese) :
		Stats(_nvaluesa*_nvaluesb*_nvaluesc*_nvaluesd*_nvaluese) 
{
	a = _a; // attribute identifier
	b = _b;
	c = _c;
	d = _d;
	e = _e;
	nvaluesa = _nvaluesa; 
	nvaluesb = _nvaluesb; 
	nvaluesc = _nvaluesc; 
	nvaluesd = _nvaluesd; 
	nvaluese = _nvaluese; 
	m = nvaluesd*nvaluese;
	n = m*nvaluesc;
	o = n*nvaluesb;
};

void Stats5::update(int *ex) {
	(*v)[ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + ex[e]]++;
	total++;
}

double Stats5::getprob(int *ex) {
	int x = ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + ex[e];
	return (double)(prior+(*(v))[x])*divisor; 
}

class KCache {
	vector<KRegion *> *regions;
	
	void prepareModel(struct KModel *m, double prior) {
		int i,j;
		
		for(i = 0; i < m->nogroups; ++i) {
			int *p = m->groups[i].l;
			KRegion *r = new KRegion;
			r->magnitude = m->groups[i].times;
			r->n = m->groups[i].n;
			switch(r->n) {
				case 1:
					r->s = (Stats*)new Stats1(p[0],d->card[p[0]]);
					break;
				case 2:
					r->s = (Stats*)new Stats2(p[0],p[1],d->card[p[0]],d->card[p[1]]);
					break;
				case 3:
					r->s = (Stats*)new Stats3(p[0],p[1],p[2],d->card[p[0]],d->card[p[1]],d->card[p[2]]);
					break;
				case 4:
					r->s = (Stats*)new Stats4(p[0],p[1],p[2],p[3],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]]);
					break;
				case 5:
					r->s = (Stats*)new Stats5(p[0],p[1],p[2],p[3],p[4],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]],d->card[p[4]]);
					break;
				default:
					fprintf(stderr, "Region of size %d is inappropriate.\n", m->groups[i].n);
					return;
			}
			// enter the data
			for (j = 0; j < d->nn; ++j) {
				r->s->update(d->data[j]);
			}
			r->s->set_prior(prior);
			r->s->model();
			regions->push_back(r);
		}
	}

public:
	KInput *d;

	void Test(int depth, int iterations) {
		// prepare the stats

	}

	double *Classify(int *ex) {
		int i, maxi;
		unsigned int j;
		double odds, prob, lprob, sum, norm, maxo, mino;
		double *results, *best; // results buffer

		results = (double *)malloc(sizeof(double)*d->card[d->na]);
		best = (double *)malloc(sizeof(double)*d->card[d->na]);

		// for all the classes
		sum = 0.0;
		for(i = 0; i < d->card[d->na]; ++i){
			ex[d->na] = i;
			odds = 0.0;
			for(j = 0; j < regions->size(); ++j) {
				prob = (*regions)[j]->s->getprob(ex);
				if(prob >= 1e-8)
					lprob = log(prob);
				else
					lprob = -18.5;
				odds += (*regions)[j]->magnitude*lprob;
			}
			results[i] = odds;
			sum += exp(results[i]);
		}

		norm = log(sum);
		maxi = 0;
		maxo = results[i];
		mino = 0.0;
		sum = 0.0;
		best[0] = 1.0;
		for(i = 0; i < d->card[d->na]; ++i) {
			if(results[i] > maxo) {
				for(j = 0; j < i; ++j)
					best[j] = 0.0;
				best[i] = 1.0;
				maxo = results[i];
				maxi = i;
			} else if(results[i] == maxo) {
				best[i] = 1.0;
			} else {
				best[i] = 0.0;
			}
			results[i] = exp(results[i]-norm);
			mino = min(mino,results[i]);
			sum += results[i];
		}
		sum += mino*d->card[d->na]; // correct for negative minimum
		if(sum < 1e-200) {
			// everything failed, send the backup
			sum = 0.0;
			for(i = 0; i < d->card[d->na]; ++i)
				sum += best[i];
			sum = 1.0/sum;
			for(i = 0; i < d->card[d->na]; ++i)
				best[i] *= sum;
			free(results);
			return best;
		} else {
			sum = 1.0/sum;
			for(i = 0; i < d->card[d->na]; ++i) {
				results[i] *= sum;
			}
			free(best);
			return results;
		}
	}

	~KCache() {
		for (unsigned int i = 0; i < regions->size(); ++i) 
			delete (*regions)[i]->s; 
		delete regions;
	}

	KCache(struct KInput *input, struct KModel *m, double prior) {
		d = input;
		regions = new vector<KRegion *>;
		this->prepareModel(m,prior);
	}
};


void Kcleanup(struct KInput *p) {
	int i;
	if (p->data != NULL) {
		for (i=0; i < p->nn; ++i)
			free(p->data[i]);
		free(p->data);
	}
	if (p->card != NULL) {
		free(p->card);
	}
	free(p);
}



// analyzes the data
void Klearn(struct KInfo *in, int samples, int depth, struct KMatrix *OutValue) {
	KCache *c;
	int s, d, columns = 0, i, totalmemory, prev, curr;
	int x1,x2,x3,x4, n, m, *ex, *tex, *indices, xoffs;
	vector< int * > *LUT = new vector<int *>(depth+1);
	vector< int > *idx; // bootstrap indices
	vector< Stats * > *stat_train = new vector <Stats *>;
	vector< Stats * > *stat_test  = new vector <Stats *>;
	double *probs, *oddss, *toddss, *tprobs, log2fix;

	c = (KCache *)(in)->c;
	log2fix = 1.0/log(2.0);
	idx = new vector<int>(c->d->nn);
	tex = new int[c->d->k]; // temporary example
	probs = new double[c->d->card[c->d->na]]; // temporary class probabilities
	oddss = new double[c->d->card[c->d->na]]; // temporary class odds
	toddss = new double[c->d->card[c->d->na]]; // true temporary class odds
	tprobs = new double[c->d->card[c->d->na]]; // true temporary class probs
	indices = new int[depth+1];

	cMersenneTwister *rng = new cMersenneTwister((long)rand()*(long)rand());

	// compute the number of combinations
	stat_train->push_back((Stats *) new Stats1(c->d->na,c->d->card[c->d->na]));
	stat_test->push_back((Stats *) new Stats1(c->d->na,c->d->card[c->d->na]));
	columns = 1;
	totalmemory = 0;
	i = c->d->na;
	prev = 1;
	for(d = 1; d <= depth; ++d) {
		int *pi;

		indices[d] = -1;

		curr = prev * i;
		prev = curr;
		totalmemory += curr;
		if (totalmemory > MAX_INTS) {
			// limit depth and prevent failure due to excessive number of combinations
			depth = d-1; 
			break;
		}
		pi = (*LUT)[d] = new int[curr];

		// combination loop
		for(x1 = 0; x1 < i; ++x1) {
			if(d == 1) {
				pi[x1] = columns++;
				stat_train->push_back((Stats *) new Stats2(i,x1,c->d->card[i],c->d->card[x1]));
				stat_test->push_back((Stats *) new Stats2(i,x1,c->d->card[i],c->d->card[x1]));
			} else 
				for(x2 = x1+1; x2 < i; ++x2) {
					if (d==2) {
						pi[x1+x2*i] = columns++;
						stat_train->push_back((Stats *) new Stats3(i,x1,x2,c->d->card[i],c->d->card[x1],c->d->card[x2]));
						stat_test->push_back((Stats *) new Stats3(i,x1,x2,c->d->card[i],c->d->card[x1],c->d->card[x2]));
					} else
						for(x3 = x2+1; x3 < i; ++x3) {
							if(d == 3) {
								pi[x1+(x2+x3*i)*i] = columns++;
								stat_train->push_back((Stats *) new Stats4(i,x1,x2,x3,c->d->card[i],c->d->card[x1],c->d->card[x2],c->d->card[x3]));
								stat_test->push_back((Stats *) new Stats4(i,x1,x2,x3,c->d->card[i],c->d->card[x1],c->d->card[x2],c->d->card[x3]));
							} else
								for(x4 = x3+1; x4 < i; ++x4) {
									if(d == 4) {
										pi[x1+(x2+(x3+x4*i)*i)*i] = columns++;
										stat_train->push_back((Stats *) new Stats5(i,x1,x2,x3,x4,c->d->card[i],c->d->card[x1],c->d->card[x2],c->d->card[x3],c->d->card[x4]));
										stat_test->push_back((Stats *) new Stats5(i,x1,x2,x3,x4,c->d->card[i],c->d->card[x1],c->d->card[x2],c->d->card[x3],c->d->card[x4]));
									} else
										assert(1==0);
								}
						}
				}
		}
	}

	// allocation
	OutValue->columns = columns;
	OutValue->rows = samples;
	OutValue->l = (double *)malloc(sizeof(double)*(OutValue->columns)*(OutValue->rows));

	// reference training
	for(n = 0; n < c->d->nn; ++n) {
		for(m = 0; m < columns; ++m) {
			(*stat_train)[m]->update(c->d->data[n]);
		}
	}
	for(m = 0; m < columns; ++m) {
		(*stat_train)[m]->set_prior(1.0);
		(*stat_train)[m]->model();
	}

	xoffs = 0;
	for(s = 0; s < samples; ++s) {
		// create the bootstrap resamples
		// reset the counts
		for(m = 0; m < columns; ++m) {
			(*stat_test)[m]->reset();
		}
		if(s > 0) {
			for(n = 0; n < c->d->nn; ++n) {
				//ex = c->d->data[rand()%c->d->nn];
				ex = c->d->data[rng->Random()%c->d->nn];
				for(m = 0; m < columns; ++m) {
					(*stat_test)[m]->update(ex);
				}
			}
		} else {
			// the first time without randomization
			for(n = 0; n < c->d->nn; ++n) {
				ex = c->d->data[n];
				for(m = 0; m < columns; ++m) {
					(*stat_test)[m]->update(ex);
				}
			}
		}
		// model the probabilities
		for(m = 0; m < columns; ++m) {
			(*stat_test)[m]->model();
		}
		if (s > 0) {
			// test all interaction models
			for(d = 0; d <= depth; ++d) {	
				if(d == 0) {
					indices[0]=i; // don't corrupt the data!
					OutValue->l[xoffs++] = log2fix*KirkwoodInt(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
				} else
					for(x1 = 0; x1 < i; ++x1) {
						indices[0] = x1;
						if(d == 1) {
							OutValue->l[xoffs++] = log2fix*KirkwoodInt(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
						} else 
							for(x2 = x1+1; x2 < i; ++x2) {
								indices[1] = x2;
								if(d == 2) {
									OutValue->l[xoffs++] = log2fix*KirkwoodInt(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
								} else 
									for(x3 = x2+1; x3 < i; ++x3) {
										indices[2] = x3;
										if(d == 3) {
											OutValue->l[xoffs++] = log2fix*KirkwoodInt(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
										} else 
											for(x4 = x3+1; x4 < i; ++x4) {
												indices[3] = x4;
												assert(d==4);
												OutValue->l[xoffs++] = log2fix*KirkwoodInt(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
											}
									}
							}
					}
			}
		} else {
			// test all part-to-whole models
			for(d = 0; d <= depth; ++d) {	
				if(d == 0) {
					indices[0]=i; // don't corrupt the data!
					OutValue->l[xoffs++] = log2fix*KirkwoodClass(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
				} else
					for(x1 = 0; x1 < i; ++x1) {
						indices[0] = x1;
						if(d == 1) {
							OutValue->l[xoffs++] = log2fix*KirkwoodClass(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
						} else 
							for(x2 = x1+1; x2 < i; ++x2) {
								indices[1] = x2;
								if(d == 2) {
									OutValue->l[xoffs++] = log2fix*KirkwoodClass(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
								} else 
									for(x3 = x2+1; x3 < i; ++x3) {
										indices[2] = x3;
										if(d == 3) {
											OutValue->l[xoffs++] = log2fix*KirkwoodClass(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
										} else 
											for(x4 = x3+1; x4 < i; ++x4) {
												indices[3] = x4;
												assert(d==4);
												OutValue->l[xoffs++] = log2fix*KirkwoodClass(d, i, indices, c->d->card, probs, tprobs, oddss, toddss, stat_train, stat_test, *LUT, tex);
											}
									}
							}
					}
			}
		}
	}

	// purge memory
	for(d = 1; d <= depth; ++d) {
		delete[] (*LUT)[d];
	}
	assert(columns == stat_train->size() && columns == stat_test->size());
	for(d = 0; d < columns; ++d) {
		delete (*stat_train)[d];
		delete (*stat_test)[d];
	}
	delete LUT, stat_train, stat_test, idx;
	delete[] tex;
	delete[] probs;
	delete[] oddss;
	delete[] toddss;
	delete[] tprobs;
	delete[] indices;
	delete rng;
}

void Ktestmodels(struct KInfo *in, struct KModels *ms, int samples, struct KMatrix *OutValue) {
	KCache *c;
	int s, columns = 0, i,j,k;
	int n, *ex, *tex, xoffs;
	vector< int > *idx; // bootstrap indices
	KCRegion *r;
	KCModel *mod;
	KInput *d;
	double *probs, *oddss, *toddss, *tprobs, log2fix;
	long N;

	c = (KCache *)(in)->c;
	if (c == NULL) {
		// no cache yet... prepare
		KModel tm;
		tm.groups = NULL;
		tm.nogroups = 0;
		in->c = new KCache(in->i,&tm,0.0);
		c = (KCache *)in->c;
	}
	d = c->d;
	log2fix = 1.0/log(2.0);
	idx = new vector<int>(c->d->nn);
	tex = new int[c->d->k]; // temporary example
	probs = new double[c->d->card[c->d->na]+10]; // temporary class probabilities
	oddss = new double[c->d->card[c->d->na]+10]; // temporary class odds
	toddss = new double[c->d->card[c->d->na]+10]; // true temporary class odds
	tprobs = new double[c->d->card[c->d->na]+10]; // true temporary class probs

	cMersenneTwister *rng = new cMersenneTwister((long)rand()*(long)rand());
	vector<KCRegion *> regions;
	vector<KCModel *> models; 

	map<long,int> lookup;
	map<long,int>::iterator found;
	N = c->d->k; // number of information sources
	for(k = 0; k < ms->nomodels; ++k) {
		KModel *m = ms->models + k;
		// COPY THE MODEL
		mod = new KCModel;
		mod->n = m->nogroups;
		mod->magnitudes = (int *)malloc(sizeof(int)*mod->n);
		mod->stats      = (int *)malloc(sizeof(int)*mod->n);
		// copy attribute indices
		mod->na = m->groups[0].n;
		mod->indices    = (int *)malloc(sizeof(int)*mod->na);
		for(i = 0; i < mod->na; ++i) {
			mod->indices[i] = m->groups[0].l[i];
		}
		assert(mod->indices[mod->na-1] == c->d->na); // last must be the label
		models.push_back(mod);

		for(i = 0; i < m->nogroups; ++i) {
			long key;
			int idx;
			int *p = m->groups[i].l;
			r = new KCRegion;
			r->n = m->groups[i].n;
			mod->magnitudes[i] = m->groups[i].times;

			// BUILD THE KEY
			key = 0; j = 0;
			do {
				assert(key < (2147483647-N)/N); // prevent overflow
				key *= N;
				key += 1+p[j];
			} while (++j < r->n);
			found = lookup.find(key);
			if (found == lookup.end()) {
				// not found
				switch(r->n) {
					case 1:
						r->train = (Stats*)new Stats1(p[0],d->card[p[0]]);
						r->test = (Stats*)new Stats1(p[0],d->card[p[0]]);
						break;
					case 2:
						r->train = (Stats*)new Stats2(p[0],p[1],d->card[p[0]],d->card[p[1]]);
						r->test = (Stats*)new Stats2(p[0],p[1],d->card[p[0]],d->card[p[1]]);
						break;
					case 3:
						r->train = (Stats*)new Stats3(p[0],p[1],p[2],d->card[p[0]],d->card[p[1]],d->card[p[2]]);
						r->test = (Stats*)new Stats3(p[0],p[1],p[2],d->card[p[0]],d->card[p[1]],d->card[p[2]]);
						break;
					case 4:
						r->train = (Stats*)new Stats4(p[0],p[1],p[2],p[3],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]]);
						r->test = (Stats*)new Stats4(p[0],p[1],p[2],p[3],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]]);
						break;
					case 5:
						r->train = (Stats*)new Stats5(p[0],p[1],p[2],p[3],p[4],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]],d->card[p[4]]);
						r->test = (Stats*)new Stats5(p[0],p[1],p[2],p[3],p[4],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]],d->card[p[4]]);
						break;
					default:
						fprintf(stderr, "Region of size %d is inappropriate.\n", m->groups[i].n);
						return;
				}
				idx = regions.size();
				lookup[key] = idx;
				regions.push_back(r);
			} else {
				// found
				idx = found->second;
				r->train = regions[idx]->train;
				r->test = regions[idx]->test;
			}
			mod->stats[i] = idx;
		}
	}

	// allocation
	OutValue->columns = ms->nomodels;
	OutValue->rows = samples;
	OutValue->l = (double *)malloc(sizeof(double)*(OutValue->columns)*(OutValue->rows));

	// reference training
	for(n = 0; n < c->d->nn; ++n) {
		// iterate over all keys
		for(found = lookup.begin(); found != lookup.end(); ++found) {
			regions[found->second]->train->update(c->d->data[n]);
			regions[found->second]->test->update(c->d->data[n]);
		}
	}
	for(found = lookup.begin(); found != lookup.end(); ++found) {
		regions[found->second]->train->set_prior(1.0);
		regions[found->second]->train->model();
		regions[found->second]->test->model();
	}

	xoffs = 0;
	for(s = 0; s < samples; ++s) {
		// create the bootstrap resamples
		if(s > 0) {
			// reset the counts
			for(found = lookup.begin(); found != lookup.end(); ++found) {
				regions[found->second]->test->reset();
			}
			for(n = 0; n < c->d->nn; ++n) {
				//ex = c->d->data[rand()%c->d->nn];
				ex = c->d->data[rng->Random()%c->d->nn];
				for(found = lookup.begin(); found != lookup.end(); ++found) {
					regions[found->second]->test->update(ex);
				}
			}
			// model the probabilities
			for(found = lookup.begin(); found != lookup.end(); ++found) {
				regions[found->second]->test->model();
			}
		}
		if (s == 0) {
			// compute divergence
			for(i = 0; i < models.size(); ++i) {
				OutValue->l[xoffs++] = log2fix*KikuchiClass(c->d->card, models[i], regions, probs, tprobs, oddss, toddss, tex);
			}
		} else {
			// compute the self-divergence
			for(i = 0; i < models.size(); ++i) {
				OutValue->l[xoffs++] = log2fix*KikuchiInt(c->d->card, models[i], regions, probs, tprobs, oddss, toddss, tex);
			}
		}
	}

	for(found = lookup.begin(); found != lookup.end(); ++found) {
		delete regions[found->second]->train;
		delete regions[found->second]->test;
	}
	for(i = 0; i < models.size(); ++i) {
		free(models[i]->indices);
		free(models[i]->magnitudes);
		free(models[i]->stats);
	}
	delete[] tex;
	delete[] probs;
	delete[] oddss;
	delete[] toddss;
	delete[] tprobs;
	delete rng;
}


// stores the model and learns the statistics
void Kprepare(struct KInfo *in, double prior, struct KModel *m) {
	in->c = (void *)new KCache(in->i,m, prior);
}

void Kdie(struct KInfo *p) {
	Kcleanup(p->i);
	delete ((KCache *)(p)->c);
	free(p);
}

void Kuse(struct KInfo *in, int *ex, struct KList *OutValue) {
	OutValue->n = ((KCache *)(in)->c)->d->card[((KCache *)(in)->c)->d->na];
	OutValue->l = ((KCache *)(in)->c)->Classify(ex);
}

// creates the data structure
struct KInfo *Kremember(struct KInput *input) {
	struct KInfo *i;

	i = (struct KInfo *)malloc(sizeof(struct KInfo));

	// gather the total counts
	i->i = input;
	i->c = NULL;

	return i;
}