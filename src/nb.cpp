#include "nb.h"

#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#define DELETEVECTOR(vv) {for (int i = 0; i < vv->size(); ++i) delete (*vv)[i]; delete vv;}
#ifndef minx
template <class T> inline T minx(T x,T y) { return (x<y)?x:y; }
#endif
#ifndef maxx
template <class T> inline T maxx(T x,T y) { return (x>y)?x:y; }
#endif

extern double alnorm(double x, bool upper);

double WilcoxonPaired(const vector<double> &a, const vector<double> &b, double &zz) {
	int i,p,n,npos,nprev,prevp,prevn;
	double bneg,bpos,prev,tiestore,z;
	vector<double> dd(a.size());
	vector<double> absdd(a.size());

	for(i = 0; i < a.size(); ++i) {
		dd[i] = b[i]-a[i];
	}
	dd.push_back(1e202); // guards
	dd.push_back(-1e202);
	sort(dd.begin(),dd.end());
	npos = 1;
	n = dd.size()-1;
	p = 0;
	prev = -1e200;
	nprev = 0;
	prevp = 0;
	prevn = 0;
	tiestore = 0.0;
	bneg = 0.0;
	bpos = 0.0;
	
	// get to the middle
	while(dd[p] < 1e-6)
		++p;
	while(dd[n] > -1e-6)
		--n;

	while(p < dd.size()-1 || n > 0) {
		// check for duplication
		if(fabs(dd[p]) - 1e-6 < prev) {
			prevp++;
			nprev++;
			++p;
		} else 
		if(fabs(dd[n]) - 1e-6 < prev) {
			prevn++;
			nprev++;
			--n;
		} else {
			if(nprev > 1) {
				// been a tie
				bneg += prevn*(npos+(nprev-1)/2.0);
				bpos += prevp*(npos+(nprev-1)/2.0);
				tiestore += (double)nprev*nprev*nprev-nprev;
				npos += nprev;
				prevn = prevp = nprev = 0;
			} else if (nprev == 1) {
				// normal diff
				bneg += prevn*npos;
				bpos += prevp*npos;
				prevp = prevn = 0;
				npos++;
				nprev = 0;
			}
			if(fabs(dd[p]) <= fabs(dd[n])) {
				// positive wins
				prev = fabs(dd[p]);
				prevp++;
				nprev++;
				++p;
			} else {
				// negative wins
				prev = fabs(dd[n]);
				prevn++;
				nprev++;
				--n;
			}
		}
	}
	if(nprev > 1) {
		// been a tie
		bneg += prevn*(npos+(nprev-1)/2.0);
		bpos += prevp*(npos+(nprev-1)/2.0);
		tiestore += (double)nprev*nprev*nprev-nprev;
		npos += nprev-1;
	} else if (nprev == 1) {
		// normal diff
		bneg += prevn*npos;
		bpos += prevp*npos;
	}
	//zz = z = (min(bneg,bpos) - npos*(npos+1.0)/4.0);
	zz = z = (bpos - npos*(npos+1.0)/4.0);
	z -= 0.5; // continuity correction: greater
	double t = sqrt(npos*(npos+1.0)*(2.0*npos+1.0)/24.0 + tiestore/48.0);
	if (t > 1e-6) // not to divide-by-zero
		z /= t;
	else
		z *= 1e6;
	if (bneg == 0 && bpos == 0) {
		return alnorm(0,true); // identical...
	} else
		return alnorm(z,true);
}




class NBStats {
public:
	virtual void update(int *ex) = 0;
	virtual double getprob(int *ex,int j) = 0;
	virtual double getprobLOO(int *ex,int j) = 0;
};


class NBStats1 : NBStats {
public:
	vector< double > *v; // frequencies
	double total;
	int a;
	int nvaluesa;

	NBStats1(int _a, int _nvalues, double ALPHA) {
		nvaluesa = _nvalues; 
		v = new vector< double >(nvaluesa,ALPHA);
		total = nvaluesa*ALPHA;
		a = _a; // attribute identifier
	};

	inline void update(int *ex) {
		(*(v))[ex[a]]++;
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		return (double)((*(v))[j])/(double)(total); 
	}

	inline double getprobLOO(int *ex,int j) {
		int x = j;

		// leave-one-out 
		if (ex[a] == j)
			return (double)((*(v))[x]-1)/(double)(total-1); 
		else
			return (double)((*(v))[x])/(double)(total-1); 
	}

	~NBStats1() {
		delete v;
	}
};

class NBStats2 : NBStats {
public:
	vector< double > *v; // frequencies
	double total;
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;

	NBStats2(int _a, int _b, int _nvaluesa, int _nvaluesb, double ALPHA) {
		a = _a; // attribute identifier
		b = _b;
		nvaluesa = _nvaluesa; 
		nvaluesb = _nvaluesb; 
		v = new vector< double >(nvaluesa*nvaluesb,ALPHA);
		total = nvaluesa*nvaluesb*ALPHA;
	};

	inline void update(int *ex) {
		(*v)[ex[a]*nvaluesb + ex[b]]++;
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		int x = ex[a]*nvaluesb + j;
		return (double)((*(v))[x])/(double)(total); 
	}

	inline double getprobLOO(int *ex, int j) {
		assert(total > 0);
		int x = ex[a]*nvaluesb + j;
		double p;
		//printf("%d %d %d %d\n",j,x,(*v)[x],total);

		if (ex[b] == j)
			p = (*v)[x]-1;
		else
			p = (*v)[x];

		return p/(total-1);
	}

	~NBStats2() {
		delete v;
	}
};

class NBStats3 : NBStats{
public:
	vector< double > *v;  // frequencies
	double total;
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;
	int c,m;
	int nvaluesc;

	NBStats3(int _a, int _b, int _c, int _nvaluesa, int _nvaluesb, int _nvaluesc, double ALPHA) {
		a = _a; // attribute identifier
		b = _b;
		c = _c;
		nvaluesa = _nvaluesa; 
		nvaluesb = _nvaluesb; 
		nvaluesc = _nvaluesc; 
		m = nvaluesb*nvaluesc;
		assert(m*nvaluesa < 256000);
		v = new vector< double >(nvaluesa*nvaluesb*nvaluesc,ALPHA);
		total = nvaluesa*nvaluesb*nvaluesc*ALPHA;
	};

	inline void update(int *ex) {
		(*v)[ex[a]*m + ex[b]*nvaluesc + ex[c]]++;
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		int x = ex[a]*m + ex[b]*nvaluesc + j;
		return (double)((*(v))[x])/(double)(total); 
	}

	inline double getprobLOO(int *ex,int j) {
		int x = ex[a]*m + ex[b]*nvaluesc + j;

		if (ex[c] == j)
			return (double)((*v)[x]-1)/(double)(total-1); 
		else
			return (double)((*v)[x])/(double)(total-1); 
	}

	~NBStats3(){
		delete v;
	}
};

class NBStats4 : NBStats{
public:
	vector< double > *v;  // frequencies
	double total;
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;
	int c,m;
	int nvaluesc;
	int d,n;
	int nvaluesd;

	NBStats4(int _a, int _b, int _c, int _d, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd, double ALPHA) {
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
		assert(n*nvaluesa < 256000);
		v = new vector< double >(nvaluesa*nvaluesb*nvaluesc*nvaluesd,ALPHA);
		total = nvaluesa*nvaluesb*nvaluesc*nvaluesd*ALPHA;
	};

	inline void update(int *ex) {
		(*v)[ex[a]*n + ex[b]*m + ex[c]*nvaluesd + ex[d]]++;
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		int x = ex[a]*n + ex[b]*m + ex[c]*nvaluesd + j;
		return (double)((*(v))[x])/(double)(total); 
	}

	inline double getprobLOO(int *ex,int j) {
		int x = ex[a]*n + ex[b]*m + ex[c]*nvaluesd + j;

		if (ex[d] == j)
			return (double)((*v)[x]-1)/(double)(total-1); 
		else
			return (double)((*v)[x])/(double)(total-1); 
	}

	~NBStats4(){
		delete v;
	}
};


class NBStats5 : NBStats{
public:
	vector< double > *v;  // frequencies
	double total;
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

	NBStats5(int _a, int _b, int _c, int _d, int _e, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd, int _nvaluese, double ALPHA) {
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
		assert(o*nvaluesa < 256000 && o*nvaluesa >= 0);
		v = new vector< double >(nvaluesa*nvaluesb*nvaluesc*nvaluesd*nvaluese,ALPHA);
		total = nvaluesa*nvaluesb*nvaluesc*nvaluesd*nvaluese*ALPHA;
	};

	inline void update(int *ex) {
		(*v)[ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + ex[e]]++;
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		int x = ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + j;
		return (double)((*(v))[x])/(double)(total); 
	}

	inline double getprobLOO(int *ex,int j) {
		int x = ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + j;

		if (ex[e] == j)
			return (double)((*v)[x]-1)/(double)(total-1); 
		else
			return (double)((*v)[x])/(double)(total-1); 
	}

	~NBStats5(){
		delete v;
	}
};

class TNBStats3 : NBStats{
public:
	vector< double > *v;  // frequencies
	double total;
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;
	int c,m;
	int nvaluesc;
	NBStats2 *sub;

	TNBStats3(int _a, int _b, int _c, int _nvaluesa, int _nvaluesb, int _nvaluesc, double ALPHA) {
		a = _a; // attribute identifier
		b = _b;
		c = _c;
		nvaluesa = _nvaluesa; 
		nvaluesb = _nvaluesb; 
		nvaluesc = _nvaluesc; 
		m = nvaluesb*nvaluesc;
		assert(m*nvaluesa < 256000);
		v = new vector< double >(nvaluesa*nvaluesb*nvaluesc,ALPHA);
		total = nvaluesa*nvaluesb*nvaluesc*ALPHA;
		sub = new NBStats2(b,c,nvaluesb,nvaluesc,ALPHA);
	};

	inline void update(int *ex) {
		(*v)[ex[a]*m + ex[b]*nvaluesc + ex[c]]++;
		sub->update(ex);
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		int x = ex[a]*m + ex[b]*nvaluesc + j;
		return (double)((*(v))[x])/(double)(total)/sub->getprob(ex,j); 
	}

	inline double getprobLOO(int *ex,int j) {
		int x = ex[a]*m + ex[b]*nvaluesc + j;

		if (ex[c] == j)
			return (double)((*v)[x]-1)/(double)(total-1)/sub->getprobLOO(ex,j); 
		else
			return (double)((*v)[x])/(double)(total-1)/sub->getprobLOO(ex,j); 
	}

	~TNBStats3(){
		delete v;
		delete sub;
	}
};

class TNBStats4 : NBStats{
public:
	vector< double  > *v;  // frequencies
	double  total;
	int a;
	int nvaluesa;
	int b;
	int nvaluesb;
	int c,m;
	int nvaluesc;
	int d,n;
	int nvaluesd;
	NBStats3 *sub;

	TNBStats4(int _a, int _b, int _c, int _d, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd, double ALPHA) {
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
		assert(n*nvaluesa < 256000);
		v = new vector< double  >(nvaluesa*nvaluesb*nvaluesc*nvaluesd,ALPHA);
		total = nvaluesa*nvaluesb*nvaluesc*nvaluesd*ALPHA;
		sub = new NBStats3(b,c,d,nvaluesb,nvaluesc,nvaluesd,ALPHA);
	};

	inline void update(int *ex) {
		(*v)[ex[a]*n + ex[b]*m + ex[c]*nvaluesd + ex[d]]++;
		sub->update(ex);
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		int x = ex[a]*n + ex[b]*m + ex[c]*nvaluesd + j;
		return (double)((*(v))[x])/(double)(total)/sub->getprob(ex,j); 
	}

	inline double getprobLOO(int *ex,int j) {
		int x = ex[a]*n + ex[b]*m + ex[c]*nvaluesd + j;

		if (ex[d] == j)
			return (double)((*v)[x]-1)/(double)(total-1)/sub->getprobLOO(ex,j); 
		else
			return (double)((*v)[x])/(double)(total-1)/sub->getprobLOO(ex,j); 
	}

	~TNBStats4(){
		delete v;
		delete sub;
	}
};


class TNBStats5 : NBStats{
public:
	vector< double  > *v;  // frequencies
	double  total;
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
	NBStats4 *sub;

	TNBStats5(int _a, int _b, int _c, int _d, int _e, int _nvaluesa, int _nvaluesb, int _nvaluesc, int _nvaluesd, int _nvaluese, double ALPHA) {
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
		assert(o*nvaluesa < 256000 && o*nvaluesa >= 0);
		v = new vector< double  >(nvaluesa*nvaluesb*nvaluesc*nvaluesd*nvaluese,ALPHA);
		total = nvaluesa*nvaluesb*nvaluesc*nvaluesd*nvaluese*ALPHA;
		sub = new NBStats4(b,c,d,e,nvaluesb,nvaluesc,nvaluesd,nvaluese,ALPHA);
	};

	inline void update(int *ex) {
		(*v)[ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + ex[e]]++;
		sub->update(ex);
		total++;
	}

	inline double getprob(int *ex,int j) {
		assert(total > 0);
		int x = ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + j;
		return (double)((*(v))[x])/(double)(total)/sub->getprob(ex,j); 
	}

	inline double getprobLOO(int *ex,int j) {
		int x = ex[a]*o + ex[b]*n + ex[c]*m + ex[d]*nvaluese + j;

		if (ex[e] == j)
			return (double)((*v)[x]-1)/(double)(total-1)/sub->getprobLOO(ex,j); 
		else
			return (double)((*v)[x])/(double)(total-1)/sub->getprobLOO(ex,j); 
	}

	~TNBStats5(){
		delete v;
		delete sub;
	}
};

class NBCache {
	vector<NBStats1 *>  *atts;
	vector<NBStats2 *> *cond;
	vector<NBStats *> *interact;
	vector<NBStats *> *tinteract;
	map<int,int>	 LUT;
	map<int,int>	 tLUT;
	vector<double> *kl_scores, *s_kl_s, *r_kl_s;
	vector<double> *brier_scores, *s_brier_s, *r_brier_s;
	vector<double> *ca_scores, *s_ca_s, *r_ca_s;
	vector<double> **probabilities;
	vector<int>		modeli;
	vector<int>		modela;
	vector<double>  weights;
	double k;
	
	void prepareCache(vector<int> &intl, struct NBModel *m, double k) {
		int i,j;
		
		this->k = k;
		intl.erase(intl.begin(),intl.end());

		for(i = 0; i < m->nogroups; ++i) {
			if (m->groups[i].n == 2) {
				NBStats3 *ss;

				int *p = m->groups[i].l;
				int val = minx(p[0],p[1])*d->k + maxx(p[0],p[1]);
				map<int,int>::iterator vv = LUT.find(val);
				if(vv == LUT.end()) {
					// the thing does not exist yet
					ss = new NBStats3(p[0],p[1],d->na,d->card[p[0]],d->card[p[1]],d->card[d->na],this->k);
					for (j = 0; j < d->nn; ++j) {
						ss->update(d->data[j]);
					}
					intl.push_back(interact->size());
					LUT[val] = interact->size();
					interact->push_back((NBStats *)ss);
				} else {
					intl.push_back(vv->second);
				}
			} else {
				if (m->groups[i].n > 2 && m->groups[i].n <= 4) {
					int key;
					
					int *p = m->groups[i].l;
					sort(p,p+m->groups[i].n);
					key = 0;
					for (j = 0; j < m->groups[i].n; ++j) {
						key *= d->k;
						key += m->groups[i].l[j];
					}
					map<int,int>::iterator vv = LUT.find(key);
					if(vv == LUT.end()) {
						// the thing does not exist yet
						if (m->groups[i].n == 3) {
							NBStats4 *ss;
							ss = new NBStats4(p[0],p[1],p[2],d->na,d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[d->na],this->k);
							for (j = 0; j < d->nn; ++j) {
								ss->update(d->data[j]);
							}
							intl.push_back(interact->size());
							LUT[key] = interact->size();
							interact->push_back((NBStats *)ss);
						} else if (m->groups[i].n == 4) {
							NBStats5 *ss;
							ss = new NBStats5(p[0],p[1],p[2],p[3],d->na,d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]],d->card[d->na],this->k);
							for (j = 0; j < d->nn; ++j) {
								ss->update(d->data[j]);
							}
							intl.push_back(interact->size());
							LUT[key] = interact->size();
							interact->push_back((NBStats *)ss);
						}
					} else {
						intl.push_back(vv->second);
					}
				}
			}
		}
	}

	void prepareTANCache(vector<int> &intl, struct NBModel *m, double k) {
		int i,j;
		this->k = k;
		
		intl.erase(intl.begin(),intl.end());

		for(i = 0; i < m->nogroups; ++i) {
			if (m->groups[i].n == 2) {
				TNBStats3 *ss;

				int *p = m->groups[i].l;
				int val = p[0]*d->k + p[1];
				map<int,int>::iterator vv = tLUT.find(val);
				if(vv == tLUT.end()) {
					// the thing does not exist yet
					ss = new TNBStats3(p[0],p[1],d->na,d->card[p[0]],d->card[p[1]],d->card[d->na],this->k);
					for (j = 0; j < d->nn; ++j) {
						ss->update(d->data[j]);
					}
					intl.push_back(tinteract->size());
					tLUT[val] = tinteract->size();
					tinteract->push_back((NBStats *)ss);
				} else {
					intl.push_back(vv->second);
				}
			} else {
				if (m->groups[i].n > 2 && m->groups[i].n <= 4) {
					int key;
					
					int *p = m->groups[i].l;
					key = 0;
					for (j = 0; j < m->groups[i].n; ++j) {
						key *= d->k;
						key += m->groups[i].l[j];
					}
					map<int,int>::iterator vv = tLUT.find(key);
					if(vv == tLUT.end()) {
						// the thing does not exist yet
						if (m->groups[i].n == 3) {
							TNBStats4 *ss;
							ss = new TNBStats4(p[0],p[1],p[2],d->na,d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[d->na],this->k);
							for (j = 0; j < d->nn; ++j) {
								ss->update(d->data[j]);
							}
							intl.push_back(tinteract->size());
							tLUT[key] = tinteract->size();
							tinteract->push_back((NBStats *)ss);
						} else if (m->groups[i].n == 4) {
							TNBStats5 *ss;
							ss = new TNBStats5(p[0],p[1],p[2],p[3],d->na,d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]],d->card[d->na],this->k);
							for (j = 0; j < d->nn; ++j) {
								ss->update(d->data[j]);
							}
							intl.push_back(tinteract->size());
							tLUT[key] = tinteract->size();
							tinteract->push_back((NBStats *)ss);
						}
					} else {
						intl.push_back(vv->second);
					}
				}
			}
		}
	}

public:
	NBInput *d;

	void Evaluate(struct NBModel *m, struct NBResult *v) {
		int i,j,k;
		vector<double> fs(d->card[d->na]);
		vector<int> intl;
		double q_brier, q_kl, prior, invp, f, sum, posterior;
		double stde_b, stde_kl, mean_b, mean_kl, stde_ca, mean_ca;

		// check if there are uncached attributes in the model
		prepareCache(intl,m,d->a);

		// perform the leave-one-out experiments
		mean_b = mean_kl = mean_ca = 0;
		for(i = 0; i < d->maxcvi; ++i) {
			int *ex = d->data[i];

			// for all the classes
			sum = 0.0;
			for(j = 0; j < d->card[d->na]; ++j) {
				prior = (*atts)[d->na]->getprobLOO(ex,j);
				f = prior;
				invp = 1.0/prior;
				// for all the attributes
				for(k = 0; k < m->nosingles; ++k) {
					f *= (*cond)[m->singles[k]]->getprobLOO(ex,j) * invp;
					//f /= (*cond)[m->singles[k]]->getprobLOO(ex,j);
				}
				// for all interactions
				for(k = 0; k < intl.size(); ++k) {
					f *= (*interact)[intl[k]]->getprobLOO(ex,j) * invp;
				}
				fs[j] = f;
				sum += f;
			}

			// compute the quality, given the true class
			// BRIER SCORE
			q_brier = 0.0;
			for(j = 0; j < d->card[d->na]; ++j) {
				f = fs[j]/sum;
				(*probabilities[j])[i] = f;
				//printf("%f ",f);
				q_brier += f*f;
			}
			//printf("\n");
			posterior = fs[ex[d->na]]/sum;
			q_brier = q_brier + 1 - 2*posterior;
			q_brier /= d->card[d->na];
			// KULLBACK-LEIBLER
			q_kl = -log(posterior);
			mean_kl += q_kl;
			mean_b += q_brier;
			(*kl_scores)[i] = q_kl;
			(*brier_scores)[i] = q_brier;

			f = fs[ex[d->na]];
			for(j = 0; j < d->card[d->na]; ++j) {
				if (f+1e-6 < fs[j]) {
					f = 0.0;
					break;
				}
			}
			if(f > 0.0) {
				(*ca_scores)[i] = 0.0;
			} else {
				(*ca_scores)[i] = 1.0;
			}
			mean_ca += (*ca_scores)[i];
		}
		// compute the standard error
		mean_kl /= d->maxcvi;
		mean_b /= d->maxcvi;
		mean_ca /= d->maxcvi;
		stde_kl = 0.0;
		stde_b = 0.0;
		stde_ca = 0.0;
		for(i = 0; i < d->maxcvi; ++i) {
			double t;
			
			t = (*kl_scores)[i]-mean_kl;
			stde_kl += t*t;
			t = (*brier_scores)[i]-mean_b;
			stde_b += t*t;
			t = (*ca_scores)[i]-mean_ca;
			stde_ca += t*t;
		}
		stde_kl /= d->maxcvi-1;
		stde_b /= d->maxcvi-1;
		stde_ca /= d->maxcvi-1;
		stde_kl = sqrt(stde_kl/d->maxcvi);
		stde_b = sqrt(stde_b/d->maxcvi);
		stde_ca = sqrt(stde_ca/d->maxcvi);


		v->kl_err = stde_kl; // KL
		v->kl_q = mean_kl;
		v->b_err = stde_b; // BRIER
		v->b_q = mean_b;
		v->er_err = stde_ca; // ERROR RATE
		v->er_q = mean_ca;
	}

	
	void TANEvaluate(struct NBModel *m, struct NBResult *v) {
		int i,j,k;
		vector<double> fs(d->card[d->na]);
		vector<int> intl;
		double q_brier, q_kl, prior, invp, f, sum, posterior;
		double stde_b, stde_kl, mean_b, mean_kl, stde_ca, mean_ca;

		// check if there are uncached attributes in the model
		prepareTANCache(intl,m,d->a);

		// perform the leave-one-out experiments
		mean_b = mean_kl = mean_ca = 0;
		for(i = 0; i < d->maxcvi; ++i) {
			int *ex = d->data[i];

			// for all the classes
			sum = 0.0;
			for(j = 0; j < d->card[d->na]; ++j) {
				prior = (*atts)[d->na]->getprobLOO(ex,j);
				f = prior;
				invp = 1.0/prior;
				// for all the attributes
				for(k = 0; k < m->nosingles; ++k) {
					f *= (*cond)[m->singles[k]]->getprobLOO(ex,j) * invp;
				}
				// for all interactions
				for(k = 0; k < intl.size(); ++k) {
					f *= (*tinteract)[intl[k]]->getprobLOO(ex,j);
				}
				fs[j] = f;
				sum += f;
			}

			// compute the quality, given the true class
			// BRIER SCORE
			q_brier = 0.0;
			for(j = 0; j < d->card[d->na]; ++j) {
				f = fs[j]/sum;
				(*probabilities[j])[i] = f;
				//printf("%f ",f);
				q_brier += f*f;
			}
			//printf("\n");
			posterior = fs[ex[d->na]]/sum;
			q_brier = q_brier + 1 - 2*posterior;
			q_brier /= d->card[d->na];
			// KULLBACK-LEIBLER
			q_kl = -log(posterior);
			mean_kl += q_kl;
			mean_b += q_brier;
			(*kl_scores)[i] = q_kl;
			(*brier_scores)[i] = q_brier;

			f = fs[ex[d->na]];
			for(j = 0; j < d->card[d->na]; ++j) {
				if (f+1e-6 < fs[j]) {
					f = 0.0;
					break;
				}
			}
			if(f > 0.0) {
				(*ca_scores)[i] = 0.0;
			} else {
				(*ca_scores)[i] = 1.0;
			}
			mean_ca += (*ca_scores)[i];
		}
		// compute the standard error
		mean_kl /= d->maxcvi;
		mean_b /= d->maxcvi;
		mean_ca /= d->maxcvi;
		stde_kl = 0.0;
		stde_b = 0.0;
		stde_ca = 0.0;
		for(i = 0; i < d->maxcvi; ++i) {
			double t;
			
			t = (*kl_scores)[i]-mean_kl;
			stde_kl += t*t;
			t = (*brier_scores)[i]-mean_b;
			stde_b += t*t;
			t = (*ca_scores)[i]-mean_ca;
			stde_ca += t*t;
		}
		stde_kl /= d->maxcvi-1;
		stde_b /= d->maxcvi-1;
		stde_ca /= d->maxcvi-1;
		stde_kl = sqrt(stde_kl/d->maxcvi);
		stde_b = sqrt(stde_b/d->maxcvi);
		stde_ca = sqrt(stde_ca/d->maxcvi);


		v->kl_err = stde_kl; // KL
		v->kl_q = mean_kl;
		v->b_err = stde_b; // BRIER
		v->b_q = mean_b;
		v->er_err = stde_ca; // ERROR RATE
		v->er_q = mean_ca;
	}

	
	void Divergence(struct NBModel *m, struct NBResult *v) {
		int i,j,k;
		vector<double> fs(d->card[d->na]);
		vector<double> as(d->card[d->na]);
		vector<int> intl;
		double q_brier, q_kl, prior, invp, f, sum, posterior, frac, tot;
		double stde_b, stde_kl, mean_b, mean_kl, stde_ca, mean_ca;
		int quit;
		int *ex;

		// check if there are uncached attributes in the model
		prepareCache(intl,m,d->a);

		//printf("%d %d\n",m->nosingles,intl.size());

		// a dummy example
		ex = (int *)calloc(d->na+1,sizeof(int));
		// a list of attributes appearing in the data
		vector<int> modelatts;
		for(k = 0; k < m->nosingles; ++k) {
			modelatts.push_back(m->singles[k]);
			ex[m->singles[k]] = 0;
		}
		for(k = 0; k < m->nogroups; ++k) {
			for (i = 0; i < (m->groups[k]).n; ++i) {
				j = (m->groups[k]).l[i];
				if(find(modelatts.begin(), modelatts.end(),j) == modelatts.end()) {
					modelatts.push_back(j);
					ex[j] = 0;
				}
			}
		}
		mean_b = mean_kl = mean_ca = 0;
		frac = 1.0/d->maxcvi;
		quit = 0;
		while(true) {
			// for all combinations of class values

			//for(j = 0; j <= d->na; ++j)
			//	printf("%d ",ex[j]);
			//printf("\n");

			tot = 0.0; // AS
			for(j = 0; j < d->card[d->na]; ++j)
				as[j] = 0.0;
			// FREQUENCY COUNT
			for(i = 0; i < d->maxcvi; ++i) {
				int *tex = d->data[i];
				int match = true;
				for(k = 0; k < modelatts.size(); ++k) {
					if(tex[modelatts[k]] != ex[modelatts[k]]) {
						match = false;
						break;
					}
				}
				if(match) {
					as[tex[d->na]] += 1.0;
					tot += 1.0;
				}
			}

			if(tot > 0.0) {
				sum = 0.0; // NB
				for(j = 0; j < d->card[d->na]; ++j) {
					// NAIVE BAYES
					prior = (*atts)[d->na]->getprob(ex,j);
					f = prior;
					invp = 1.0/prior;
					// for all the attributes
					for(k = 0; k < m->nosingles; ++k) {
						f *= (*cond)[m->singles[k]]->getprob(ex,j) * invp;
					}
					// for all interactions
					for(k = 0; k < intl.size(); ++k) {
						f *= (*interact)[intl[k]]->getprob(ex,j) * invp;
					}
					fs[j] = f;
					sum += f;
				}

				for(j = 0; j < d->card[d->na]; ++j) {
					double pa = frac*as[j];
					if(pa > 0) {
						// P(c|ab) -> p(a,b,c)/p(a,b) = f(abc)/f(ab)
						double ac = as[j]/tot;
						// P(c|ab) -> p'(c|a,b)
						double bc = fs[j]/sum;
						//printf("%d %f %f\n",j,ac,bc);

						mean_kl += pa*(log(ac)-log(bc));
						mean_b  += pa*(1.0/d->card[d->na])*(ac-bc)*(ac-bc);
					}
					
				}
			}

			while(quit < modelatts.size() && ex[modelatts[quit]] >= (d->card[modelatts[quit]]-1))
				++quit;
			if(quit < modelatts.size()) {
				ex[modelatts[quit]]++;
				for(k = 0; k < quit; ++k) {
					ex[modelatts[k]] = 0;
					quit = 0;
				}
			} else break;
		}

		v->kl_err = 0; // KL
		v->kl_q = mean_kl;
		v->b_err = 0; // BRIER
		v->b_q = mean_b;
		v->er_err = 0; // ERROR RATE
		v->er_q = mean_ca;

		free(ex);
	}


	void EvaluateW(struct NBModel *m, double *w, struct NBResult *v) {
		int i,j,k;
		vector<double> fs(d->card[d->na]);
		vector<int> intl;
		double q_brier, q_kl, prior, invp, f, sum, posterior;
		double stde_b, stde_kl, mean_b, mean_kl, stde_ca, mean_ca;

		// check if there are uncached attributes in the model
		prepareCache(intl, m,d->a);

		// perform the leave-one-out experiments
		mean_b = mean_kl = mean_ca = 0;
		for(i = 0; i < d->maxcvi; ++i) {
			int *ex = d->data[i];

			// for all the classes
			sum = 0.0;
			for(j = 0; j < d->card[d->na]; ++j) {
				int weight_c = 0;
				prior = (*atts)[d->na]->getprobLOO(ex,j);
				f = prior;
				invp = 1.0/prior;
				// for all the attributes
				for(k = 0; k < m->nosingles; ++k) {
					f *= pow((*cond)[m->singles[k]]->getprobLOO(ex,j) * invp,w[weight_c++]);
				}
				// for all interactions
				for(k = 0; k < intl.size(); ++k) {
					f *= pow((*interact)[intl[k]]->getprobLOO(ex,j) * invp,w[weight_c++]);
				}
				fs[j] = f;
				sum += f;
			}

			// compute the quality, given the true class
			// BRIER SCORE
			q_brier = 0.0;
			for(j = 0; j < d->card[d->na]; ++j) {
				f = fs[j]/sum;
				(*probabilities[j])[i] = f;
				//printf("%f ",f);
				q_brier += f*f;
			}
			//printf("\n");
			posterior = fs[ex[d->na]]/sum;
			q_brier = q_brier + 1 - 2*posterior;
			q_brier /= d->card[d->na];
			// KULLBACK-LEIBLER
			q_kl = -log(posterior);
			mean_kl += q_kl;
			mean_b += q_brier;
			(*kl_scores)[i] = q_kl;
			(*brier_scores)[i] = q_brier;

			f = fs[ex[d->na]];
			for(j = 0; j < d->card[d->na]; ++j) {
				if (f < fs[j]) {
					f = 0.0;
					break;
				}
			}
			if(f > 0.0) {
				(*ca_scores)[i] = 0.0;
			} else {
				(*ca_scores)[i] = 1.0;
			}
			mean_ca += (*ca_scores)[i];
		}
		// compute the standard error
		mean_kl /= d->maxcvi;
		mean_b /= d->maxcvi;
		mean_ca /= d->maxcvi;
		stde_kl = 0.0;
		stde_b = 0.0;
		stde_ca = 0.0;
		for(i = 0; i < d->maxcvi; ++i) {
			double t;
			
			t = (*kl_scores)[i]-mean_kl;
			stde_kl += t*t;
			t = (*brier_scores)[i]-mean_b;
			stde_b += t*t;
			t = (*ca_scores)[i]-mean_ca;
			stde_ca += t*t;
		}
		stde_kl /= d->maxcvi-1;
		stde_b /= d->maxcvi-1;
		stde_ca /= d->maxcvi-1;
		stde_kl = sqrt(stde_kl/d->maxcvi);
		stde_b = sqrt(stde_b/d->maxcvi);
		stde_ca = sqrt(stde_ca/d->maxcvi);


		v->kl_err = stde_kl; // KL
		v->kl_q = mean_kl;
		v->b_err = stde_b; // BRIER
		v->b_q = mean_b;
		v->er_err = stde_ca; // ERROR RATE
		v->er_q = mean_ca;
	}

	void StoreModel(struct NBModel *m, double *w) {
		int i,j;

		prepareCache(modeli, m,d->a);
		weights.erase(weights.begin(), weights.end());
		modela.erase(modela.begin(), modela.end());
		for(i = 0; i < m->nosingles; ++i) {
			weights.push_back(w[i]);
			modela.push_back(m->singles[i]);
		}
		for(j = 0; j < m->nogroups; ++j) {
			weights.push_back(w[i+j]);
		}
	}

	
	double *ClassifyW(int *ex) {
		int j,k;
		vector<int> intl;
		double prior, invp, f, sum;
		double *results; // results buffer

		results = (double *)malloc(sizeof(double)*d->card[d->na]);

		// for all the classes
		sum = 0.0;
		for(j = 0; j < d->card[d->na]; ++j) {
			int weight_c = 0;
			prior = (*atts)[d->na]->getprob(ex,j);
			f = prior;
			invp = 1.0/prior;
			// for all the attributes
			for(k = 0; k < modela.size(); ++k) {
				f *= pow((*cond)[modela[k]]->getprob(ex,j) * invp,weights[weight_c++]);
			}
			// for all interactions
			for(k = 0; k < modeli.size(); ++k) {
				f *= pow((*interact)[modeli[k]]->getprob(ex,j) * invp,weights[weight_c++]);
			}
			results[j] = f;
			sum += f;
		}
		if (sum < 1e-9) {
			for(j = 0; j < d->card[d->na]; ++j)
				// put in priors
				results[j] = (*atts)[d->na]->getprob(ex,j);
		} else
			for(j = 0; j < d->card[d->na]; ++j)
				results[j] /= sum;

		return results;
	}

	double *Classify(int *ex) {
		int j,k;
		vector<int> intl;
		double prior, invp, f, sum;
		double *results; // results buffer

		results = (double *)malloc(sizeof(double)*d->card[d->na]);

		// for all the classes
		sum = 0.0;
		for(j = 0; j < d->card[d->na]; ++j) {
			prior = (*atts)[d->na]->getprob(ex,j);
			f = prior;
			invp = 1.0/prior;
			// for all the attributes
			for(k = 0; k < modela.size(); ++k) {
				f *= (*cond)[modela[k]]->getprob(ex,j) * invp;
			}
			// for all interactions
			for(k = 0; k < modeli.size(); ++k) {
				f *= (*interact)[modeli[k]]->getprob(ex,j) * invp;
			}
			results[j] = f;
			sum += f;
		}

		//assert(sum > 1e-30);
		for(j = 0; j < d->card[d->na]; ++j)
			results[j] /= sum;

		return results;
	}

	void SaveScores() {
		for(int i = 0; i < d->nn; ++i) {
			(*s_kl_s)[i] = (*kl_scores)[i];
			(*s_brier_s)[i] = (*brier_scores)[i];
			(*s_ca_s)[i] = (*ca_scores)[i];
		}
	}

	void RememberScores() {
		for(int i = 0; i < d->nn; ++i) {
			(*r_kl_s)[i] = (*s_kl_s)[i];
			(*r_brier_s)[i] = (*s_brier_s)[i];
			(*r_ca_s)[i] = (*s_ca_s)[i];
		}
	}

	void ExportScores(int m, struct NBList *OutValue) {
		vector<double> *t;

		OutValue->n = d->nn;
		OutValue->l = (double *)malloc(sizeof(double)*d->nn);
		switch (m) {
		case 0: t = s_kl_s; break;
		case 1: t = s_brier_s; break;
		case 2: t = s_ca_s; break;
		}
		
		for(int i = 0; i < d->nn; ++i)
			OutValue->l[i] = (*t)[i];
	}

	void ExportProbabilities(int m, struct NBList *OutValue) {
		vector<double> *t;

		OutValue->n = d->nn;
		OutValue->l = (double *)malloc(sizeof(double)*d->nn);
		for(int i = 0; i < d->nn; ++i)
			OutValue->l[i] = (*probabilities[m])[i];
	}

	void TestScores(NBResult *v) {
		// perform the wilcoxon matched-pairs signed-rank test...
		v->kl_q = WilcoxonPaired(*r_kl_s, *s_kl_s, v->kl_err);
		v->b_q = WilcoxonPaired(*r_brier_s, *s_brier_s, v->b_err);
		v->er_q = WilcoxonPaired(*r_ca_s, *s_ca_s, v->er_err);
	}


	~NBCache() {
		DELETEVECTOR(atts);
		DELETEVECTOR(cond);
		DELETEVECTOR(interact);
		DELETEVECTOR(tinteract);
		delete kl_scores;
		delete brier_scores;
		delete ca_scores;
		delete s_ca_s;
		delete s_brier_s;
		delete s_kl_s;
		delete r_ca_s;
		delete r_brier_s;
		delete r_kl_s;
		for(int i = 0; i < d->card[d->na]; ++i)
			delete probabilities[i];
		delete[] probabilities;
	}

	void Grab(int a, int card, int *values) {
		int i;
		// replaces the attribute with new values
		delete (*atts)[a];
		delete (*cond)[a];
		d->card[a] = card;
		(*atts)[a] = new NBStats1(a,d->card[a],d->a);						// attribute frequencies
		(*cond)[a] = new NBStats2(a,d->na,d->card[a],d->card[d->na],d->a);	// conditionals
		for(i = 0; i < d->nn; ++i) {
			d->data[i][a] = values[i];
			(*atts)[a]->update(d->data[i]);
			(*cond)[a]->update(d->data[i]);
		}
	}
	
	NBCache(struct NBInput *input) {
		int i,j;

		d = input;
		s_kl_s = new vector<double>(d->nn);
		s_brier_s = new vector<double>(d->nn);
		s_ca_s = new vector<double>(d->nn);
		r_kl_s = new vector<double>(d->nn);
		r_brier_s = new vector<double>(d->nn);
		r_ca_s = new vector<double>(d->nn);
		kl_scores = new vector<double>(d->nn);
		brier_scores = new vector<double>(d->nn);
		ca_scores = new vector<double>(d->nn);
		atts = new vector<NBStats1 *>(d->k);
		cond = new vector<NBStats2 *>(d->na);
		interact = new vector<NBStats *>;
		tinteract = new vector<NBStats *>;
		for (i = 0; i < d->na; ++i) {
			(*atts)[i] = new NBStats1(i,d->card[i],d->a);							// attribute frequencies
			(*cond)[i] = new NBStats2(i,d->na,d->card[i],d->card[d->na],d->a);	// conditionals
		}
		(*atts)[d->na] = new NBStats1(d->na,d->card[d->na],d->a);				// class freq
		probabilities = new vector<double>*[d->card[d->na]];
		for(i = 0; i < d->card[d->na]; ++i)
			probabilities[i] = new vector<double>(d->nn);
	

/*
		for (j = 0; j < d->nn; ++j) {
			printf("%5d ",j);
			for (i = 0; i <= d->na; ++i) {
				printf("%d ",d->data[j][i]);
			}
			printf("\n");
		}
*/
		// update the statistics
		for (j = 0; j < d->nn; ++j) {
			(*atts)[d->na]->update(d->data[j]);
			for (i = 0; i < d->na; ++i) {
				(*atts)[i]->update(d->data[j]);
				(*cond)[i]->update(d->data[j]);
			}
		}
	}
};


void NBcleanup(struct NBInput *p) {
	int i;
	if (p->data != NULL) {
		for (i=0; i <= p->nn; ++i)
			free(p->data[i]);
		free(p->data);
	}
	if (p->l != NULL) {
		free(p->l);
	}
	if (p->cvi != NULL) {
		free(p->cvi);
	}
	free(p);
}



void NBkill(struct NBInfo *p) {
	NBcleanup(p->i);
	delete ((NBCache *)(p)->c);
	free(p);
}

void NBsaveScores(struct NBInfo *in) {
	((NBCache *)(in)->c)->SaveScores();
}
void NBrememberScores(struct NBInfo *in) {
	((NBCache *)(in)->c)->RememberScores();
}

void NBexportScores(struct NBInfo *in, int mode, struct NBList *OutValue) {
	((NBCache *)(in)->c)->ExportScores(mode, OutValue); // get the scores out
}

void NBexportProbabilities(struct NBInfo *in, int mode, struct NBList *OutValue) {
	((NBCache *)(in)->c)->ExportProbabilities(mode, OutValue); // get the scores out
}

void NBstoreModel(struct NBInfo *in, double *w, struct NBModel *m) {
	((NBCache *)(in)->c)->StoreModel(m,w); // remember the model
}

void NBcompareScores(struct NBInfo *in, struct NBResult *OutValue) {
	((NBCache *)(in)->c)->TestScores(OutValue);
	//OutValue->b_err = OutValue->kl_err = OutValue->er_err = 0.0; // no error estimation yet
}

void NBclassify(struct NBInfo *in, int *ex, struct NBList *OutValue) {
	OutValue->n = ((NBCache *)(in)->c)->d->card[((NBCache *)(in)->c)->d->na];
	OutValue->l = ((NBCache *)(in)->c)->Classify(ex);
}

void NBclassifyW(struct NBInfo *in, int *ex, struct NBList *OutValue) {
	OutValue->n = ((NBCache *)(in)->c)->d->card[((NBCache *)(in)->c)->d->na];
	OutValue->l = ((NBCache *)(in)->c)->ClassifyW(ex);
}

void NBquality(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue) {
	((NBCache *)(in)->c)->Evaluate(m,OutValue);
}

void TANquality(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue) {
	((NBCache *)(in)->c)->TANEvaluate(m,OutValue);
}

void NBqualityW(struct NBInfo *in, double *w, struct NBModel *m, struct NBResult *OutValue) {
	((NBCache *)(in)->c)->EvaluateW(m,w,OutValue);
}

void NBdivergence(struct NBInfo *in, struct NBModel *m, struct NBResult *OutValue) {
	((NBCache *)(in)->c)->Divergence(m,OutValue);
}


void NBupdate(struct NBInfo *in, int attribute, int card, int *values) {
	((NBCache *)(in)->c)->Grab(attribute, card, values);
}

// prepares the statistics
struct NBInfo *NBprepare(struct NBInput *input) {
	struct NBInfo *i;

	i = (struct NBInfo *)malloc(sizeof(struct NBInfo));

	// gather the total counts
	i->i = input;
	i->c = (void *)new NBCache(input);
	
	return i;
}

double NBcompareLists(int n, double *a, double *b) {
	int i;
	double t;
	vector<double> aa(n);
	vector<double> bb(n);

	for (i = 0; i < n; ++i) {
		aa[i] = a[i];
		bb[i] = b[i];
	}
	return WilcoxonPaired(aa,bb,t);
}
