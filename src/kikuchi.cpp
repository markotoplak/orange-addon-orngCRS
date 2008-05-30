#include "kikuchi.h"
#include "models.h"
#include "../statistics/pr_statistics.h"
#include "cMersenneTwister.h"
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



KStat *KCache::makeaStat(int *p, int n) {
	Stats *s;
	KStat *ks;
	int j, k, pv;
	double prob;

	ks = new KStat;

	switch(n) {
		case 1:
			s = (Stats*)new Stats1(p[0],d->card[p[0]]); break;
		case 2:
			s = (Stats*)new Stats2(p[0],p[1],d->card[p[0]],d->card[p[1]]); break;
		case 3:
			s = (Stats*)new Stats3(p[0],p[1],p[2],d->card[p[0]],d->card[p[1]],d->card[p[2]]); break;
		case 4:
			s = (Stats*)new Stats4(p[0],p[1],p[2],p[3],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]]); break;
		case 5:
			s = (Stats*)new Stats5(p[0],p[1],p[2],p[3],p[4],d->card[p[0]],d->card[p[1]],d->card[p[2]],d->card[p[3]],d->card[p[4]]); break;
		default:
			fprintf(stderr, "Region of size %d is inappropriate.\n", n);
			return NULL;
	}
	// enter the data
	for (j = 0; j < d->nn; ++j) {
		s->update(d->data[j]);
	}
	s->set_prior(0);
	s->model();

	// get the predictions for all instances in the data
	ks->n = n;
	ks->s = s;
	ks->pred = (double*)malloc(sizeof(double)*d->nn*d->card[d->na]);
	for (j = 0; j < d->nn; ++j) {
		pv = d->data[j][d->na]; // store the original value
		for(k = 0; k < d->card[d->na]; ++k) { // for all label values
			d->data[j][d->na] = k;
			prob = s->getprob(d->data[j]);
			if(prob >= 1e-8)
				prob = log(prob);
			else
				prob = -18.5;
			ks->pred[j*d->card[d->na] + k] = prob;
		}
		d->data[j][d->na] = pv; // restore the original label value
	}
	return ks;
};

KStat *KCache::findaStat(int *p, int n) {
	int N, key, j, pp, idx;
	KStat *ks;

	N = d->k; // number of information sources

	// BUILD THE KEY
	key = 0; j = 0;
	pp = -1;
	do {
		assert(key < (2147483647-N)/N); // prevent overflow
		key *= N;
		key += 1+p[j];

		// assure that the ordering is right
		assert(p[j]>pp);
		pp = p[j];
	} while (++j < n);
	found = lookup.find(key);
	if (found == lookup.end()) {
		// not found
		ks = makeaStat(p,n);

		idx = region_bank->size();
		lookup[key] = idx;
		ks->idx = idx; // remember the index...
		region_bank->push_back(ks);
		return ks;
	} else {
		// found
		idx = found->second;
		return (*region_bank)[idx];
	}
}

void KCache::addModel(struct KModel *m) {
	int i,j,k,idx;
	double sum;
	
	// build the regions
	for(i = 0; i < m->nogroups; ++i) {
		int *p = m->groups[i].l;
		KRegion *r = new KRegion;
		r->magnitude = m->groups[i].times;
		r->ks = findaStat(p,m->groups[i].n);
		if(m->groups[i].l[m->groups[i].n-1] != d->na){ // check that the last attribute is the label
			//printf("unlabelled!\n\n");
			r->ks->s->set_unlabelled();
		}
		regions->push_back(r);

		// update degrees of freedom
		dof += r->magnitude * r->ks->s->getDOF();
		adof += r->magnitude * r->ks->s->getDOFsup(d->card[d->na]);
	}

	// build the classifications
	loss = 0.0;
	for (j = 0; j < d->nn; ++j) {
		for(k = 0; k < d->card[d->na]; ++k) { // for all label values
			idx = j*d->card[d->na] + k;
			sum = 0.0;
			for(i = 0; i < regions->size(); ++i) { // for all regions
				sum += (*regions)[i]->magnitude * (*regions)[i]->ks->pred[idx];
			}
			results[k] = saved[idx] = predictions[idx] = sum;
		}
		deLogize(temp);
		loss -= log(temp[d->data[j][d->na]]);
	}
}

void KCache::testModel(struct KModel *m, double *result, double *outdof) {
	int i,j,k, idx, tadof, tdof;
	double sum,error;
	vector<KRegion *> myregs;
	
	// prepare the temporary regions
	for(i = 0; i < m->nogroups; ++i) {
		int *p = m->groups[i].l;
		KRegion *r = new KRegion;
		assert(m->groups[i].l[m->groups[i].n-1] == d->na); // check that the last attribute is the label
		r->magnitude = m->groups[i].times;
		r->ks = findaStat(p,m->groups[i].n);
		myregs.push_back(r);
	}

	// build the classifications
	error = 0.0;
	for (j = 0; j < d->nn; ++j) {
		for(k = 0; k < d->card[d->na]; ++k) { // for all label values
			idx = j*d->card[d->na] + k;
			sum = saved[idx];
			for(i = 0; i < myregs.size(); ++i) { // for all new regions
				sum += myregs[i]->magnitude * myregs[i]->ks->pred[idx];
			}
			results[k] = sum;
		}
		deLogize(temp);
		// compute the log-likelihood, etc.
		error -= log(temp[d->data[j][d->na]]);
	}

	tdof = 0;
	tadof = 0;
	for(i = 0; i < myregs.size(); ++i) { // for all regions
		//printf("%d %d\n",myregs[i]->ks->s->getDOF(),myregs[i]->ks->s->getDOFsup(d->card[d->na]));
		tdof += myregs[i]->magnitude * myregs[i]->ks->s->getDOF();
		tadof += myregs[i]->magnitude * myregs[i]->ks->s->getDOFsup(d->card[d->na]);
	}
	//printf("-%d %d-\n",tdof,tadof);

	*result = error;
	*outdof = tdof-tadof;
}

void KCache::emptyModel() {
	int i;
	dof = 0;
	adof = 0;
	if(regions != NULL) {
		for (i = 0; i < regions->size(); ++i)
			delete (*regions)[i];
		delete regions;
	}
	regions = new vector<KRegion *>;
}

void KCache::emptyEnsemble() {
	int i;

	if (ensemble != NULL) {
		for (i = 0; i < ensemble->regions->size(); ++i) {
			free(ensemble->predictions[i]);
			delete (*ensemble->regions)[i];
		}
		free(ensemble->regions);
		free(ensemble->predictions);
		for (i = 0; i < ensemble->n; ++i) {
			delete[] ((*ensemble->indices)[i]).indices;
			delete[] ((*ensemble->indices)[i]).magnitudes;
		}
		delete ensemble->indices;
		delete ensemble->weights;
		delete ensemble;
		ensemble = NULL;
	}
}

void KCache::setEnsemble(struct KModels *m, struct KArray *weights) {
	int i,j,k,idx;
	double sum;
	map<long,int> LUT;
	map<long,int>::iterator LUTi;
	
	emptyEnsemble();
	ensemble = new KEnsemble;
	ensemble->regions = new vector<KRegion *>;
	
	// size it up
	ensemble->n = m->nomodels;
	ensemble->weights = new vector<double>(ensemble->n);
	ensemble->indices = new vector<KModelCompact>(ensemble->n);
	assert(ensemble->n == weights->n);
	
	// build the regions, but only storing the region once.
	for(k = 0; k < ensemble->n; ++k) {
		int ng;
		ng = m->models[k].nogroups;
		(*ensemble->weights)[k] = weights->l[k];
		(*ensemble->indices)[k].n = ng;
		(*ensemble->indices)[k].indices = new int[ng];
		(*ensemble->indices)[k].magnitudes = new int[ng];
		for(i = 0; i < ng; ++i) {
			KStat *ks;

			ks = findaStat(m->models[k].groups[i].l,m->models[k].groups[i].n);
			if(m->models[k].groups[i].l[m->models[k].groups[i].n-1] != d->na){ // check that the last attribute is the label
				ks->s->set_unlabelled();
			}
			// find the region
			LUTi = LUT.find(ks->idx);
			if (LUTi == LUT.end()) {
				// new region in this ensemble
				KRegion *r = new KRegion;
				r->magnitude = 0;
				r->ks = ks;
				idx = ensemble->regions->size();
				LUT[ks->idx] = idx;
				ensemble->regions->push_back(r);
			} else {
				// found
				idx = found->second;
			}
			// form the model
			(*ensemble->indices)[k].indices[i] = idx;
			(*ensemble->indices)[k].magnitudes[i] = m->models[k].groups[i].times;
		}
	}
	// allocate the memory for the predictions of individual regions
	ensemble->predictions = (double **)malloc(sizeof(double *)*ensemble->regions->size());
	for(k = 0; k < ensemble->regions->size(); ++k)
		ensemble->predictions[k] = (double *)malloc(sizeof(double)*d->card[d->na]); // for every class value, a prediction

}

void KCache::ClassifyEnsemble(int *ex, double *outresult) {
	int i, idx;
	unsigned int j,k;
	double sum, prob, lprob, pred;

	// for all the classes
	for(i = 0; i < d->card[d->na]; ++i){
		ex[d->na] = i;
		for(j = 0; j < ensemble->regions->size(); ++j) {
			prob = (*ensemble->regions)[j]->ks->s->getprob(ex);
			if(prob >= 1e-8)
				lprob = log(prob);
			else
				lprob = -18.5;
			ensemble->predictions[j][i] = lprob;
		}
		outresult[i] = 0.0; // the running average
	}
	// for all models
	for(i = 0; i < ensemble->n; ++i) {
		// for all values
		for(j = 0; j < d->card[d->na]; ++j) {
			// start with zero
			results[j] = 0.0;
			// for all groups
			for (k = 0; k < (*ensemble->indices)[i].n; ++k) {
				idx = (*ensemble->indices)[i].indices[k];
				pred = ensemble->predictions[idx][j];
				results[j] +=  pred * (*ensemble->indices)[i].magnitudes[k];
			}
		}
		deLogize(temp);
		//printf("\n** %d: (%f)",i,(*ensemble->weights)[i]);
		for(j = 0; j < d->card[d->na]; ++j) {
			//printf("%f ",temp[j]);
			outresult[j] += (*ensemble->weights)[i]*temp[j];
		}
	}
	// re-normalize
	sum = 0.0;
	for(i = 0; i < d->card[d->na]; ++i)
		sum += outresult[i];
	sum = 1.0/sum;
	//printf("\ntotal: ");
	for(i = 0; i < d->card[d->na]; ++i) {
		outresult[i] *= sum;
		//printf("%f ",outresult[i]);
	}
	//printf("\n");
}

// convert the odds into probabilities for a given case
// gets the input (odds) in results[], writes the probabilities in outresults
void KCache::deLogize(double *outresult) {
	int i,j,maxi;
	double sum, norm, maxo;

	maxo = -1e200;
	for(i = 0; i < d->card[d->na]; ++i){
		maxo = maxx(maxo,results[i]);
	}

	sum = 0.0;
	for(i = 0; i < d->card[d->na]; ++i){
		best[i] = 1.0;
		norm = exp(results[i]-maxo);
		sum += norm;
		outresult[i] = norm;
	}

	if (sum < 1e-200) {
		// normalize in log-space
		// cutting ties
		best[0] = 1.0;
		for(i = 1; i < d->card[d->na]; ++i) {
			best[i] = 0;
		}
		norm = log(sum);
		maxi = 0;
		maxo = results[0];
		sum = 0.0;
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
			outresult[i] = exp(results[i]-norm);
			sum += outresult[i];
		}
	}
	if(sum < 1e-200) {
		// everything failed, send the backup
		sum = 0.0;
		for(i = 0; i < d->card[d->na]; ++i)
			sum += best[i];
		sum = 1.0/sum;
		for(i = 0; i < d->card[d->na]; ++i)
			outresult[i] = best[i] * sum;
	} else {
		sum = 1.0/sum;
		for(i = 0; i < d->card[d->na]; ++i) {
			outresult[i] = outresult[i]*sum;
		}
	}
}

void KCache::Classify(int *ex, double *outresult) {
	int i;
	unsigned int j;
	double odds, prob, lprob;

	// for all the classes
	for(i = 0; i < d->card[d->na]; ++i){
		ex[d->na] = i;
		odds = 0.0;
		for(j = 0; j < regions->size(); ++j) {
			prob = (*regions)[j]->ks->s->getprob(ex);
			if(prob >= 1e-8)
				lprob = log(prob);
			else
				lprob = -18.5;
			odds += (*regions)[j]->magnitude*lprob;
			//printf("(%d %f) ",(*regions)[j]->magnitude,lprob);
		}
		results[i] = odds;
		//printf("%d => %f\n",i,odds);
	}
	deLogize(outresult);
	//for(i = 0; i < d->card[d->na]; ++i)
	//	printf("%d => %f\n",i,outresult[i]);
}

// get the log-loss
double KCache::getLoss(int *ex) {
	int i, pv;
	unsigned int j;
	double odds, prob, lprob;

	// for all the classes
	pv = ex[d->na];
	for(i = 0; i < d->card[d->na]; ++i){
		ex[d->na] = i;
		odds = 0.0;
		for(j = 0; j < regions->size(); ++j) {
			prob = (*regions)[j]->ks->s->getprob(ex);
			if(prob >= 1e-8)
				lprob = log(prob);
			else
				lprob = -18.5;
			odds += (*regions)[j]->magnitude*lprob;
		}
		results[i] = odds;
	}
	deLogize(temp);
	return -log(temp[pv]);
}

KCache::~KCache() {
	unsigned int i;
	for (i = 0; i < region_bank->size(); ++i) {
		delete (*region_bank)[i]->s; 
		free((*region_bank)[i]->pred);
	}
	delete region_bank;
	delete regions;
	emptyEnsemble();
	free(results); free(best); free(temp); free(predictions); free(saved);
}

KCache::KCache(struct KInput *input, double prior) {
	d = input;
	regions = NULL;
	region_bank = new vector<KStat *>;

	results = (double *)malloc(sizeof(double)*d->card[d->na]);
	best = (double *)malloc(sizeof(double)*d->card[d->na]);
	temp = (double *)malloc(sizeof(double)*d->card[d->na]);
	predictions = (double *)malloc(sizeof(double)*d->nn*d->card[d->na]);
	saved = (double *)malloc(sizeof(double)*d->nn*d->card[d->na]);
	ensemble = NULL;

	dof = 0;
	adof = 0;
	this->prior = prior;
}


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

// stores the model and learns the statistics
void Ksetmodel(struct KInfo *in, struct KModel *m) {
	((KCache *)in->c)->emptyModel();
	((KCache *)in->c)->addModel(m);

}

// stores the model and learns the statistics
void Kaddmodel(struct KInfo *in, struct KModel *m) {
	((KCache *)in->c)->addModel(m);
}

// stores the model and learns the statistics
void Ktestaddition(struct KInfo *in, struct KModel *m, struct KList *OutValue) {
	OutValue->n = 2;
	OutValue->l = (double *)malloc(sizeof(double)*2);
	((KCache *)in->c)->testModel(m,OutValue->l,OutValue->l+1);
}

// get DOF
void KgetDOF(struct KInfo *in, struct KList *OutValue) {
	OutValue->n = 3;
	OutValue->l = (double *)malloc(sizeof(double)*3);
	OutValue->l[0] = ((KCache *)in->c)->dof;
	OutValue->l[1] = ((KCache *)in->c)->adof;
	OutValue->l[2] = ((KCache *)in->c)->loss;
}

// cleanup
void Kdie(struct KInfo *p) {
	Kcleanup(p->i);
	delete ((KCache *)(p)->c);
	free(p);
}

void Ksetensemble(struct KInfo *in, struct KModels *ms, struct KArray *weights) {
	assert(ms->nomodels == weights->n);
	((KCache *)in->c)->setEnsemble(ms,weights);
}

void Kuseensemble(struct KInfo *in, int *ex, struct KList *OutValue) {
	int card = ((KCache *)(in)->c)->d->card[((KCache *)(in)->c)->d->na];
	OutValue->n = card;
	OutValue->l = (double *)malloc(sizeof(double)*card);
	((KCache *)(in)->c)->ClassifyEnsemble(ex,OutValue->l);		
}
	
// classify
void Kuse(struct KInfo *in, int *ex, struct KList *OutValue) {
	int card = ((KCache *)(in)->c)->d->card[((KCache *)(in)->c)->d->na];
	OutValue->n = card;
	OutValue->l = (double *)malloc(sizeof(double)*card);
	((KCache *)(in)->c)->Classify(ex,OutValue->l);
}

double Kvalidate(struct KInfo *in, int *ex) {
	return ((KCache *)(in)->c)->KCache::getLoss(ex);
}

void Kcheckreversal(struct KInfo *in, struct KList *OutValue) {
	OutValue->n = 2;
	OutValue->l = (double *)malloc(sizeof(double)*2);
	((KCache *)(in)->c)->CheckReversal(OutValue->l);
}

// creates the data structure
struct KInfo *Kremember(struct KInput *input, double prior) {
	struct KInfo *i;

	i = (struct KInfo *)malloc(sizeof(struct KInfo));

	i->i = input;
	i->c = (void *)new KCache(i->i,prior);

	return i;
}

void KCache::CheckReversal(double *outresult) {
	int j,k,idx;
	double error,error_r,sum;

	error = error_r = 0.0;
	for (j = 0; j < d->nn; ++j) {
		for(k = 0; k < d->card[d->na]; ++k) { // for all label values
			idx = j*d->card[d->na] + k;
			results[k] = predictions[idx];
		}
		deLogize(temp);
		error -= log(temp[d->data[j][d->na]]);
		// revert
		sum = 0.0;
		for(k = 0; k < d->card[d->na]; ++k) {
			temp[k] = 1.0-temp[k];
			sum += temp[k];
		}
		error_r -= log(temp[d->data[j][d->na]]/sum);
	}
	outresult[0] = error;
	outresult[1] = error_r;
}
