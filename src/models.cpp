#include "models.h"

// compare the trained classifier using a degree-1 model with "reality"
double KirkwoodClass(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex) {
	int idx,j, x1, x2, x3;
	double odds, prob, sum, tprob, todds, codds, tcodds, tsum, divergence;

	assert(n_att >= 0 && n_att <= 4);
	
	// iterate through all the combinations of attribute values
	for(idx = 0; idx < n_att; ++idx) {
		ex[I[idx]] = 0;
	}

	divergence = 0.0;
	do {
		idx = n_att-1;
		// iterate through label values
		sum = 0.0;
		for(j = 0; j < cards[C]; ++j) {
			ex[C] = j;
			// obtain the probabilities
			odds = 0.0;
			odds += safelog((*stat_train)[0]->getprob(ex));
			if(n_att > 1) {
				for(x1 = 0; x1 < n_att; ++x1) {
					odds -= safelog((*stat_train)[LUT[1][I[x1]]]->getprob(ex));
					if(n_att > 2) {
						for(x2 = x1+1; x2 < n_att; ++x2) {
							odds += safelog((*stat_train)[LUT[2][I[x1]+I[x2]*C]]->getprob(ex));
							if(n_att > 3) {
								for(x3 = x2+1; x3 < n_att; ++x3) {
									odds -= safelog((*stat_train)[LUT[3][I[x1]+(I[x2]+I[x3]*C)*C]]->getprob(ex));
								}
							}
						}
					}
				}
			}
			// fix polarity appropriately
			if (n_att == 2 || n_att == 4) {
				odds = -odds;
			}

			// get the probability
			prob = exp(odds);
			sum += prob;
			probs[j] = prob;
			oddss[j] = odds;
		}
		assert(sum > 1e-20);
		
		// normalize
		tsum = 1.0/sum;
		for(j = 0; j < cards[C]; ++j)
			probs[j] *= tsum;
		codds = safelog(sum);

		// verify the quality
		tsum = 0.0;
		for(j = 0; j < cards[C]; ++j) {
			ex[C] = j;
			// obtain the 'true' probability
			switch(n_att) {
				case 0: tprob = ((*stat_test)[0]->getprob(ex)); break;
				case 1: tprob = ((*stat_test)[LUT[1][I[0]]]->getprob(ex)); break;
				case 2:	tprob = ((*stat_test)[LUT[2][I[0]+I[1]*C]]->getprob(ex)); break;
				case 3: tprob = ((*stat_test)[LUT[3][I[0]+(I[1]+I[2]*C)*C]]->getprob(ex)); break;
				case 4: tprob = ((*stat_test)[LUT[4][I[0]+(I[1]+(I[2]+I[3]*C)*C)*C]]->getprob(ex)); break;
			}
			tsum += tprob;
			todds = safelog(tprob);
			tprobs[j] = tprob;
			toddss[j] = todds;
		}
		tcodds = safelog(tsum);
		for(j = 0; j < cards[C]; ++j) {
			// compute the divergence
			divergence += tprobs[j]*(toddss[j]-oddss[j] + codds - tcodds);
		}

		// move to next
		if(idx < 0)
			break;
		++ex[I[idx]];
		while(idx > 0 && ex[I[idx]] >= cards[I[idx]]) {
			ex[I[idx]] = 0;
			++ex[I[--idx]];
		}
	} while(ex[I[idx]] < cards[I[idx]]);
	return divergence;
}

// compare the trained non-simplified model with reality
double KirkwoodInt(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex) {
	int idx,j;
	double odds, prob, sum, tprob, todds, codds, tcodds, tsum, divergence;

	assert(n_att >= 0 && n_att <= 4);
	
	// iterate through all the combinations of attribute values
	for(idx = 0; idx < n_att; ++idx) {
		ex[I[idx]] = 0;
	}

	divergence = 0.0;
	do {
		// iterate through label values
		idx = n_att-1;
		sum = 0.0;
		tsum = 0.0;
		for(j = 0; j < cards[C]; ++j) {
			ex[C] = j;
			// obtain the 'true' probability
			switch(n_att) {
				case 0: 
					tprob = ((*stat_test)[0]->getprob(ex)); 
					prob = ((*stat_train)[0]->getprob(ex)); 
					break;
				case 1: 
					tprob = ((*stat_test)[LUT[1][I[0]]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[1][I[0]]]->getprob(ex)); 
					break;
				case 2:	
					tprob = ((*stat_test)[LUT[2][I[0]+I[1]*C]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[2][I[0]+I[1]*C]]->getprob(ex)); 
					break;
				case 3: 
					tprob = ((*stat_test)[LUT[3][I[0]+(I[1]+I[2]*C)*C]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[3][I[0]+(I[1]+I[2]*C)*C]]->getprob(ex)); 
					break;
				case 4: 
					tprob = ((*stat_test)[LUT[4][I[0]+(I[1]+(I[2]+I[3]*C)*C)*C]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[4][I[0]+(I[1]+(I[2]+I[3]*C)*C)*C]]->getprob(ex)); 
					break;
			}
			tsum += tprob;
			sum += prob;
			odds = safelog(prob);
			todds = safelog(tprob);
			tprobs[j] = tprob;
			probs[j] = prob;
			toddss[j] = todds;
			oddss[j] = odds;
		}
		// normalization factors due to conditionalization
		tcodds = safelog(tsum);
		codds = safelog(sum);
		for(j = 0; j < cards[C]; ++j) {
			// compute the divergence
			divergence += tprobs[j]*(toddss[j]-oddss[j]+codds - tcodds);
		}

		// move to next
		if(idx < 0)
			break;
		++ex[I[idx]];
		while(idx > 0 && ex[I[idx]] >= cards[I[idx]]) {
			ex[I[idx]] = 0;
			++ex[I[--idx]];
		}
	} while(ex[I[idx]] < cards[I[idx]]);
	return divergence;
}

// compare the trained non-simplified model with reality
double KirkwoodJoint(const int n_att, const int C, const int *I, const int *cards, double *probs, double *tprobs, double *oddss, double *toddss, vector <Stats *> *stat_train, vector <Stats *> *stat_test, vector<int *> &LUT, int *ex) {
	int idx,j;
	double odds, prob, tprob, todds, divergence;

	assert(n_att >= 0 && n_att <= 4);
	
	// iterate through all the combinations of attribute values
	for(idx = 0; idx < n_att; ++idx) {
		ex[I[idx]] = 0;
	}

	divergence = 0.0;
	do {
		idx = n_att-1;
		for(j = 0; j < cards[C]; ++j) {
			ex[C] = j;
			// obtain the 'true' probability
			switch(n_att) {
				case 0: 
					tprob = ((*stat_test)[0]->getprob(ex)); 
					prob = ((*stat_train)[0]->getprob(ex)); 
					break;
				case 1: 
					tprob = ((*stat_test)[LUT[1][I[0]]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[1][I[0]]]->getprob(ex)); 
					break;
				case 2:	
					tprob = ((*stat_test)[LUT[2][I[0]+I[1]*C]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[2][I[0]+I[1]*C]]->getprob(ex)); 
					break;
				case 3: 
					tprob = ((*stat_test)[LUT[3][I[0]+(I[1]+I[2]*C)*C]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[3][I[0]+(I[1]+I[2]*C)*C]]->getprob(ex)); 
					break;
				case 4: 
					tprob = ((*stat_test)[LUT[4][I[0]+(I[1]+(I[2]+I[3]*C)*C)*C]]->getprob(ex)); 
					prob = ((*stat_train)[LUT[4][I[0]+(I[1]+(I[2]+I[3]*C)*C)*C]]->getprob(ex)); 
					break;
			}
			odds = safelog(prob);
			todds = safelog(tprob);
			tprobs[j] = tprob;
			toddss[j] = todds;
			oddss[j] = odds;
		}
		for(j = 0; j < cards[C]; ++j) {
			// compute the divergence
			divergence += tprobs[j]*(toddss[j]-oddss[j]);
		}

		// move to next
		if(idx < 0)
			break;
		++ex[I[idx]];
		while(idx > 0 && ex[I[idx]] >= cards[I[idx]]) {
			ex[I[idx]] = 0;
			++ex[I[--idx]];
		}
	} while(ex[I[idx]] < cards[I[idx]]);
	return divergence;
}

// compare the trained classifier using a degree-1 model with "reality"
double KikuchiClass(const int *cards, KCModel *m, vector<KCRegion *> &regions, double *probs, double *tprobs, double *oddss, double *toddss, int *ex) {
	int idx,i,j, n_att, C;
	int *I;
	double odds, prob, sum, tprob, todds, codds, tcodds, tsum, divergence, lprob;

	// iterate through all the combinations of attribute values
	n_att = m->na-1;
	I = m->indices;
	C = I[n_att];
	for(idx = 0; idx < n_att; ++idx)
		ex[I[idx]] = 0;

	divergence = 0.0;
	do {
		idx = n_att-1;

		// iterate through label values
		sum = 0.0;
		for(j = 0; j < cards[C]; ++j) {
			ex[C] = j;
			// obtain the probabilities
			odds = 0.0;
			for(i = 1; i < m->n; ++i) {
				// skip first -- original
				prob = regions[m->stats[i]]->train->getprob(ex);
				if(prob >= 1e-8)
					lprob = log(prob);
				else
					lprob = -18.5;
				odds += m->magnitudes[i]*lprob;
			}

			// get the probability
			prob = exp(odds);
			sum += prob;
			probs[j] = prob;
			oddss[j] = odds;
		}
		assert(sum > 1e-20);
		
		// normalize
		tsum = 1.0/sum;
		for(j = 0; j < cards[C]; ++j)
			probs[j] *= tsum;
		codds = safelog(sum);

		// verify the quality
		tsum = 0.0;
		for(j = 0; j < cards[C]; ++j) {
			ex[C] = j;
			tprob = regions[m->stats[0]]->test->getprob(ex);
			tsum += tprob;
			todds = safelog(tprob);
			tprobs[j] = tprob;
			toddss[j] = todds;
		}
		tcodds = safelog(tsum);
		for(j = 0; j < cards[C]; ++j) {
			// compute the divergence
			divergence += tprobs[j]*(toddss[j]-oddss[j] + codds - tcodds);
		}

		// move to next
		if(idx < 0)
			break;
		++ex[I[idx]];
		while(idx > 0 && ex[I[idx]] >= cards[I[idx]]) {
			ex[I[idx]] = 0;
			++ex[I[--idx]];
		}
	} while(ex[I[idx]] < cards[I[idx]]);
	return divergence;
}

// compare the trained non-simplified model with reality
double KikuchiInt(const int *cards, KCModel *m, vector<KCRegion *> &regions, double *probs, double *tprobs, double *oddss, double *toddss, int *ex) {
	int idx,j, n_att, C;
	int *I;
	double odds, prob, sum, tprob, todds, codds, tcodds, tsum, divergence;

	// iterate through all the combinations of attribute values
	n_att = m->na-1;
	I = m->indices;
	C = I[n_att];
	for(idx = 0; idx < n_att; ++idx)
		ex[I[idx]] = 0;

	divergence = 0.0;
	do {
		idx = n_att-1; // do not modify the class...
		// verify the quality
		sum = 0.0;
		tsum = 0.0;
		for(j = 0; j < cards[C]; ++j) {
			ex[C] = j;
			prob = regions[m->stats[0]]->train->getprob(ex);
			tprob = regions[m->stats[0]]->test->getprob(ex);
			tsum += tprob;
			sum += prob;
			odds = safelog(prob);
			todds = safelog(tprob);
			tprobs[j] = tprob;
			probs[j] = prob;
			toddss[j] = todds;
			oddss[j] = odds;
		}
		// normalization factors due to conditionalization
		tcodds = safelog(tsum);
		codds = safelog(sum);
		for(j = 0; j < cards[C]; ++j) {
			// compute the divergence
			divergence += tprobs[j]*(toddss[j]-oddss[j]+codds - tcodds);
		}

		// move to next
		if(idx < 0)
			break;
		++ex[I[idx]];
		while(idx > 0 && ex[I[idx]] >= cards[I[idx]]) {
			ex[I[idx]] = 0;
			++ex[I[--idx]];
		}
	} while(ex[I[idx]] < cards[I[idx]]);
	return divergence;
}