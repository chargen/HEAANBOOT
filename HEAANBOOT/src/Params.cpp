#include "Params.h"


Params::Params(long logN, long logq, double sigma, long h) : logN(logN), logq(logq), sigma(sigma), h(h) {
	N = 1 << logN;
}

long Params::suggestlogN(long lambda, long logq) {
	long res = ceil(logq * (lambda + 110) / 3.6);
	double logres = log2((double)res);
	return (long)ceil(logres);
}
