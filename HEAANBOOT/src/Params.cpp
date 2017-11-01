#include "Params.h"

#include <cmath>

Params::Params(long logN, long logq, double sigma, long h) :
			logN(logN), logq(logq), sigma(sigma), h(h) {
	//-----------------------------------------
	N = 1 << logN;
	logqq = 2 * logq;
	q = power2_ZZ(logq);
	qq = power2_ZZ(logqq);

	rotGroup = new long[N / 2];

	long tmp = 1;
	long M = N << 1;
	for (long i = 0; i < N / 2; ++i) {
		rotGroup[i] = tmp;
		tmp *= 5;
		tmp %= M;
	}
}

long Params::suggestlogN(long lambda, long logq) {
	long res = ceil(logq * (lambda + 110) / 3.6);
	double logres = log2((double)res);
	return (long)ceil(logres);
}
