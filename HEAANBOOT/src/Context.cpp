#include "Context.h"


Context::Context(Params& params) :
	logN(params.logN), logq(params.logq), sigma(params.sigma), h(params.h), N(params.N) {

	M = N << 1;
	logqq = 2 * logq;
	q = power2_ZZ(logq);
	qq = power2_ZZ(logqq);

	rotGroup = new long[N / 2];
	long tmp = 1;
	for (long i = 0; i < N / 2; ++i) {
		rotGroup[i] = tmp;
		tmp *= 5;
		tmp %= M;
	}

	ksiPowsr = new RR[M + 1];
	ksiPowsi = new RR[M + 1];

	for (long j = 0; j < M; ++j) {
		RR angle = 2.0 * Pi * j / M;
		ksiPowsr[j] = cos(angle);
		ksiPowsi[j] = sin(angle);
	}

	ksiPowsr[M] = ksiPowsr[0];
	ksiPowsi[M] = ksiPowsi[0];

	taylorCoeffsMap.insert(pair<string, double*>(LOGARITHM, new double[11]{0,1,-0.5,1./3,-1./4,1./5,-1./6,1./7,-1./8,1./9,-1./10}));
	taylorCoeffsMap.insert(pair<string, double*>(EXPONENT, new double[11]{1,1,0.5,1./6,1./24,1./120,1./720,1./5040, 1./40320,1./362880,1./3628800}));
	taylorCoeffsMap.insert(pair<string, double*>(SIGMOID, new double[11]{1./2,1./4,0,-1./48,0,1./480,0,-17./80640,0,31./1451520,0}));
}

Context::~Context() {
	delete[] rotGroup;
	delete[] ksiPowsi;
	delete[] ksiPowsr;
}

