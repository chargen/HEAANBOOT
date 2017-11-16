#ifndef HEAAN_CONTEXT_H_
#define HEAAN_CONTEXT_H_

#include <NTL/RR.h>

#include "Common.h"
#include "Params.h"

using namespace std;
using namespace NTL;

static RR const Pi = ComputePi_RR();

static string LOGARITHM = "Logarithm"; ///< log(x)
static string EXPONENT  = "Exponent"; ///< exp(x)
static string SIGMOID   = "Sigmoid"; ///< sigmoid(x) = exp(x) / (1 + exp(x))

class Context {
public:

	long logN; ///< N is a power of 2 that corresponds to the ring Z[X] / (X^N + 1)
	long logq; ///< q corresponds to the highest modulus
	double sigma; ///< sigma corresponds to standard deviation for error and secret key coefficients generation from Gaussian distribution
	long h; ///< hamming weight of secret key

	long N;
	long M;
	long logqq; ///< qq = q * q

	ZZ q;
	ZZ qq;

	long* rotGroup; ///< auxiliary information about rotation group indexes for batch encoding
	RR* ksiPowsr; ///< storing ksi pows for fft calculation
	RR* ksiPowsi; ///< storing ksi pows for fft calculation
	map<string, double*> taylorCoeffsMap; ///< storing taylor coefficients for function calculation

	Context(Params& params);

	virtual ~Context();
};

#endif /* CONTEXT_H_ */
