#ifndef HEAAN_SCHEMEAUX_H_
#define HEAAN_SCHEMEAUX_H_

#include "CZZ.h"
#include "Params.h"
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <iostream>
#include <map>

using namespace std;
using namespace NTL;

static RR const Pi = ComputePi_RR();

class SchemeAux {
public:

	long M;
	RR* ksiPowsr; ///< storing ksi pows for fft calculation
	RR* ksiPowsi; ///< storing ksi pows for fft calculation
	map<string, double*> taylorCoeffsMap; ///< storing taylor coefficients for function calculation

	//-----------------------------------------

	SchemeAux(long& logN);

	//-----------------------------------------

};

static string LOGARITHM = "Logarithm"; ///< log(x)
static string EXPONENT  = "Exponent"; ///< exp(x)
static string SIGMOID   = "Sigmoid"; ///< sigmoid(x) = exp(x) / (1 + exp(x))

#endif
