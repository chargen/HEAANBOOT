#include "EvaluatorUtils.h"


ZZ EvaluatorUtils::evalZZ(const double& x, const long& bits) {
	return evalZZ(to_RR(x), bits);
}

ZZ EvaluatorUtils::evalZZ(const RR& x, const long& bits) {
	RR xp = MakeRR(x.x, x.e + bits);
	return RoundToZZ(xp);
}

CZZ EvaluatorUtils::evalCZZ(const double& xr, const double& xi, const long& bits) {
	return evalCZZ(to_RR(xr), to_RR(xi), bits);
}

CZZ EvaluatorUtils::evalCZZ(const RR& xr, const RR& xi, const long& bits) {
	RR xrp = MakeRR(xr.x, xr.e + bits);
	RR xip = MakeRR(xi.x, xi.e + bits);
	return CZZ(RoundToZZ(xrp), RoundToZZ(xip));
}

CZZ EvaluatorUtils::evalRandCZZ(const long& bits) {
	return CZZ(RandomBits_ZZ(bits), RandomBits_ZZ(bits));
}

CZZ EvaluatorUtils::evalRandCZZ0(const long& bits) {
	return CZZ(RandomBits_ZZ(bits));
}

CZZ EvaluatorUtils::evalRandCZZCircle(const long& bits) {
	RR angle = random_RR();
	RR mr = cos(angle * 2 * Pi);
	RR mi = sin(angle * 2 * Pi);
	return evalCZZ(mr, mi, bits);
}

ZZ* EvaluatorUtils::evalRandZZArray(const long& size, const long& bits) {
	ZZ* res = new ZZ[size];
	for (long i = 0; i < size; i++) {
		res[i] = RandomBits_ZZ(bits);
	}
	return res;
}

CZZ* EvaluatorUtils::evalRandCZZArray(const long& size, const long& bits) {
	CZZ* res = new CZZ[size];
	for (long i = 0; i < size; i++) {
		res[i].r = RandomBits_ZZ(bits);
		res[i].i = RandomBits_ZZ(bits);
	}
	return res;
}

CZZ* EvaluatorUtils::evalRandCZZ0Array(const long& size, const long& bits) {
	CZZ* res = new CZZ[size];
	for (long i = 0; i < size; i++) {
		res[i].r = RandomBits_ZZ(bits);
	}
	return res;
}

CZZ EvaluatorUtils::evalCZZPow(const double& xr, const double& xi, const long& degree, const long& bits) {
	long logDegree = log2(degree);
	long po2Degree = 1 << logDegree;
	CZZ res = evalCZZPow2(xr, xi, logDegree, bits);
	long remDegree = degree - po2Degree;
	if(remDegree > 0) {
		CZZ tmp = evalCZZPow(xr, xi, remDegree, bits);
		res *= tmp;
		res >>= bits;
	}
	return res;
}

CZZ EvaluatorUtils::evalCZZPow(const RR& xr, const RR& xi, const long& degree, const long& bits) {
	long logDegree = log2(degree);
	long po2Degree = 1 << logDegree;
	CZZ res = evalCZZPow2(xr, xi, logDegree, bits);
	long remDegree = degree - po2Degree;
	if(remDegree > 0) {
		CZZ tmp = evalCZZPow(xr, xi, remDegree, bits);
		res *= tmp;
		res >>= bits;
	}
	return res;
}

CZZ EvaluatorUtils::evalCZZPow2(const double& xr, const double& xi, const long& logDegree, const long& bits) {
	return evalCZZPow2(to_RR(xr), to_RR(xi), logDegree, bits);
}

CZZ EvaluatorUtils::evalCZZPow2(const RR& xr, const RR& xi, const long& logDegree, const long& bits) {
	CZZ res = evalCZZ(xr, xi, bits);
	for (int i = 0; i < logDegree; ++i) {
		res *= res;
		res >>= bits;
	}
	return res;
}

CZZ* EvaluatorUtils::evalCZZPowArray(const double& xr, const double& xi, const long& degree, const long& bits) {
	return  evalCZZPowArray(to_RR(xr), to_RR(xi), degree, bits);
}

CZZ* EvaluatorUtils::evalCZZPowArray(const RR& xr, const RR& xi, const long& degree, const long& bits) {
	CZZ* res = new CZZ[degree];
	CZZ m = evalCZZ(xr, xi, bits);
	res[0] = m;
	for (long i = 0; i < degree - 1; ++i) {
		res[i + 1] = (res[i] * m) >> bits;
	}
	return res;
}

CZZ* EvaluatorUtils::evalCZZPow2Array(const double& xr, const double& xi, const long& logDegree, const long& bits) {
	return evalCZZPow2Array(to_RR(xr), to_RR(xi), logDegree, bits);
}

CZZ* EvaluatorUtils::evalCZZPow2Array(const RR& xr, const RR& xi, const long& logDegree, const long& bits) {
	CZZ* res = new CZZ[logDegree + 1];
	CZZ m = evalCZZ(xr, xi, bits);
	res[0] = m;
	for (long i = 0; i < logDegree; ++i) {
		res[i + 1] = (res[i] * res[i]) >> bits;
	}
	return res;
}

CZZ EvaluatorUtils::evalCZZInv(const double& xr, const double& xi, const long& bits) {
	return evalCZZInv(to_RR(xr), to_RR(xi), bits);
}

CZZ EvaluatorUtils::evalCZZInv(const RR& xr, const RR& xi, const long& bits) {
	RR xinvr = xr / (xr * xr + xi * xi);
	RR xinvi = -xi / (xr * xr + xi * xi);

	return evalCZZ(xinvr, xinvi, bits);
}

CZZ EvaluatorUtils::evalCZZLog(const double& xr, const double& xi, const long& bits) {
	double xlogr = log(xr * xr + xi * xi) / 2;
	double xlogi = atan(xi / xr);

	return evalCZZ(xlogr, xlogi, bits);
}

CZZ EvaluatorUtils::evalCZZExp(const double& xr, const double& xi, const long& bits) {
	double xrexp = exp(xr);
	double xexpr = xrexp * cos(xi);
	double xexpi = xrexp * sin(xi);

	return evalCZZ(xexpr, xexpi, bits);
}

CZZ EvaluatorUtils::evalCZZExp(const RR& xr, const RR& xi, const long& bits) {
	RR xrexp = exp(xr);
	RR xexpr = xrexp * cos(xi);
	RR xexpi = xrexp * sin(xi);

	return evalCZZ(xexpr, xexpi, bits);
}

CZZ EvaluatorUtils::evalCZZSigmoid(const double& xr, const double& xi, const long& bits) {
	double xrexp = exp(xr);
	double xexpr = xrexp * cos(xi);
	double xexpi = xrexp * sin(xi);

	double xsigmoidr = (xexpr * (xexpr + 1) + (xexpi * xexpi)) / ((xexpr + 1) * (xexpr + 1) + (xexpi * xexpi));
	double xsigmoidi = xexpi / ((xexpr + 1) * (xexpr + 1) + (xexpi * xexpi));

	return evalCZZ(xsigmoidr, xsigmoidi, bits);
}

CZZ EvaluatorUtils::evalCZZSigmoid(const RR& xr, const RR& xi, const long& bits) {
	RR xrexp = exp(xr);
	RR xexpr = xrexp * cos(xi);
	RR xexpi = xrexp * sin(xi);

	RR xsigmoidr = (xexpr * (xexpr + 1) + (xexpi * xexpi)) / ((xexpr + 1) * (xexpr + 1) + (xexpi * xexpi));
	RR xsigmoidi = xexpi / ((xexpr + 1) * (xexpr + 1) + (xexpi * xexpi));

	return evalCZZ(xsigmoidr, xsigmoidi, bits);
}

void EvaluatorUtils::leftShiftAndEqual(CZZ*& vals, const long& size, const long& bits) {
	for (long i = 0; i < size; ++i) {
		vals[i] <<= bits;
	}
}

void EvaluatorUtils::leftRotateAndEqual(CZZ*& vals, const long& size, const long& rotSize) {
	long remrotSize = rotSize % size;
	if(remrotSize != 0) {
		long divisor = GCD(remrotSize, size);
		long steps = size / divisor;
		for (long i = 0; i < divisor; ++i) {
			CZZ tmp = vals[i];
			long idx = i;
			for (long j = 0; j < steps - 1; ++j) {
				vals[idx] = vals[(idx + remrotSize) % size];
				idx = (idx + remrotSize) % size;
			}
			vals[idx] = tmp;
		}
	}
}

void EvaluatorUtils::rightRotateAndEqual(CZZ*& vals, const long& size, const long& rotSize) {
	long remrotSize = rotSize % size;
	long leftremrotSize = (size - remrotSize) % size;
	leftRotateAndEqual(vals, size, leftremrotSize);
}
