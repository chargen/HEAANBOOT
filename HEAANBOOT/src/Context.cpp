#include "Context.h"
#include "Ring2Utils.h"
#include "EvaluatorUtils.h"

Context::Context(Params& params) :
	logN(params.logN), logQ(params.logQ), sigma(params.sigma), h(params.h), N(params.N) {

	Nh = N >> 1;
	M = N << 1;
	logPQ = logQ << 1;
	Q = power2_ZZ(logQ);
	PQ = power2_ZZ(logPQ);

	rotGroup = new long[Nh];
	long fivePows = 1;
	for (long i = 0; i < Nh; ++i) {
		rotGroup[i] = fivePows;
		fivePows *= 5;
		fivePows %= M;
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

ZZX Context::encodeLarge(CZZ* vals, long slots) {
	CZZ* uvals = new CZZ[slots];
	long i, j, idx;
	for (i = 0; i < slots; ++i) {
		uvals[i] = vals[i] << logQ;
	}
	ZZX mx;
	mx.SetLength(N);

	long gap = Nh / slots;

	fftSpecialInv(uvals, slots);

	for (i = 0, j = Nh, idx = 0; i < slots; ++i, j += gap, idx += gap) {
		mx[idx] = uvals[i].r;
		mx[j] = uvals[i].i;
	}
	delete[] uvals;
	return mx;
}

ZZX Context::encodeSmall(CZZ* vals, long slots) {
	CZZ* uvals = new CZZ[slots];
	long i, jdx, idx;
	for (i = 0; i < slots; ++i) {
		uvals[i] = vals[i] << logQ;
	}
	ZZX mx;
	mx.SetLength(N);

	long gap = Nh / slots;

	fftSpecialInv(uvals, slots);

	for (i = 0, jdx = Nh, idx = 0; i < slots; ++i, jdx += gap, idx += gap) {
		mx[idx] = uvals[i].r >> logQ;
		mx[jdx] = uvals[i].i >> logQ;
	}
	delete[] uvals;
	return mx;
}

void Context::addBootContext(long logSlots, long logp) {
	if(bootContextMap.find(logSlots) == bootContextMap.end()) {
		ZZ p = power2_ZZ(logp);

		long dlogp = logp << 1;
		long logSlotsh = logSlots >> 1;
		long slots = 1 << logSlots;
		long gap = Nh >> logSlots;

		long lk = 1 << logSlotsh;

		long i, j, jdx, idx, deg;

		ZZX* pvec = new ZZX[slots];
		ZZX* pvecInv = new ZZX[slots];

		CZZ* pvals = new CZZ[slots];
		for (j = 0; j < slots; ++j) {
			for (i = j; i < slots; ++i) {
				deg = ((M - rotGroup[i]) * (i - j) * gap) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], dlogp);
			}
			for (i = 0; i < j; ++i) {
				deg =((M - rotGroup[i]) * (slots + i - j) * gap) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], dlogp);
			}
			fftSpecialInv(pvals, slots);
			pvec[j].SetLength(N);
			for (i = 0, jdx = Nh, idx = 0; i < slots; ++i, jdx += gap, idx += gap) {
				pvec[j][idx] = pvals[i].r >> logp;
				pvec[j][jdx] = pvals[i].i >> logp;
			}
		}

		for (j = 0; j < slots; ++j) {
			for (i = j; i < slots; ++i) {
				deg = (rotGroup[i - j] * i * gap) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], dlogp);
			}
			for (i = 0; i < j; ++i) {
				deg = (rotGroup[slots + i - j] * i * gap) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], dlogp);
			}
			fftSpecialInv(pvals, slots);
			pvecInv[j].SetLength(N);
			for (i = 0, jdx = Nh, idx = 0; i < slots; ++i, jdx += gap, idx += gap) {
				pvecInv[j][idx] = pvals[i].r >> logp;
				pvecInv[j][jdx] = pvals[i].i >> logp;
			}
		}

		delete[] pvals;

		for (i = lk; i < slots; i += lk) {
			for (j = 0; j < lk; ++j) {
				 Ring2Utils::inpowerAndEqual(pvec[j + i], rotGroup[slots - i], p, N);
				 Ring2Utils::inpowerAndEqual(pvecInv[j + i], rotGroup[slots - i], p, N);
			}
		}

		bootContextMap.insert(pair<long, BootContext>(logSlots, BootContext(pvec, pvecInv, logp)));
	}
}

void Context::bitReverse(CZZ* vals, const long size) {
	for (long i = 1, j = 0; i < size; ++i) {
		long bit = size >> 1;
		for (; j >= bit; bit>>=1) {
			j -= bit;
		}
		j += bit;
		if(i < j) {
			swap(vals[i], vals[j]);
		}
	}
}

void Context::fft(CZZ* vals, const long size) {
	bitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		long MoverLen = M / len;
		long lenh = len >> 1;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = j * MoverLen;
				CZZ u = vals[i + j];
				CZZ v = vals[i + j + lenh];
				RR tmp1 = to_RR(v.r) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.r + v.i) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.i - v.r) * ksiPowsr[idx];
				v.r = to_ZZ(tmpr);
				v.i = to_ZZ(tmpi);
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftInvLazy(CZZ* vals, const long size) {
	bitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		long MoverLen = M / len;
		long lenh = len >> 1;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = (len - j) * MoverLen;
				CZZ u = vals[i + j];
				CZZ v = vals[i + j + lenh];
				RR tmp1 = to_RR(v.r) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.r + v.i) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.i - v.r) * ksiPowsr[idx];
				v.r = to_ZZ(tmpr);
				v.i = to_ZZ(tmpi);
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftInv(CZZ* vals, const long size) {
	fftInvLazy(vals, size);
	for (long i = 0; i < size; ++i) {
		vals[i] /= size;
	}
}

void Context::fftSpecial(CZZ* vals, const long size) {
	bitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		for (long i = 0; i < size; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			for (long j = 0; j < lenh; ++j) {
				long idx = ((rotGroup[j] % lenq)) * M / lenq;
				CZZ u = vals[i + j];
				CZZ v = vals[i + j + lenh];
				RR tmp1 = to_RR(v.r) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.r + v.i) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.i - v.r) * ksiPowsr[idx];
				v.r = to_ZZ(tmpr);
				v.i = to_ZZ(tmpi);

				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftSpecialInvLazy(CZZ* vals, const long size) {
	for (long len = size; len >= 1; len >>= 1) {
		for (long i = 0; i < size; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			for (long j = 0; j < lenh; ++j) {
				long idx = (lenq - (rotGroup[j] % lenq)) * M / lenq;
				CZZ u = vals[i + j] + vals[i + j + lenh];
				CZZ v = vals[i + j] - vals[i + j + lenh];
				RR tmp1 = to_RR(v.r) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.r + v.i) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.i - v.r) * ksiPowsr[idx];
				v.r = to_ZZ(tmpr);
				v.i = to_ZZ(tmpi);

				vals[i + j] = u;
				vals[i + j + lenh] = v;
			}
		}
	}
	bitReverse(vals, size);
}

void Context::fftSpecialInv(CZZ* vals, const long size) {
	fftSpecialInvLazy(vals, size);
	for (long i = 0; i < size; ++i) {
		vals[i] /= size;
	}
}

Context::~Context() {
	delete[] rotGroup;
	delete[] ksiPowsi;
	delete[] ksiPowsr;
}

