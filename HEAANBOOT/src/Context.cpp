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

void Context::addBootContext(long logl, long logp) {
	if(bootKeyMap.find(logl) == bootKeyMap.end()) {
		ZZ p = power2_ZZ(logp);

		long loglh = logl/2;
		long l = 1 << logl;
		long Noverl = N >> (logl + 1);

		long lk = 1 << loglh;
		long lm = 1 << (logl - loglh);

		long i, j, idx, jdx;

		ZZX* pvec = new ZZX[l];
		ZZX* pvecInv = new ZZX[l];

		for (i = 0; i < l; ++i) {
			pvec[i].SetLength(N);
			pvecInv[i].SetLength(N);
		}

		CZZ* pvals = new CZZ[l];
		for (j = 0; j < l; ++j) {
			for (i = j; i < l; ++i) {
				long deg =((M - rotGroup[i]) * (i - j) * Noverl) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], logp);
			}
			for (i = 0; i < j; ++i) {
				long deg =((M - rotGroup[i]) * (l + i - j) * Noverl) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], logp);
			}

			fftSpecialInv(pvals, l);

			for (i = 0, idx = 0, jdx = Nh; i < l; ++i, idx += Noverl, jdx += Noverl) {
				pvec[j].rep[idx] = pvals[i].r;
				pvec[j].rep[jdx] = pvals[i].i;
			}
		}

		for (j = 0; j < l; ++j) {
			for (i = j; i < l; ++i) {
				long deg = (rotGroup[i - j] * i * Noverl) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], logp);
			}
			for (i = 0; i < j; ++i) {
				long deg = (rotGroup[l + i - j] * i * Noverl) % M;
				pvals[i] = EvaluatorUtils::evalCZZ(ksiPowsr[deg], ksiPowsi[deg], logp);
			}

			fftSpecialInv(pvals, l);

			for (i = 0, idx = 0, jdx = Nh; i < l; ++i, idx += Noverl, jdx += Noverl) {
				pvecInv[j].rep[idx] = pvals[i].r;
				pvecInv[j].rep[jdx] = pvals[i].i;
			}
		}

		delete[] pvals;

		for (i = 1; i < lm; ++i) {
			for (j = 0; j < lk; ++j) {
				 Ring2Utils::inpowerAndEqual(pvec[j + lk * i], rotGroup[l - lk * i], p, N);
				 Ring2Utils::inpowerAndEqual(pvecInv[j + lk * i], rotGroup[l - lk * i], p, N);
			}
		}

		cout << "x" << endl;
		BootContext bootContext(pvec, pvecInv, logp);
		cout << "y" << endl;
		bootKeyMap.insert(pair<long, BootContext>(logl, bootContext));
		cout << "z" << endl;
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
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < len / 2; ++j) {
				CZZ u = vals[i + j];
				CZZ v = vals[i + j + len / 2];
				RR tmp1 = to_RR(v.r) * (ksiPowsr[j * MoverLen] + ksiPowsi[j * MoverLen]);
				RR tmpr = tmp1 - to_RR(v.r + v.i) * ksiPowsi[j * MoverLen];
				RR tmpi = tmp1 + to_RR(v.i - v.r) * ksiPowsr[j * MoverLen];
				v.r = to_ZZ(tmpr);
				v.i = to_ZZ(tmpi);
				vals[i + j] = u + v;
				vals[i + j + len / 2] = u - v;
			}
		}
	}
}

void Context::fftInvLazy(CZZ* vals, const long size) {
	bitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		long MoverLen = M / len;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < len / 2; ++j) {
				CZZ u = vals[i + j];
				CZZ v = vals[i + j + len / 2];
				RR tmp1 = to_RR(v.r) * (ksiPowsr[(len - j) * MoverLen] + ksiPowsi[(len - j) * MoverLen]);
				RR tmpr = tmp1 - to_RR(v.r + v.i) * ksiPowsi[(len - j) * MoverLen];
				RR tmpi = tmp1 + to_RR(v.i - v.r) * ksiPowsr[(len - j) * MoverLen];
				v.r = to_ZZ(tmpr);
				v.i = to_ZZ(tmpi);
				vals[i + j] = u + v;
				vals[i + j + len / 2] = u - v;
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

