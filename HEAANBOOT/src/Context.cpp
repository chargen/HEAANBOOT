#include "Context.h"
#include "Ring2Utils.h"
#include "EvaluatorUtils.h"

#include "StringUtils.h"

Context::Context(Params& params) :
	logN(params.logN), logQ(params.logQ), sigma(params.sigma), h(params.h), N(params.N) {

	Nh = N >> 1;
	logNh = logN - 1;
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

	qpowvec = new ZZ[logPQ + 1];

	qpowvec[0] = ZZ(1);
	for (long i = 1; i < logPQ + 1; ++i) {
		qpowvec[i] = qpowvec[i - 1] << 1;
	}

	taylorCoeffsMap.insert(pair<string, double*>(LOGARITHM, new double[11]{0,1,-0.5,1./3,-1./4,1./5,-1./6,1./7,-1./8,1./9,-1./10}));
	taylorCoeffsMap.insert(pair<string, double*>(EXPONENT, new double[11]{1,1,0.5,1./6,1./24,1./120,1./720,1./5040, 1./40320,1./362880,1./3628800}));
	taylorCoeffsMap.insert(pair<string, double*>(SIGMOID, new double[11]{1./2,1./4,0,-1./48,0,1./480,0,-17./80640,0,31./1451520,0}));
}

ZZX Context::encode(complex<double>* vals, long slots, long logp) {
	complex<double>* uvals = new complex<double>[slots];
	long i, jdx, idx;
	copy(vals, vals + slots, uvals);

	ZZX mx;
	mx.SetLength(N);
	long gap = Nh / slots;
	fftSpecialInv(uvals, slots);
	for (i = 0, jdx = Nh, idx = 0; i < slots; ++i, jdx += gap, idx += gap) {
		mx.rep[idx] = EvaluatorUtils::evalZZ(uvals[i].real(), logp);
		mx.rep[jdx] = EvaluatorUtils::evalZZ(uvals[i].imag(), logp);
	}
	delete[] uvals;
	return mx;
}

ZZX Context::encode(double* vals, long slots, long logp) {
	complex<double>* uvals = new complex<double>[slots];
	long i, jdx, idx;
	for (i = 0; i < slots; ++i) {
		uvals[i].real(vals[i]);
	}

	ZZX mx;
	mx.SetLength(N);

	long gap = Nh / slots;

	fftSpecialInv(uvals, slots);

	for (i = 0, jdx = Nh, idx = 0; i < slots; ++i, jdx += gap, idx += gap) {
		mx.rep[idx] = EvaluatorUtils::evalZZ(uvals[i].real(), logp);
		mx.rep[jdx] = EvaluatorUtils::evalZZ(uvals[i].imag(), logp);
	}
	delete[] uvals;
	return mx;
}

void Context::addBootContext(long logSlots, long logp) {
	if(bootContextMap.find(logSlots) == bootContextMap.end()) {
		long slots = 1 << logSlots;
		long dslots = slots << 1;
		long logk = logSlots >> 1;

		long k = 1 << logk;
		long i, pos,  ki, jdx, idx, deg;
		long gap = Nh >> logSlots;

		ZZX* pvecInv = new ZZX[slots];
		ZZX* pvec = new ZZX[slots];

		complex<double>* pvals = new complex<double>[dslots];

		ZZX p1, p2;
		double c = 0.25/M_PI;

		if(logSlots < logNh) {
			long dgap = gap >> 1;
			for (ki = 0; ki < slots; ki += k) {
				for (pos = ki; pos < ki + k; ++pos) {
					for (i = 0; i < slots - pos; ++i) {
						deg = ((M - rotGroup[i + pos]) * i * gap) % M;
						pvals[i].real(to_double(ksiPowsr[deg]));
						pvals[i].imag(to_double(ksiPowsi[deg]));
						pvals[i + slots].real(-pvals[i].imag());
						pvals[i + slots].imag(pvals[i].real());
					}
					for (i = slots - pos; i < slots; ++i) {
						deg =((M - rotGroup[i + pos - slots]) * i * gap) % M;
						pvals[i].real(to_double(ksiPowsr[deg]));
						pvals[i].imag(to_double(ksiPowsi[deg]));
						pvals[i + slots].real(-pvals[i].imag());
						pvals[i + slots].imag(pvals[i].real());
					}
					EvaluatorUtils::rightRotateAndEqual(pvals, dslots, ki);
					fftSpecialInv(pvals, dslots);
					pvec[pos].SetLength(N);
					for (i = 0, jdx = Nh, idx = 0; i < dslots; ++i, jdx += dgap, idx += dgap) {
						pvec[pos].rep[idx] = EvaluatorUtils::evalZZ(pvals[i].real(), logp);
						pvec[pos].rep[jdx] = EvaluatorUtils::evalZZ(pvals[i].imag(), logp);
					}
				}
			}
			for (i = 0; i < slots; ++i) {
				pvals[i] = 0.0;
				pvals[i + slots].real(0);
				pvals[i + slots].imag(-c);
			}
			p1.SetLength(N);
			fftSpecialInv(pvals, dslots);
			for (i = 0, jdx = Nh, idx = 0; i < dslots; ++i, jdx += dgap, idx += dgap) {
				p1.rep[idx] = EvaluatorUtils::evalZZ(pvals[i].real(), logp);
				p1.rep[jdx] = EvaluatorUtils::evalZZ(pvals[i].imag(), logp);
			}

			for (i = 0; i < slots; ++i) {
				pvals[i] = c;
				pvals[i + slots] = 0;
			}

			p2.SetLength(N);
			fftSpecialInv(pvals, dslots);
			for (i = 0, jdx = Nh, idx = 0; i < dslots; ++i, jdx += dgap, idx += dgap) {
				p2.rep[idx] = EvaluatorUtils::evalZZ(pvals[i].real(), logp);
				p2.rep[jdx] = EvaluatorUtils::evalZZ(pvals[i].imag(), logp);
			}
		} else {
			for (ki = 0; ki < slots; ki += k) {
				for (pos = ki; pos < ki + k; ++pos) {
					for (i = 0; i < slots - pos; ++i) {
						deg = ((M - rotGroup[i + pos]) * i * gap) % M;
						pvals[i].real(to_double(ksiPowsr[deg]));
						pvals[i].imag(to_double(ksiPowsi[deg]));
					}
					for (i = slots - pos; i < slots; ++i) {
						deg =((M - rotGroup[i + pos - slots]) * i * gap) % M;
						pvals[i].real(to_double(ksiPowsr[deg]));
						pvals[i].imag(to_double(ksiPowsi[deg]));
					}
					EvaluatorUtils::rightRotateAndEqual(pvals, slots, ki);
					fftSpecialInv(pvals, slots);
					pvec[pos].SetLength(N);
					for (i = 0, jdx = Nh, idx = 0; i < slots; ++i, jdx += gap, idx += gap) {
						pvec[pos].rep[idx] = EvaluatorUtils::evalZZ(pvals[i].real(), logp);
						pvec[pos].rep[jdx] = EvaluatorUtils::evalZZ(pvals[i].imag(), logp);
					}
				}
			}
		}

		for (ki = 0; ki < slots; ki += k) {
			for (pos = ki; pos < ki + k; ++pos) {

				for (i = 0; i < slots - pos; ++i) {
					deg = (rotGroup[i] * (i + pos) * gap) % M;
					pvals[i].real(to_double(ksiPowsr[deg]));
					pvals[i].imag(to_double(ksiPowsi[deg]));
				}
				for (i = slots - pos; i < slots; ++i) {
					deg = (rotGroup[i] * (i + pos - slots) * gap) % M;
					pvals[i].real(to_double(ksiPowsr[deg]));
					pvals[i].imag(to_double(ksiPowsi[deg]));
				}
				EvaluatorUtils::rightRotateAndEqual(pvals, slots, ki);
				fftSpecialInv(pvals, slots);
				pvecInv[pos].SetLength(N);
				for (i = 0, jdx = Nh, idx = 0; i < slots; ++i, jdx += gap, idx += gap) {
					pvecInv[pos].rep[idx] = EvaluatorUtils::evalZZ(pvals[i].real(), logp);
					pvecInv[pos].rep[jdx] = EvaluatorUtils::evalZZ(pvals[i].imag(), logp);
				}
			}
		}

		delete[] pvals;
		bootContextMap.insert(pair<long, BootContext>(logSlots, BootContext(pvec, pvecInv, p1, p2, logp)));
	}
}

void Context::bitReverse(complex<double>* vals, const long size) {
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

void Context::fft(complex<double>* vals, const long size) {
	bitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		long MoverLen = M / len;
		long lenh = len >> 1;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = j * MoverLen;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				RR tmp1 = to_RR(v.real()) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.real() + v.imag()) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.imag() - v.real()) * ksiPowsr[idx];
				v.real(to_double(tmpr));
				v.imag(to_double(tmpi));
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftInvLazy(complex<double>* vals, const long size) {
	bitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		long MoverLen = M / len;
		long lenh = len >> 1;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = (len - j) * MoverLen;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				RR tmp1 = to_RR(v.real()) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.real() + v.imag()) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.imag() - v.real()) * ksiPowsr[idx];
				v.real(to_double(tmpr));
				v.imag(to_double(tmpi));
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftInv(complex<double>* vals, const long size) {
	fftInvLazy(vals, size);
	for (long i = 0; i < size; ++i) {
		vals[i] /= size;
	}
}

void Context::fftSpecial(complex<double>* vals, const long size) {
	bitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		for (long i = 0; i < size; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			for (long j = 0; j < lenh; ++j) {
				long idx = ((rotGroup[j] % lenq)) * M / lenq;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				RR tmp1 = to_RR(v.real()) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.real() + v.imag()) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.imag() - v.real()) * ksiPowsr[idx];
				v.real(to_double(tmpr));
				v.imag(to_double(tmpi));
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftSpecialInvLazy(complex<double>* vals, const long size) {
	for (long len = size; len >= 1; len >>= 1) {
		for (long i = 0; i < size; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			for (long j = 0; j < lenh; ++j) {
				long idx = (lenq - (rotGroup[j] % lenq)) * M / lenq;
				complex<double> u = vals[i + j] + vals[i + j + lenh];
				complex<double> v = vals[i + j] - vals[i + j + lenh];
				RR tmp1 = to_RR(v.real()) * (ksiPowsr[idx] + ksiPowsi[idx]);
				RR tmpr = tmp1 - to_RR(v.real() + v.imag()) * ksiPowsi[idx];
				RR tmpi = tmp1 + to_RR(v.imag() - v.real()) * ksiPowsr[idx];
				v.real(to_double(tmpr));
				v.imag(to_double(tmpi));
				vals[i + j] = u;
				vals[i + j + lenh] = v;
			}
		}
	}
	bitReverse(vals, size);
}

void Context::fftSpecialInv(complex<double>* vals, const long size) {
	fftSpecialInvLazy(vals, size);
	for (long i = 0; i < size; ++i) {
		vals[i] /= size;
	}
}

Context::~Context() {
	delete[] rotGroup;
	delete[] ksiPowsr;
	delete[] ksiPowsi;
}

