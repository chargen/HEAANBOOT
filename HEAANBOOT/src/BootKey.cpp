#include "BootKey.h"

BootKey::BootKey(Context& context, long logp, long logl) : logp(logp) {
	ZZ p = power2_ZZ(logp);

	ZZX ex;
	long loglh = logl/2;
	long l = 1 << logl;
	long dl = l << 1;
	long Noverl = context.N >> logl;

	long lk = 1 << loglh;
	long lm = 1 << (logl - loglh);

	pvec = new ZZX[l];
	pvecInv = new ZZX[l];

	CZZ* pdvals = new CZZ[dl];
	for (long j = 0; j < l; ++j) {
		for (long i = j; i < l; ++i) {
			long deg =((2 * context.N - context.rotGroup[i]) * (i - j) * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], logp);
			long idx = (context.rotGroup[i-j] % (2 * dl)  - 1) / 2;
					pdvals[idx] = tmp;
					pdvals[dl - idx - 1] = tmp.conjugate();
		}
		for (long i = 0; i < j; ++i) {
			long deg =((2 * context.N - context.rotGroup[i]) * (l + i - j) * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], logp);
			long idx = (context.rotGroup[l + i - j] % (2 * dl) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[dl - idx - 1] = tmp.conjugate();
		}

		NumUtils::fftSpecialInv(pdvals, dl, context.ksiPowsr, context.ksiPowsi, context.M);

		pvec[j].SetLength(context.N);
		long idx = 0;
		long gap = context.N / dl;
		for (long i = 0; i < dl; ++i) {
			pvec[j].rep[idx] = pdvals[i].r;
			idx += gap;
		}
	}

	for (long j = 0; j < l; ++j) {
		for (long i = j; i < l; ++i) {
			long deg =(context.rotGroup[i-j] * i * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], logp);
			long idx = (context.rotGroup[i-j] % (2 * dl) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[dl - idx - 1] = tmp.conjugate();
		}
		for (long i = 0; i < j; ++i) {
			long deg = (context.rotGroup[l + i - j] * i * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], logp);
			long idx = (context.rotGroup[l + i - j] % (2 * dl) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[dl - idx - 1] = tmp.conjugate();
		}
		NumUtils::fftSpecialInv(pdvals, dl, context.ksiPowsr, context.ksiPowsi, context.M);

		pvecInv[j].SetLength(context.N);
		long idx = 0;
		long gap = context.N / dl;
		for (long i = 0; i < dl; ++i) {
			pvecInv[j].rep[idx] = pdvals[i].r;
			idx += gap;
		}
	}
	delete[] pdvals;

	for (long i = 1; i < lm; ++i) {
		for (long j = 0; j < lk; ++j) {
			pvec[j + lk * i] = Ring2Utils::inpower(pvec[j + lk * i], context.rotGroup[l - lk * i], p, context.N);
			pvecInv[j + lk * i] = Ring2Utils::inpower(pvecInv[j + lk * i], context.rotGroup[l - lk * i], p, context.N);
		}
	}
}
