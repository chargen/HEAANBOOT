#include "BootKey.h"

BootKey::BootKey(Context& context, long pBits, long l) : pBits(pBits) {
	ZZ pmod = power2_ZZ(pBits);

	ZZX ex;
	long lh = l/2;
	long lpow = 1 << l;
	long lpow2 = lpow << 1;
	long Noverl = context.N >> l;

	long lk = 1 << lh;
	long lm = 1 << (l - lh);

	pvec = new ZZX[lpow];
	pvecInv = new ZZX[lpow];

	CZZ* pdvals = new CZZ[lpow2];
	for (long j = 0; j < lpow; ++j) {
		for (long i = j; i < lpow; ++i) {
			long deg =((2 * context.N - context.rotGroup[i]) * (i - j) * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], pBits);
			long idx = (context.rotGroup[i-j] % (2 * lpow2)  - 1) / 2;
					pdvals[idx] = tmp;
					pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}
		for (long i = 0; i < j; ++i) {
			long deg =((2 * context.N - context.rotGroup[i]) * (lpow + i - j) * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], pBits);
			long idx = (context.rotGroup[lpow + i - j] % (2 * lpow2) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}

		NumUtils::fftSpecialInv(pdvals, lpow2, context.ksiPowsr, context.ksiPowsi, context.M);

		pvec[j].SetLength(context.N);
		long idx = 0;
		long gap = context.N / lpow2;
		for (long i = 0; i < lpow2; ++i) {
			pvec[j].rep[idx] = pdvals[i].r;
			idx += gap;
		}
	}

	for (long j = 0; j < lpow; ++j) {
		for (long i = j; i < lpow; ++i) {
			long deg =(context.rotGroup[i-j] * i * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], pBits);
			long idx = (context.rotGroup[i-j] % (2 * lpow2) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}
		for (long i = 0; i < j; ++i) {
			long deg = (context.rotGroup[lpow + i - j] * i * Noverl) % (2 * context.N);
			CZZ tmp = EvaluatorUtils::evalCZZ(context.ksiPowsr[deg], context.ksiPowsi[deg], pBits);
			long idx = (context.rotGroup[lpow + i - j] % (2 * lpow2) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}
		NumUtils::fftSpecialInv(pdvals, lpow2, context.ksiPowsr, context.ksiPowsi, context.M);

		pvecInv[j].SetLength(context.N);
		long idx = 0;
		long gap = context.N / lpow2;
		for (long i = 0; i < lpow2; ++i) {
			pvecInv[j].rep[idx] = pdvals[i].r;
			idx += gap;
		}
	}
	delete[] pdvals;

	for (long i = 1; i < lm; ++i) {
		for (long j = 0; j < lk; ++j) {
			pvec[j + lk * i] = Ring2Utils::inpower(pvec[j + lk * i], context.rotGroup[lpow - lk * i], pmod, context.N);
			pvecInv[j + lk * i] = Ring2Utils::inpower(pvecInv[j + lk * i], context.rotGroup[lpow - lk * i], pmod, context.N);
		}
	}
}
