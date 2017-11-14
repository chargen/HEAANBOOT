#include "PublicKey.h"

#include "NumUtils.h"
#include "Params.h"
#include "Ring2Utils.h"

PublicKey::PublicKey(Params& params, SecretKey& secretKey) {
	addEncKey(params, secretKey);
	addMultKey(params, secretKey);
}

void PublicKey::addEncKey(Params& params, SecretKey& secretKey) {
	ZZX ex, ax, bx;

	NumUtils::sampleUniform2(ax, params.N, params.logqq);
	NumUtils::sampleGauss(ex, params.N, params.sigma);
	Ring2Utils::mult(bx, secretKey.sx, ax, params.qq, params.N);
	Ring2Utils::sub(bx, ex, bx, params.qq, params.N);

	keyMap.insert(pair<long, RLWE>(ENCRYPTION, RLWE(ax, bx)));
}

void PublicKey::addMultKey(Params& params, SecretKey& secretKey) {
	ZZX ex, ax, bx, sxsx;

	Ring2Utils::mult(sxsx, secretKey.sx, secretKey.sx, params.q, params.N);
	Ring2Utils::leftShiftAndEqual(sxsx, params.logq, params.qq, params.N);
	NumUtils::sampleUniform2(ax, params.N, params.logqq);
	NumUtils::sampleGauss(ex, params.N, params.sigma);
	Ring2Utils::addAndEqual(ex, sxsx, params.qq, params.N);
	Ring2Utils::mult(bx, secretKey.sx, ax, params.qq, params.N);
	Ring2Utils::sub(bx, ex, bx, params.qq, params.N);

	keyMap.insert(pair<long, RLWE>(MULTIPLICATION, RLWE(ax, bx)));
}

void PublicKey::addConjKey(Params& params, SecretKey& secretKey) {
	ZZX ex, ax, bx, sxconj;
	Ring2Utils::conjugate(sxconj, secretKey.sx, params.N);
	Ring2Utils::leftShiftAndEqual(sxconj, params.logq, params.qq, params.N);
	NumUtils::sampleUniform2(ax, params.N, params.logqq);
	NumUtils::sampleGauss(ex, params.N, params.sigma);
	Ring2Utils::addAndEqual(ex, sxconj, params.qq, params.N);
	Ring2Utils::mult(bx, secretKey.sx, ax, params.qq, params.N);
	Ring2Utils::sub(bx, ex, bx, params.qq, params.N);

	keyMap.insert(pair<long, RLWE>(CONJUGATION, RLWE(ax, bx)));
}

void PublicKey::addLeftRotKey(Params& params, SecretKey& secretKey, long rot) {
	ZZX ex, sxsx, ax, bx, spow;
	Ring2Utils::inpower(spow, secretKey.sx, params.rotGroup[rot], params.q, params.N);
	Ring2Utils::leftShiftAndEqual(spow, params.logq, params.qq, params.N);
	NumUtils::sampleUniform2(ax, params.N, params.logqq);
	NumUtils::sampleGauss(ex, params.N, params.sigma);
	Ring2Utils::addAndEqual(ex, spow, params.qq, params.N);
	Ring2Utils::mult(bx, secretKey.sx, ax, params.qq, params.N);
	Ring2Utils::sub(bx, ex, bx, params.qq, params.N);
	leftRotKeyMap.insert(pair<long, RLWE>(rot, RLWE(ax, bx)));
}

void PublicKey::addLeftRotKeys(Params& params, SecretKey& secretKey) {
	for (long i = 0; i < params.logN - 1; ++i) {
		long idx = 1 << i;
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(params, secretKey, idx);
		}
	}
}

void PublicKey::addRightRotKeys(Params& params, SecretKey& secretKey) {
	for (long i = 0; i < params.logN - 1; ++i) {
		long idx = params.N/2 - (1 << i);
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(params, secretKey, idx);
		}
	}
}

void PublicKey::addBootKeys(Params& params, SecretKey& secretKey, SchemeAux& aux, long lkey, long pBits) {

	if(bootKeyMap.find(lkey) == bootKeyMap.end()) {
		bootKeyMap.insert(pair<long, BootKey>(lkey, BootKey(params, aux, pBits, lkey)));
	}

	long lkeyh = lkey/2;
	long k = 1 << lkeyh;
	long m = 1 << (lkey - lkeyh);

	for (long i = 1; i < k; ++i) {
		if(leftRotKeyMap.find(i) == leftRotKeyMap.end()) {
			addLeftRotKey(params, secretKey, i);
		}
	}

	for (long i = 1; i < m; ++i) {
		long idx = i * k;
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(params, secretKey, idx);
		}
	}
}

void PublicKey::addSortKeys(Params& params, SecretKey& secretKey, long size) {
	for (long i = 1; i < size; ++i) {
		if(leftRotKeyMap.find(i) == leftRotKeyMap.end()) {
			addLeftRotKey(params, secretKey, i);
		}
	}
}

