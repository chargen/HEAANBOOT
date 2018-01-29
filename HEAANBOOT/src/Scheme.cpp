#include "Scheme.h"

#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include "EvaluatorUtils.h"
#include "NumUtils.h"
#include "Ring2Utils.h"
#include "StringUtils.h"

//-----------------------------------------

Scheme::Scheme(SecretKey& secretKey, Context& context) : context(context) {
	addEncKey(secretKey);
	addMultKey(secretKey);
};

void Scheme::addEncKey(SecretKey& secretKey) {
	ZZX ex, ax, bx;

	NumUtils::sampleUniform2(ax, context.N, context.logPQ);
	NumUtils::sampleGauss(ex, context.N, context.sigma);
	Ring2Utils::mult(bx, secretKey.sx, ax, context.PQ, context.N);
	Ring2Utils::sub(bx, ex, bx, context.PQ, context.N);

	keyMap.insert(pair<long, Key>(ENCRYPTION, Key(ax, bx)));
}

void Scheme::addMultKey(SecretKey& secretKey) {
	ZZX ex, ax, bx, sxsx;

	Ring2Utils::mult(sxsx, secretKey.sx, secretKey.sx, context.Q, context.N);
	Ring2Utils::leftShiftAndEqual(sxsx, context.logQ, context.PQ, context.N);
	NumUtils::sampleUniform2(ax, context.N, context.logPQ);
	NumUtils::sampleGauss(ex, context.N, context.sigma);
	Ring2Utils::addAndEqual(ex, sxsx, context.PQ, context.N);
	Ring2Utils::mult(bx, secretKey.sx, ax, context.PQ, context.N);
	Ring2Utils::sub(bx, ex, bx, context.PQ, context.N);

	keyMap.insert(pair<long, Key>(MULTIPLICATION, Key(ax, bx)));
}

void Scheme::addConjKey(SecretKey& secretKey) {
	ZZX ex, ax, bx, sxconj;

	Ring2Utils::conjugate(sxconj, secretKey.sx, context.N);
	Ring2Utils::leftShiftAndEqual(sxconj, context.logQ, context.PQ, context.N);
	NumUtils::sampleUniform2(ax, context.N, context.logPQ);
	NumUtils::sampleGauss(ex, context.N, context.sigma);
	Ring2Utils::addAndEqual(ex, sxconj, context.PQ, context.N);
	Ring2Utils::mult(bx, secretKey.sx, ax, context.PQ, context.N);
	Ring2Utils::sub(bx, ex, bx, context.PQ, context.N);

	keyMap.insert(pair<long, Key>(CONJUGATION, Key(ax, bx)));
}

void Scheme::addLeftRotKey(SecretKey& secretKey, long rot) {
	ZZX ex, ax, bx, sxrot;

	Ring2Utils::inpower(sxrot, secretKey.sx, context.rotGroup[rot], context.Q, context.N);
	Ring2Utils::leftShiftAndEqual(sxrot, context.logQ, context.PQ, context.N);
	NumUtils::sampleUniform2(ax, context.N, context.logPQ);
	NumUtils::sampleGauss(ex, context.N, context.sigma);
	Ring2Utils::addAndEqual(ex, sxrot, context.PQ, context.N);
	Ring2Utils::mult(bx, secretKey.sx, ax, context.PQ, context.N);
	Ring2Utils::sub(bx, ex, bx, context.PQ, context.N);

	leftRotKeyMap.insert(pair<long, Key>(rot, Key(ax, bx)));
}

void Scheme::addLeftRotKeys(SecretKey& secretKey) {
	for (long i = 0; i < context.logNh; ++i) {
		long idx = 1 << i;
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx);
		}
	}
}

void Scheme::addRightRotKeys(SecretKey& secretKey) {
	for (long i = 0; i < context.logNh; ++i) {
		long idx = context.N/2 - (1 << i);
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx);
		}
	}
}

void Scheme::addBootKey(SecretKey& secretKey, long logSlots, long logp) {
	context.addBootContext(logSlots, logp);

	addConjKey(secretKey);
	addLeftRotKeys(secretKey);

	long logk = logSlots / 2;
	long k = 1 << logk;
	long m = 1 << (logSlots - logk);

	for (long i = 1; i < k; ++i) {
		if(leftRotKeyMap.find(i) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, i);
		}
	}

	for (long i = 1; i < m; ++i) {
		long idx = i * k;
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx);
		}
	}
}

void Scheme::addSortKeys(SecretKey& secretKey, long size) {
	for (long i = 1; i < size; ++i) {
		if(leftRotKeyMap.find(i) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, i);
		}
	}
}

Plaintext Scheme::encode(CZZ* vals, long slots, long logq, bool isComplex) {
	ZZX mx = context.encodeLarge(vals, slots);
	ZZ q = power2_ZZ(logq);
	return Plaintext(mx, q, logq, slots, isComplex);
}

CZZ* Scheme::decode(Plaintext& msg) {
	long i, idx;
	long slots = msg.slots;
	long gap = context.Nh / slots;
	CZZ* res = new CZZ[slots];

	ZZ tmp;
	for (i = 0, idx = 0; i < slots; ++i, idx += gap) {
		rem(tmp, msg.mx[idx], msg.q);
		if(NumBits(tmp) == msg.logq) tmp -= msg.q;
		res[i].r = tmp;
		rem(tmp, msg.mx[idx + context.Nh], msg.q);
		if(NumBits(tmp) == msg.logq) tmp -= msg.q;
		res[i].i = tmp;
	}
	context.fftSpecial(res, slots);
	return res;
}

Plaintext Scheme::encodeSingle(CZZ& val, long logq, bool isComplex) {
	ZZX mx;
	mx.SetLength(context.N);

	ZZ q = power2_ZZ(logq);
	mx.rep[0] = val.r << context.logQ;
	if(isComplex) {
		mx.rep[context.Nh] = val.i << context.logQ;
	}
	return Plaintext(mx, q, logq, 1, isComplex);
}

CZZ Scheme::decodeSingle(Plaintext& msg) {
	CZZ res;
	ZZ tmp = msg.mx.rep[0] % msg.q;
	if(NumBits(tmp) == msg.logq) tmp -= msg.q;
	res.r = tmp;

	if(msg.isComplex) {
		tmp = msg.mx.rep[context.Nh] % msg.q;
		if(NumBits(tmp) == msg.logq) tmp -= msg.q;
		res.i = tmp;
	}
	return res;
}

Ciphertext Scheme::encryptMsg(Plaintext& msg) {
	ZZX ax, bx, vx, ex;
	Key key = keyMap.at(ENCRYPTION);
	ZZ qQ = msg.q << context.logQ;

	NumUtils::sampleZO(vx, context.N);
	Ring2Utils::mult(ax, vx, key.ax, qQ, context.N);
	NumUtils::sampleGauss(ex, context.N, context.sigma);
	Ring2Utils::addAndEqual(ax, ex, qQ, context.N);

	Ring2Utils::mult(bx, vx, key.bx, qQ, context.N);
	NumUtils::sampleGauss(ex, context.N, context.sigma);
	Ring2Utils::addAndEqual(bx, ex, qQ, context.N);

	Ring2Utils::addAndEqual(bx, msg.mx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(ax, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(bx, context.logQ, context.N);

	return Ciphertext(ax, bx, msg.q, msg.logq, msg.slots, msg.isComplex);
}

Plaintext Scheme::decryptMsg(SecretKey& secretKey, Ciphertext& cipher) {
	ZZX mx;

	Ring2Utils::mult(mx, cipher.ax, secretKey.sx, cipher.q, context.N);
	Ring2Utils::addAndEqual(mx, cipher.bx, cipher.q, context.N);

	return Plaintext(mx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

Ciphertext Scheme::encrypt(CZZ* vals, long slots, long logq, bool isComplex) {
	Plaintext msg = encode(vals, slots, logq, isComplex);
	return encryptMsg(msg);
}

CZZ* Scheme::decrypt(SecretKey& secretKey, Ciphertext& cipher) {
	Plaintext msg = decryptMsg(secretKey, cipher);
	return decode(msg);
}

Ciphertext Scheme::encryptSingle(CZZ& val, long logq, bool isComplex) {
	Plaintext msg = encodeSingle(val, logq, isComplex);
	return encryptMsg(msg);
}

CZZ Scheme::decryptSingle(SecretKey& secretKey, Ciphertext& cipher) {
	Plaintext msg = decryptMsg(secretKey, cipher);
	return decodeSingle(msg);
}

//-----------------------------------------

Ciphertext Scheme::add(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZX ax, bx;

	Ring2Utils::add(ax, cipher1.ax, cipher2.ax, cipher1.q, context.N);
	Ring2Utils::add(bx, cipher1.bx, cipher2.bx, cipher1.q, context.N);

	return Ciphertext(ax, bx, cipher1.q, cipher1.logq, cipher1.slots, cipher1.isComplex);
}

void Scheme::addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	Ring2Utils::addAndEqual(cipher1.ax, cipher2.ax, cipher1.q, context.N);
	Ring2Utils::addAndEqual(cipher1.bx, cipher2.bx, cipher1.q, context.N);
}

//-----------------------------------------

Ciphertext Scheme::addConst(Ciphertext& cipher, ZZ& cnst) {
	ZZX ax = cipher.ax;
	ZZX bx = cipher.bx;

	AddMod(bx.rep[0], cipher.bx.rep[0], cnst, cipher.q);
	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

Ciphertext Scheme::addConst(Ciphertext& cipher, CZZ& cnst) {
	ZZX ax = cipher.ax;
	ZZX bx = cipher.bx;

	AddMod(bx.rep[0], cipher.bx.rep[0], cnst.r, cipher.q);
	AddMod(bx.rep[context.Nh], cipher.bx.rep[context.Nh], cnst.i, cipher.q);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::addConstAndEqual(Ciphertext& cipher, ZZ& cnst) {
	AddMod(cipher.bx.rep[0], cipher.bx.rep[0], cnst, cipher.q);
}


void Scheme::addConstAndEqual(Ciphertext& cipher, CZZ& cnst) {
	AddMod(cipher.bx.rep[0], cipher.bx.rep[0], cnst.r, cipher.q);
	AddMod(cipher.bx.rep[context.Nh], cipher.bx.rep[context.Nh], cnst.i, cipher.q);
}

//-----------------------------------------

Ciphertext Scheme::sub(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZX ax, bx;

	Ring2Utils::sub(ax, cipher1.ax, cipher2.ax, cipher1.q, context.N);
	Ring2Utils::sub(bx, cipher1.bx, cipher2.bx, cipher1.q, context.N);

	return Ciphertext(ax, bx, cipher1.q, cipher1.logq, cipher1.slots, cipher1.isComplex);
}

void Scheme::subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	Ring2Utils::subAndEqual(cipher1.ax, cipher2.ax, cipher1.q, context.N);
	Ring2Utils::subAndEqual(cipher1.bx, cipher2.bx, cipher1.q, context.N);
}

void Scheme::subAndEqual2(Ciphertext& cipher1, Ciphertext& cipher2) {
	Ring2Utils::subAndEqual2(cipher1.ax, cipher2.ax, cipher1.q, context.N);
	Ring2Utils::subAndEqual2(cipher1.bx, cipher2.bx, cipher1.q, context.N);
}

Ciphertext Scheme::imult(Ciphertext& cipher) {
	ZZX ax, bx;

	Ring2Utils::multByMonomial(ax, cipher.ax, context.Nh, context.N);
	Ring2Utils::multByMonomial(bx, cipher.bx, context.Nh, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

Ciphertext Scheme::idiv(Ciphertext& cipher) {
	ZZX ax, bx;

	Ring2Utils::multByMonomial(ax, cipher.ax, 3 * context.Nh, context.N);
	Ring2Utils::multByMonomial(bx, cipher.bx, 3 * context.Nh, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::imultAndEqual(Ciphertext& cipher) {
	Ring2Utils::multByMonomialAndEqual(cipher.ax, context.Nh, context.N);
	Ring2Utils::multByMonomialAndEqual(cipher.bx, context.Nh, context.N);
}

void Scheme::idivAndEqual(Ciphertext& cipher) {
	Ring2Utils::multByMonomialAndEqual(cipher.ax, 3 * context.Nh, context.N);
	Ring2Utils::multByMonomialAndEqual(cipher.bx, 3 * context.Nh, context.N);
}

Ciphertext Scheme::mult(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ qQ = cipher1.q << context.logQ;
	ZZX axbx1, axbx2, axax, bxbx, axmult, bxmult;
	Key key = keyMap.at(MULTIPLICATION);

	Ring2Utils::add(axbx1, cipher1.ax, cipher1.bx, cipher1.q, context.N);
	Ring2Utils::add(axbx2, cipher2.ax, cipher2.bx, cipher1.q, context.N);
	Ring2Utils::multAndEqual(axbx1, axbx2, cipher1.q, context.N);

	Ring2Utils::mult(axax, cipher1.ax, cipher2.ax, cipher1.q, context.N);
	Ring2Utils::mult(bxbx, cipher1.bx, cipher2.bx, cipher1.q, context.N);

	Ring2Utils::mult(axmult, axax, key.ax, qQ, context.N);
	Ring2Utils::mult(bxmult, axax, key.bx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(axmult, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(bxmult, context.logQ, context.N);

	Ring2Utils::addAndEqual(axmult, axbx1, cipher1.q, context.N);
	Ring2Utils::subAndEqual(axmult, bxbx, cipher1.q, context.N);
	Ring2Utils::subAndEqual(axmult, axax, cipher1.q, context.N);
	Ring2Utils::addAndEqual(bxmult, bxbx, cipher1.q, context.N);

	return Ciphertext(axmult, bxmult, cipher1.q, cipher1.logq, cipher1.slots, cipher1.isComplex);
}

void Scheme::multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	ZZ qQ = cipher1.q << context.logQ;
	ZZX axbx1, axbx2, axax, bxbx;
	Key key = keyMap.at(MULTIPLICATION);

	Ring2Utils::add(axbx1, cipher1.ax, cipher1.bx, cipher1.q, context.N);
	Ring2Utils::add(axbx2, cipher2.ax, cipher2.bx, cipher1.q, context.N);
	Ring2Utils::multAndEqual(axbx1, axbx2, cipher1.q, context.N);

	Ring2Utils::mult(axax, cipher1.ax, cipher2.ax, cipher1.q, context.N);
	Ring2Utils::mult(bxbx, cipher1.bx, cipher2.bx, cipher1.q, context.N);

	Ring2Utils::mult(cipher1.ax, axax, key.ax, qQ, context.N);
	Ring2Utils::mult(cipher1.bx, axax, key.bx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(cipher1.ax, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(cipher1.bx, context.logQ, context.N);

	Ring2Utils::addAndEqual(cipher1.ax, axbx1, cipher1.q, context.N);
	Ring2Utils::subAndEqual(cipher1.ax, bxbx, cipher1.q, context.N);
	Ring2Utils::subAndEqual(cipher1.ax, axax, cipher1.q, context.N);
	Ring2Utils::addAndEqual(cipher1.bx, bxbx, cipher1.q, context.N);
}

//-----------------------------------------

Ciphertext Scheme::square(Ciphertext& cipher) {
	ZZ qQ = cipher.q << context.logQ;
	ZZX axax, axbx, bxbx, bxmult, axmult;
	Key key = keyMap.at(MULTIPLICATION);

	Ring2Utils::square(bxbx, cipher.bx, cipher.q, context.N);
	Ring2Utils::mult(axbx, cipher.ax, cipher.bx, cipher.q, context.N);
	Ring2Utils::addAndEqual(axbx, axbx, cipher.q, context.N);
	Ring2Utils::square(axax, cipher.ax, cipher.q, context.N);

	Ring2Utils::mult(axmult, axax, key.ax, qQ, context.N);
	Ring2Utils::mult(bxmult, axax, key.bx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(axmult, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(bxmult, context.logQ, context.N);

	Ring2Utils::addAndEqual(axmult, axbx, cipher.q, context.N);
	Ring2Utils::addAndEqual(bxmult, bxbx, cipher.q, context.N);

	return Ciphertext(axmult, bxmult, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::squareAndEqual(Ciphertext& cipher) {
	ZZ qQ = cipher.q << context.logQ;
	ZZX bxbx, axbx, axax;
	Key key = keyMap.at(MULTIPLICATION);

	Ring2Utils::square(bxbx, cipher.bx, cipher.q, context.N);
	Ring2Utils::mult(axbx, cipher.bx, cipher.ax, cipher.q, context.N);
	Ring2Utils::addAndEqual(axbx, axbx, cipher.q, context.N);
	Ring2Utils::square(axax, cipher.ax, cipher.q, context.N);

	Ring2Utils::mult(cipher.ax, axax, key.ax, qQ, context.N);
	Ring2Utils::mult(cipher.bx, axax, key.bx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(cipher.ax, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(cipher.bx, context.logQ, context.N);

	Ring2Utils::addAndEqual(cipher.ax, axbx, cipher.q, context.N);
	Ring2Utils::addAndEqual(cipher.bx, bxbx, cipher.q, context.N);
}

//-----------------------------------------

Ciphertext Scheme::multByConst(Ciphertext& cipher, ZZ& cnst) {
	ZZX ax, bx;

	Ring2Utils::multByConst(ax, cipher.ax, cnst, cipher.q, context.N);
	Ring2Utils::multByConst(bx, cipher.bx, cnst, cipher.q, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

Ciphertext Scheme::multByConst(Ciphertext& cipher, CZZ& cnst) {
	ZZX axr, bxr, axi, bxi;

	Ring2Utils::multByMonomial(axi, cipher.ax, context.Nh, context.N);
	Ring2Utils::multByMonomial(bxi, cipher.bx, context.Nh, context.N);

	Ring2Utils::multByConst(axr, cipher.ax, cnst.r, cipher.q, context.N);
	Ring2Utils::multByConst(bxr, cipher.bx, cnst.r, cipher.q, context.N);

	Ring2Utils::multByConstAndEqual(axi, cnst.i, cipher.q, context.N);
	Ring2Utils::multByConstAndEqual(bxi, cnst.i, cipher.q, context.N);

	Ring2Utils::addAndEqual(axr, axi, cipher.q, context.N);
	Ring2Utils::addAndEqual(bxr, bxi, cipher.q, context.N);

	return Ciphertext(axr, bxr, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::multByConstAndEqual(Ciphertext& cipher, ZZ& cnst) {
	Ring2Utils::multByConstAndEqual(cipher.ax, cnst, cipher.q, context.N);
	Ring2Utils::multByConstAndEqual(cipher.bx, cnst, cipher.q, context.N);
}

void Scheme::multByConstAndEqual(Ciphertext& cipher, CZZ& cnst) {
	ZZX axi, bxi;

	Ring2Utils::multByMonomial(axi, cipher.ax, context.Nh, context.N);
	Ring2Utils::multByMonomial(bxi, cipher.bx, context.Nh, context.N);

	Ring2Utils::multByConstAndEqual(cipher.ax, cnst.r, cipher.q, context.N);
	Ring2Utils::multByConstAndEqual(cipher.bx, cnst.r, cipher.q, context.N);

	Ring2Utils::multByConstAndEqual(axi, cnst.i, cipher.q, context.N);
	Ring2Utils::multByConstAndEqual(bxi, cnst.i, cipher.q, context.N);

	Ring2Utils::addAndEqual(cipher.ax, axi, cipher.q, context.N);
	Ring2Utils::addAndEqual(cipher.bx, bxi, cipher.q, context.N);
}

Ciphertext Scheme::multByPoly(Ciphertext& cipher, ZZX& poly) {
	ZZX ax, bx;

	Ring2Utils::mult(ax, cipher.ax, poly, cipher.q, context.N);
	Ring2Utils::mult(bx, cipher.bx, poly, cipher.q, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::multByPolyAndEqual(Ciphertext& cipher, ZZX& poly) {
	Ring2Utils::multAndEqual(cipher.ax, poly, cipher.q, context.N);
	Ring2Utils::multAndEqual(cipher.bx, poly, cipher.q, context.N);
}

//-----------------------------------------

Ciphertext Scheme::multByMonomial(Ciphertext& cipher, const long degree) {
	ZZX ax, bx;

	Ring2Utils::multByMonomial(ax, cipher.ax, degree, context.N);
	Ring2Utils::multByMonomial(bx, cipher.bx, degree, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::multByMonomialAndEqual(Ciphertext& cipher, const long degree) {
	Ring2Utils::multByMonomialAndEqual(cipher.ax, degree, context.N);
	Ring2Utils::multByMonomialAndEqual(cipher.bx, degree, context.N);
}

//-----------------------------------------

Ciphertext Scheme::leftShift(Ciphertext& cipher, long bits) {
	ZZX ax, bx;

	Ring2Utils::leftShift(ax, cipher.ax, bits, cipher.q, context.N);
	Ring2Utils::leftShift(bx, cipher.bx, bits, cipher.q, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::leftShiftAndEqual(Ciphertext& cipher, long bits) {
	Ring2Utils::leftShiftAndEqual(cipher.ax, bits, cipher.q, context.N);
	Ring2Utils::leftShiftAndEqual(cipher.bx, bits, cipher.q, context.N);
}

void Scheme::doubleAndEqual(Ciphertext& cipher) {
	Ring2Utils::doubleAndEqual(cipher.ax, cipher.q, context.N);
	Ring2Utils::doubleAndEqual(cipher.bx, cipher.q, context.N);
}

//-----------------------------------------

Ciphertext Scheme::reScaleBy(Ciphertext& cipher, long bitsDown) {
	ZZX ax, bx;
	long newlogq = cipher.logq - bitsDown;
	ZZ newq = cipher.q >> bitsDown;

	Ring2Utils::rightShift(ax, cipher.ax, bitsDown, context.N);
	Ring2Utils::rightShift(bx, cipher.bx, bitsDown, context.N);

	return Ciphertext(ax, bx, newq, newlogq, cipher.slots, cipher.isComplex);
}

Ciphertext Scheme::reScaleTo(Ciphertext& cipher, long newlogq) {
	ZZX ax, bx;
	long bitsDown = cipher.logq - newlogq;
	ZZ newq = power2_ZZ(newlogq);

	Ring2Utils::rightShift(ax, cipher.ax, bitsDown, context.N);
	Ring2Utils::rightShift(bx, cipher.bx, bitsDown, context.N);

	return Ciphertext(ax, bx, newq, newlogq, cipher.slots, cipher.isComplex);
}

void Scheme::reScaleByAndEqual(Ciphertext& cipher, long bitsDown) {
	Ring2Utils::rightShiftAndEqual(cipher.ax, bitsDown, context.N);
	Ring2Utils::rightShiftAndEqual(cipher.bx, bitsDown, context.N);

	cipher.logq -= bitsDown;
	cipher.q >>= bitsDown;
}

void Scheme::reScaleToAndEqual(Ciphertext& cipher, long newlogq) {
	long bitsDown = cipher.logq - newlogq;
	cipher.logq = newlogq;
	cipher.q = power2_ZZ(newlogq);

	Ring2Utils::rightShiftAndEqual(cipher.ax, bitsDown, context.N);
	Ring2Utils::rightShiftAndEqual(cipher.bx, bitsDown, context.N);
}

Ciphertext Scheme::modDownBy(Ciphertext& cipher, long bitsDown) {
	ZZX bx, ax;
	ZZ newq = cipher.q >> bitsDown;
	long newlogq = cipher.logq - bitsDown;

	Ring2Utils::mod(ax, cipher.ax, newq, context.N);
	Ring2Utils::mod(bx, cipher.bx, newq, context.N);

	return Ciphertext(ax, bx, newq, newlogq, cipher.slots, cipher.isComplex);
}

void Scheme::modDownByAndEqual(Ciphertext& cipher, long bitsDown) {
	cipher.q >>= bitsDown;
	cipher.logq -= bitsDown;

	Ring2Utils::modAndEqual(cipher.ax, cipher.q, context.N);
	Ring2Utils::modAndEqual(cipher.bx, cipher.q, context.N);
}

Ciphertext Scheme::modDownTo(Ciphertext& cipher, long newlogq) {
	ZZX bx, ax;
	ZZ newq = power2_ZZ(newlogq);

	Ring2Utils::mod(ax, cipher.ax, newq, context.N);
	Ring2Utils::mod(bx, cipher.bx, newq, context.N);
	return Ciphertext(ax, bx, newq, newlogq, cipher.slots);
}

void Scheme::modDownToAndEqual(Ciphertext& cipher, long newlogq) {
	if(cipher.logq < newlogq) {
		invalid_argument("logq of cipher should be larger than newlogq");
	}
	if(cipher.logq == newlogq) {
		return;
	}

	cipher.q = power2_ZZ(newlogq);
	cipher.logq = newlogq;

	Ring2Utils::modAndEqual(cipher.ax, cipher.q, context.N);
	Ring2Utils::modAndEqual(cipher.bx, cipher.q, context.N);
}

Ciphertext Scheme::leftRotateFast(Ciphertext& cipher, long rotSlots) {
	ZZ qQ = cipher.q << context.logQ;
	ZZX bxrot, ax, bx;
	Key key = leftRotKeyMap.at(rotSlots);

	Ring2Utils::inpower(bxrot, cipher.bx, context.rotGroup[rotSlots], context.Q, context.N);
	Ring2Utils::inpower(bx, cipher.ax, context.rotGroup[rotSlots], context.Q, context.N);

	Ring2Utils::mult(ax, bx, key.ax, qQ, context.N);
	Ring2Utils::multAndEqual(bx, key.bx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(ax, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(bx, context.logQ, context.N);

	Ring2Utils::addAndEqual(bx, bxrot, cipher.q, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::leftRotateAndEqualFast(Ciphertext& cipher, long rotSlots) {
	ZZ qQ = cipher.q << context.logQ;
	ZZX bxrot;
	Key key = leftRotKeyMap.at(rotSlots);

	Ring2Utils::inpower(bxrot, cipher.bx, context.rotGroup[rotSlots], context.Q, context.N);
	Ring2Utils::inpower(cipher.bx, cipher.ax, context.rotGroup[rotSlots], context.Q, context.N);

	Ring2Utils::mult(cipher.ax, cipher.bx, key.ax, qQ, context.N);
	Ring2Utils::multAndEqual(cipher.bx, key.bx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(cipher.ax, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(cipher.bx, context.logQ, context.N);

	Ring2Utils::addAndEqual(cipher.bx, bxrot, cipher.q, context.N);
}

Ciphertext Scheme::leftRotateByPo2(Ciphertext& cipher, long logrotSlots) {
	long rotSlots = (1 << logrotSlots);
	return leftRotateFast(cipher, rotSlots);
}

void Scheme::leftRotateByPo2AndEqual(Ciphertext& cipher, long logrotSlots) {
	long rotSlots = (1 << logrotSlots);
	leftRotateAndEqualFast(cipher, rotSlots);
}

Ciphertext Scheme::rightRotateByPo2(Ciphertext& cipher, long logrotSlots) {
	long rotSlots = context.Nh - (1 << logrotSlots);
	return leftRotateFast(cipher, rotSlots);
}

void Scheme::rightRotateByPo2AndEqual(Ciphertext& cipher, long logrotSlots) {
	long rotSlots = context.Nh - (1 << logrotSlots);
	leftRotateAndEqualFast(cipher, rotSlots);
}

Ciphertext Scheme::leftRotate(Ciphertext& cipher, long rotSlots) {
	Ciphertext res = cipher;
	leftRotateAndEqual(res, rotSlots);
	return res;
}

void Scheme::leftRotateAndEqual(Ciphertext& cipher, long rotSlots) {
	long remrotSlots = rotSlots % cipher.slots;
	long logrotSlots = log2((double)remrotSlots) + 1;
	for (long i = 0; i < logrotSlots; ++i) {
		if(bit(remrotSlots, i)) {
			leftRotateByPo2AndEqual(cipher, i);
		}
	}
}

Ciphertext Scheme::rightRotate(Ciphertext& cipher, long rotSlots) {
	Ciphertext res = cipher;
	rightRotateAndEqual(res, rotSlots);
	return res;
}

void Scheme::rightRotateAndEqual(Ciphertext& cipher, long rotSlots) {
	long remrotSlots = rotSlots % cipher.slots;
	long logrotSlots = log2((double)remrotSlots) + 1;
	for (long i = 0; i < logrotSlots; ++i) {
		if(bit(remrotSlots, i)) {
			rightRotateByPo2AndEqual(cipher, i);
		}
	}
}

Ciphertext Scheme::conjugate(Ciphertext& cipher) {
	ZZ qQ = cipher.q << context.logQ;
	ZZX bxconj, ax, bx;
	Key key = keyMap.at(CONJUGATION);

	Ring2Utils::conjugate(bxconj, cipher.bx, context.N);
	Ring2Utils::conjugate(bx, cipher.ax, context.N);

	Ring2Utils::mult(ax, bx, key.ax, qQ, context.N);
	Ring2Utils::multAndEqual(bx, key.bx, qQ, context.N);

	Ring2Utils::rightShiftAndEqual(ax, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(bx, context.logQ, context.N);

	Ring2Utils::addAndEqual(bx, bxconj, cipher.q, context.N);

	return Ciphertext(ax, bx, cipher.q, cipher.logq, cipher.slots, cipher.isComplex);
}

void Scheme::conjugateAndEqual(Ciphertext& cipher) {
	ZZ Pq = cipher.q << context.logQ;
	ZZX bxconj;
	Key key = keyMap.at(CONJUGATION);

	Ring2Utils::conjugate(bxconj, cipher.bx, context.N);
	Ring2Utils::conjugate(cipher.bx, cipher.ax, context.N);

	Ring2Utils::mult(cipher.ax, cipher.bx, key.ax, Pq, context.N);
	Ring2Utils::multAndEqual(cipher.bx, key.bx, Pq, context.N);

	Ring2Utils::rightShiftAndEqual(cipher.ax, context.logQ, context.N);
	Ring2Utils::rightShiftAndEqual(cipher.bx, context.logQ, context.N);

	Ring2Utils::addAndEqual(cipher.bx, bxconj, cipher.q, context.N);
}

void Scheme::normalizeAndEqual(Ciphertext& cipher) {
	for (int i = 0; i < context.N; ++i) {
		if(NumBits(cipher.ax.rep[i]) == cipher.logq) cipher.ax.rep[i] -= cipher.q;
		if(NumBits(cipher.bx.rep[i]) == cipher.logq) cipher.bx.rep[i] -= cipher.q;
	}
}

void Scheme::linTransformAndEqual(Ciphertext& cipher) {
	long slots = cipher.slots;
	long logSlots = log2(slots);
	long logk = logSlots / 2;
	long k = 1 << logk;

	long ki, j;
	Ciphertext* rotvec = new Ciphertext[k];
	rotvec[0] = cipher;

	NTL_EXEC_RANGE(k - 1, first, last);
	for (j = first; j < last; ++j) {
		rotvec[j + 1] = leftRotateFast(rotvec[0], j + 1);
	}
	NTL_EXEC_RANGE_END;

	BootContext bootContext = context.bootContextMap.at(logSlots);

	Ciphertext* tmpvec = new Ciphertext[k];

	NTL_EXEC_RANGE(k, first, last);
	for (j = first; j < last; ++j) {
		tmpvec[j] = multByPoly(rotvec[j], bootContext.pvec[j]);
	}
	NTL_EXEC_RANGE_END;

	for (j = 1; j < k; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	cipher = tmpvec[0];

	for (ki = k; ki < slots; ki += k) {
		NTL_EXEC_RANGE(k, first, last);
		for (j = first; j < last; ++j) {
			tmpvec[j] = multByPoly(rotvec[j], bootContext.pvec[j + ki]);
		}
		NTL_EXEC_RANGE_END;
		for (j = 1; j < k; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}
		leftRotateAndEqualFast(tmpvec[0], ki);
		addAndEqual(cipher, tmpvec[0]);
	}
	delete[] rotvec;
	delete[] tmpvec;
}

void Scheme::linTransformInvAndEqual(Ciphertext& cipher) {
	long slots = cipher.slots;
	long logSlots = log2(slots);
	long logk = logSlots / 2;
	long k = 1 << logk;

	long ki, j;
	Ciphertext* rotvec = new Ciphertext[k];
	rotvec[0] = cipher;

	NTL_EXEC_RANGE(k-1, first, last);
	for (j = first; j < last; ++j) {
		rotvec[j + 1] = leftRotateFast(rotvec[0], j + 1);
	}
	NTL_EXEC_RANGE_END;

	BootContext bootContext = context.bootContextMap.at(logSlots);

	Ciphertext* tmpvec = new Ciphertext[k];

	NTL_EXEC_RANGE(k, first, last);
	for (j = first; j < last; ++j) {
		tmpvec[j] = multByPoly(rotvec[j], bootContext.pvecInv[j]);
	}
	NTL_EXEC_RANGE_END;

	for (j = 1; j < k; ++j) {
		addAndEqual(tmpvec[0], tmpvec[j]);
	}
	cipher = tmpvec[0];

	for (ki = k; ki < slots; ki+=k) {
		NTL_EXEC_RANGE(k, first, last);
		for (j = first; j < last; ++j) {
			tmpvec[j] = multByPoly(rotvec[j], bootContext.pvecInv[j + ki]);
		}
		NTL_EXEC_RANGE_END;

		for (j = 1; j < k; ++j) {
			addAndEqual(tmpvec[0], tmpvec[j]);
		}

		leftRotateAndEqualFast(tmpvec[0], ki);
		addAndEqual(cipher, tmpvec[0]);
	}
	delete[] rotvec;
	delete[] tmpvec;
}

void Scheme::evaluateExp2piAndEqual(Ciphertext& cipher, long logp) {
	Ciphertext cipher2 = square(cipher);
	reScaleByAndEqual(cipher2, logp); // cipher2.logq : logq - logp

	Ciphertext cipher4 = square(cipher2);
	reScaleByAndEqual(cipher4, logp); // cipher4.logq : logq -2logp

	RR c = 1/(2*Pi);
	ZZ pc = EvaluatorUtils::evalZZ(c, logp);

	Ciphertext cipher01 = addConst(cipher, pc); // cipher01.logq : logq

	c = 2*Pi;
	pc = EvaluatorUtils::evalZZ(c, logp);

	multByConstAndEqual(cipher01, pc);
	reScaleByAndEqual(cipher01, logp); // cipher01.logq : logq - logp

	c = 3/(2*Pi);
	pc = EvaluatorUtils::evalZZ(c, logp);
	Ciphertext cipher23 = addConst(cipher, pc); // cipher23.logq : logq

	c = 4*Pi*Pi*Pi/3;
	pc = EvaluatorUtils::evalZZ(c, logp);
	multByConstAndEqual(cipher23, pc);
	reScaleByAndEqual(cipher23, logp); // cipher23.logq : logq - logp

	multAndEqual(cipher23, cipher2);
	reScaleByAndEqual(cipher23, logp); // cipher23.logq : logq - 2logp

	addAndEqual(cipher23, cipher01); // cipher23.logq : logq - 2logp

	c = 5/(2*Pi);
	pc = EvaluatorUtils::evalZZ(c, logp);
	Ciphertext cipher45 = addConst(cipher, pc); // cipher45.logq : logq

	c = 4*Pi*Pi*Pi*Pi*Pi/15;
	pc = EvaluatorUtils::evalZZ(c, logp);
	multByConstAndEqual(cipher45, pc);
	reScaleByAndEqual(cipher45, logp); // cipher45.logq : logq - logp

	c = 7/(2*Pi);
	pc = EvaluatorUtils::evalZZ(c, logp);
	addConstAndEqual(cipher, pc); // cipher.logq : logq

	c = 8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
	pc = EvaluatorUtils::evalZZ(c, logp);
	multByConstAndEqual(cipher, pc);
	reScaleByAndEqual(cipher, logp); // cipher.logq : logq - logp

	multAndEqual(cipher, cipher2);
	reScaleByAndEqual(cipher, logp); // cipher.logq : logq - 2logp

	modDownByAndEqual(cipher45, logp); // cipher45.logq : logq - 2logp
	addAndEqual(cipher, cipher45); // cipher.logq : logq - 2logp

	multAndEqual(cipher, cipher4);
	reScaleByAndEqual(cipher, logp); // cipher.logq : logq - 3logp

	modDownByAndEqual(cipher23, logp);
	addAndEqual(cipher, cipher23); // cipher.logq : logq - 3logp
}

void Scheme::removeIpartAndEqual(Ciphertext& cipher, long logq, long logT, long logI) {
	long slots = cipher.slots;
	long logSlots = log2(slots);
	BootContext bootContext = context.bootContextMap.at(logSlots);
	if(logSlots == 0 && !cipher.isComplex) {
		imultAndEqual(cipher);
		reScaleByAndEqual(cipher, logT); // bitDown: logT
		evaluateExp2piAndEqual(cipher, logq + logI); // bitDown: logT + 3(logq + logI)
		for (long i = 0; i < logI + logT; ++i) {
			squareAndEqual(cipher);
			reScaleByAndEqual(cipher, logq + logI);
		}
		Ciphertext tmp = conjugate(cipher);
		subAndEqual(cipher, tmp);
		idivAndEqual(cipher);
		multByConstAndEqual(cipher, bootContext.c);
		// bitDown: logT + 3(logq + logI) + (logI + logT)(logq + logI)
	} else if(logSlots < context.logNh) {
		Ciphertext tmp = conjugate(cipher);
		subAndEqual(cipher, tmp);
		reScaleByAndEqual(cipher, logT + 1); // bitDown: logT + 1
		evaluateExp2piAndEqual(cipher, logq + logI); // bitDown: logT + 1 + 3(logq + logI)
		for (long i = 0; i < logI + logT; ++i) {
			squareAndEqual(cipher);
			reScaleByAndEqual(cipher, logq + logI);
		}
		tmp = conjugate(cipher);
		subAndEqual(cipher, tmp);

		tmp = multByPoly(cipher, bootContext.p1);
		Ciphertext tmprot = leftRotateFast(tmp, slots);
		addAndEqual(tmp, tmprot);
		multByPolyAndEqual(cipher, bootContext.p2);
		tmprot = leftRotateFast(cipher, slots);
		addAndEqual(cipher, tmprot);
		addAndEqual(cipher, tmp);
		// bitDown: logT + 1 + 3(logq + logI) + (logI + logT)(logq + logI)
	} else {
		Ciphertext tmp = conjugate(cipher);
		Ciphertext c2 = sub(cipher, tmp);
		addAndEqual(cipher, tmp);
		imultAndEqual(cipher);
		reScaleByAndEqual(cipher, logT + 1); // cipher bitDown: logT + 1
		reScaleByAndEqual(c2, logT + 1); // c2 bitDown: logT + 1
		evaluateExp2piAndEqual(cipher, logq + logI); // cipher bitDown: logT + 1 + 3(logq + logI)
		evaluateExp2piAndEqual(c2, logq + logI); // c2 bitDown: logT + 1 + 3(logq + logI)
		for (long i = 0; i < logI + logT; ++i) {
			squareAndEqual(c2);
			squareAndEqual(cipher);
			reScaleByAndEqual(c2, logq + logI);
			reScaleByAndEqual(cipher, logq + logI);
		}
		tmp = conjugate(c2);
		subAndEqual(c2, tmp);
		tmp = conjugate(cipher);
		subAndEqual(cipher, tmp);
		imultAndEqual(cipher);
		subAndEqual2(c2, cipher);
		multByConstAndEqual(cipher, bootContext.c);
		// bitDown: logT + 1 + 3(logq + logI) + (logI + logT)(logq + logI)
	}
	reScaleByAndEqual(cipher, logq + 2 * logI);
	// if (logSlots == 0 && !cipher.isComplex) bitDown: logT + 3(logq + logI) + (logI + logT)(logq + logI) + logq + 2logI
	// else bitDown: logT + 1 + 3(logq + logI) + (logI + logT)(logq + logI) + logq + 2logI
}

void Scheme::bootstrapAndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI) {
	long logSlots = log2(cipher.slots);

	modDownToAndEqual(cipher, logq);
	normalizeAndEqual(cipher);

	cipher.logq = logQ;
	cipher.q = power2_ZZ(logQ);

	for (long i = logSlots; i < context.logNh; ++i) {
		Ciphertext rot = leftRotateByPo2(cipher, i);
		addAndEqual(cipher, rot);
	}

	if (logSlots == 0 && !cipher.isComplex) {
			Ciphertext cconj = conjugate(cipher);
			addAndEqual(cipher, cconj);
			reScaleByAndEqual(cipher, context.logN - logSlots); // bitDown: context.logN - logSlots
			removeIpartAndEqual(cipher, logq, logT, logI); // bitDown: context.logN - logSlots + (logq + logI + 4) * logq + (logq + logI + 5) * logI + logT
	} else {
		reScaleByAndEqual(cipher, context.logNh - logSlots); // bitDown: context.logNh - logSlots
		linTransformAndEqual(cipher);
		reScaleByAndEqual(cipher, logq + logI + logSlots); // bitDown: context.logNh + logq + logI
		removeIpartAndEqual(cipher, logq, logT, logI); // bitDown: context.logNh + (logI + logT + 5) * logq + (logI + logT + 6) * logI + logT + 1
		linTransformInvAndEqual(cipher);
		reScaleByAndEqual(cipher, logq + logI); // bitDown: context.logNh + (logI + logT + 6) * logq + (logI + logT + 7) * logI + logT + 1
	}
}
