#include "Scheme.h"

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <cmath>

#include "EvaluatorUtils.h"
#include "NumUtils.h"
#include "Params.h"
#include "PubKey.h"
#include "Ring2Utils.h"

using namespace std;
using namespace NTL;

//-----------------------------------------

Message Scheme::encode(CZZ*& vals, long slots, long cbits, bool isComplex) {
	long doubleslots = slots << 1;
	ZZ mod = power2_ZZ(cbits);
	CZZ* gvals = new CZZ[doubleslots];
	for (long i = 0; i < slots; ++i) {
		long idx = (params.rotGroup[i] % (slots << 2) - 1) / 2;
		gvals[idx] = vals[i] << params.logq;
		gvals[doubleslots - idx - 1] = vals[i].conjugate() << params.logq;
	}

	ZZX mx;
	mx.SetLength(params.N);
	long idx = 0;
	long gap = params.N / doubleslots;

	NumUtils::fftSpecialInv(gvals, doubleslots, aux);

	for (long i = 0; i < doubleslots; ++i) {
		mx.rep[idx] = gvals[i].r;
		idx += gap;
	}
	delete[] gvals;
	return Message(mx, mod, cbits, slots);
}

CZZ* Scheme::decode(Message& msg) {
	long doubleslots = msg.slots * 2;
	CZZ* fftinv = new CZZ[doubleslots];

	long idx = 0;
	long gap = params.N / doubleslots;
	for (long i = 0; i < doubleslots; ++i) {
		ZZ tmp = msg.mx.rep[idx] % msg.mod;
		if(NumBits(tmp) == msg.cbits) tmp -= msg.mod;
		fftinv[i] = CZZ(tmp, ZZ(0));
		idx += gap;
	}
	NumUtils::fftSpecial(fftinv, doubleslots, aux);
	CZZ* res = new CZZ[msg.slots];
	for (long i = 0; i < msg.slots; ++i) {
		long idx = (params.rotGroup[i] % (msg.slots << 2) - 1) / 2;
		res[i] = fftinv[idx];
	}
	delete[] fftinv;
	return res;
}

Message Scheme::encodeSingle(CZZ& val, long cbits, bool isComplex) {
	ZZX mx;
	mx.SetLength(params.N);
	ZZ mod = power2_ZZ(cbits);
	if(isComplex) {
		CZZ gval = val << params.logq;
		mx.rep[0] = gval.r;
		mx.rep[params.N / 2] = gval.i;
		return Message(mx, mod, cbits);
	} else {
		mx.rep[0] = val.r << params.logq;
		return Message(mx, mod, cbits, 1, false);
	}
}

CZZ Scheme::decodeSingle(Message& msg) {
	CZZ res;
	res.r = msg.mx.rep[0];
	if(msg.isComplex) res.i = msg.mx.rep[params.N / 2];
	return res;
}

Cipher Scheme::encryptMsg(Message& msg) {
	ZZX ax, bx, vx, eax, ebx;
	NumUtils::sampleZO(vx, params.N);
	RLWE key = publicKey.keyMap.at(ENCRYPTION);

	ZZ Pmod = msg.mod << params.logq;

	Ring2Utils::mult(ax, vx, key.ax, Pmod, params.N);
	NumUtils::sampleGauss(eax, params.N, params.sigma);
	Ring2Utils::addAndEqual(ax, eax, Pmod, params.N);

	Ring2Utils::mult(bx, vx, key.bx, Pmod, params.N);
	NumUtils::sampleGauss(ebx, params.N, params.sigma);
	Ring2Utils::addAndEqual(bx, ebx, Pmod, params.N);

	Ring2Utils::addAndEqual(bx, msg.mx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(ax, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bx, params.logq, params.N);

	return Cipher(ax, bx, msg.mod, msg.cbits, msg.slots);
}

Message Scheme::decryptMsg(SecKey& secretKey, Cipher& cipher) {
	ZZX mx;
	Ring2Utils::mult(mx, cipher.ax, secretKey.sx, cipher.mod, params.N);
	Ring2Utils::addAndEqual(mx, cipher.bx, cipher.mod, params.N);
	return Message(mx, cipher.mod, cipher.cbits, cipher.slots, cipher.isComplex);
}

Cipher Scheme::encrypt(CZZ*& vals, long slots, long cbits, bool isComplex) {
	Message msg = encode(vals, slots, cbits, isComplex);
	return encryptMsg(msg);
}

CZZ* Scheme::decrypt(SecKey& secretKey, Cipher& cipher) {
	Message msg = decryptMsg(secretKey, cipher);
	return decode(msg);
}

Cipher Scheme::encryptSingle(CZZ& val, long cbits, bool isComplex) {
	Message msg = encodeSingle(val, cbits, isComplex);
	return encryptMsg(msg);
}

CZZ Scheme::decryptSingle(SecKey& secretKey, Cipher& cipher) {
	Message msg = decryptMsg(secretKey, cipher);
	return decodeSingle(msg);
}

//-----------------------------------------

void Scheme::normalizeAndEqual(Cipher& cipher) {
	for (int i = 0; i < params.N; ++i) {
		if(NumBits(cipher.ax.rep[i]) == cipher.cbits) cipher.ax.rep[i] -= cipher.mod;
		if(NumBits(cipher.bx.rep[i]) == cipher.cbits) cipher.bx.rep[i] -= cipher.mod;
	}
}

Cipher Scheme::add(Cipher& cipher1, Cipher& cipher2) {
	ZZX ax, bx;

	Ring2Utils::add(ax, cipher1.ax, cipher2.ax, cipher1.mod, params.N);
	Ring2Utils::add(bx, cipher1.bx, cipher2.bx, cipher1.mod, params.N);

	return Cipher(ax, bx, cipher1.mod, cipher1.cbits, cipher1.slots);
}

void Scheme::addAndEqual(Cipher& cipher1, Cipher& cipher2) {
	Ring2Utils::addAndEqual(cipher1.ax, cipher2.ax, cipher1.mod, params.N);
	Ring2Utils::addAndEqual(cipher1.bx, cipher2.bx, cipher1.mod, params.N);
}

//-----------------------------------------

Cipher Scheme::addConst(Cipher& cipher, ZZ& cnst) {
	ZZX ax = cipher.ax;
	ZZX bx = cipher.bx;

	AddMod(bx.rep[0], cipher.bx.rep[0], cnst, cipher.mod);
	return Cipher(ax, bx, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::addConstAndEqual(Cipher& cipher, ZZ& cnst) {
	ZZ mod = power2_ZZ(cipher.cbits);
	AddMod(cipher.bx.rep[0], cipher.bx.rep[0], cnst, mod);
}

//-----------------------------------------

Cipher Scheme::sub(Cipher& cipher1, Cipher& cipher2) {
	ZZX ax, bx;

	Ring2Utils::sub(ax, cipher1.ax, cipher2.ax, cipher1.mod, params.N);
	Ring2Utils::sub(bx, cipher1.bx, cipher2.bx, cipher1.mod, params.N);

	return Cipher(ax, bx, cipher1.mod, cipher1.cbits, cipher1.slots);
}

void Scheme::subAndEqual(Cipher& cipher1, Cipher& cipher2) {
	Ring2Utils::subAndEqual(cipher1.ax, cipher2.ax, cipher1.mod, params.N);
	Ring2Utils::subAndEqual(cipher1.bx, cipher2.bx, cipher1.mod, params.N);
}

void Scheme::subAndEqual2(Cipher& cipher1, Cipher& cipher2) {
	Ring2Utils::subAndEqual2(cipher1.ax, cipher2.ax, cipher1.mod, params.N);
	Ring2Utils::subAndEqual2(cipher1.bx, cipher2.bx, cipher1.mod, params.N);
}

Cipher Scheme::conjugate(Cipher& cipher) {
	ZZ Pmod = cipher.mod << params.logq;

	ZZX bxconj, bxres, axres;

	Ring2Utils::conjugate(bxconj, cipher.bx, params.N);
	Ring2Utils::conjugate(bxres, cipher.ax, params.N);

	RLWE key = publicKey.keyMap.at(CONJUGATION);

	Ring2Utils::mult(axres, bxres, key.ax, Pmod, params.N);
	Ring2Utils::multAndEqual(bxres, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(axres, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bxres, params.logq, params.N);

	Ring2Utils::addAndEqual(bxres, bxconj, cipher.mod, params.N);
	return Cipher(axres, bxres, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::conjugateAndEqual(Cipher& cipher) {
	ZZ Pmod = cipher.mod << params.logq;

	ZZX bxconj, bxres, axres;

	Ring2Utils::conjugate(bxconj, cipher.bx, params.N);
	Ring2Utils::conjugate(bxres, cipher.ax, params.N);

	RLWE key = publicKey.keyMap.at(CONJUGATION);
	Ring2Utils::mult(axres, bxres, key.ax, Pmod, params.N);
	Ring2Utils::multAndEqual(bxres, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(axres, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bxres, params.logq, params.N);

	Ring2Utils::addAndEqual(bxres, bxconj, cipher.mod, params.N);

	cipher.ax = axres;
	cipher.bx = bxres;
}

Cipher Scheme::imult(Cipher& cipher, const long precisionBits) {
	ZZX bxres, axres;
	Ring2Utils::multByMonomial(axres, cipher.ax, params.N / 2, params.N);
	Ring2Utils::multByMonomial(bxres, cipher.bx, params.N / 2, params.N);
	Cipher res(axres, bxres, cipher.mod, cipher.cbits, cipher.slots);
	return res;
}

void Scheme::imultAndEqual(Cipher& cipher, const long precisionBits) {
	Ring2Utils::multByMonomialAndEqual(cipher.ax, params.N / 2, params.N);
	Ring2Utils::multByMonomialAndEqual(cipher.bx, params.N / 2, params.N);
}

Cipher Scheme::mult(Cipher& cipher1, Cipher& cipher2) {
	ZZ Pmod = cipher1.mod << params.logq;

	ZZX axbx1 = Ring2Utils::add(cipher1.ax, cipher1.bx, cipher1.mod, params.N);
	ZZX axbx2 = Ring2Utils::add(cipher2.ax, cipher2.bx, cipher1.mod, params.N);
	Ring2Utils::multAndEqual(axbx1, axbx2, cipher1.mod, params.N);

	ZZX bxbx = Ring2Utils::mult(cipher1.bx, cipher2.bx, cipher1.mod, params.N);
	ZZX axax = Ring2Utils::mult(cipher1.ax, cipher2.ax, cipher1.mod, params.N);

	RLWE key = publicKey.keyMap.at(MULTIPLICATION);
	ZZX axmult = Ring2Utils::mult(axax, key.ax, Pmod, params.N);
	ZZX bxmult = Ring2Utils::mult(axax, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(axmult, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bxmult, params.logq, params.N);

	Ring2Utils::addAndEqual(axmult, axbx1, cipher1.mod, params.N);
	Ring2Utils::subAndEqual(axmult, bxbx, cipher1.mod, params.N);
	Ring2Utils::subAndEqual(axmult, axax, cipher1.mod, params.N);
	Ring2Utils::addAndEqual(bxmult, bxbx, cipher1.mod, params.N);

	return Cipher(axmult, bxmult, cipher1.mod, cipher1.cbits, cipher1.slots);
}

void Scheme::multAndEqual(Cipher& cipher1, Cipher& cipher2) {
	ZZ Pmod = cipher1.mod << params.logq;

	ZZX axbx1 = Ring2Utils::add(cipher1.ax, cipher1.bx, cipher1.mod, params.N);
	ZZX axbx2 = Ring2Utils::add(cipher2.ax, cipher2.bx, cipher1.mod, params.N);
	Ring2Utils::multAndEqual(axbx1, axbx2, cipher1.mod, params.N);

	ZZX bxbx = Ring2Utils::mult(cipher1.bx, cipher2.bx, cipher1.mod, params.N);
	ZZX axax = Ring2Utils::mult(cipher1.ax, cipher2.ax, cipher1.mod, params.N);

	RLWE key = publicKey.keyMap.at(MULTIPLICATION);
	cipher1.ax = Ring2Utils::mult(axax, key.ax, Pmod, params.N);
	cipher1.bx = Ring2Utils::mult(axax, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(cipher1.ax, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(cipher1.bx, params.logq, params.N);

	Ring2Utils::addAndEqual(cipher1.ax, axbx1, cipher1.mod, params.N);
	Ring2Utils::subAndEqual(cipher1.ax, bxbx, cipher1.mod, params.N);
	Ring2Utils::subAndEqual(cipher1.ax, axax, cipher1.mod, params.N);
	Ring2Utils::addAndEqual(cipher1.bx, bxbx, cipher1.mod, params.N);
}

//-----------------------------------------

Cipher Scheme::square(Cipher& cipher) {
	ZZ Pmod = cipher.mod << params.logq;

	ZZX axax, axbx, bxbx, bxmult, axmult;

	Ring2Utils::square(bxbx, cipher.bx, cipher.mod, params.N);
	Ring2Utils::mult(axbx, cipher.ax, cipher.bx, cipher.mod, params.N);
	Ring2Utils::addAndEqual(axbx, axbx, cipher.mod, params.N);
	Ring2Utils::square(axax, cipher.ax, cipher.mod, params.N);

	RLWE key = publicKey.keyMap.at(MULTIPLICATION);
	Ring2Utils::mult(axmult, axax, key.ax, Pmod, params.N);
	Ring2Utils::mult(bxmult, axax, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(axmult, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bxmult, params.logq, params.N);

	Ring2Utils::addAndEqual(axmult, axbx, cipher.mod, params.N);
	Ring2Utils::addAndEqual(bxmult, bxbx, cipher.mod, params.N);

	return Cipher(axmult, bxmult, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::squareAndEqual(Cipher& cipher) {
	ZZ Pmod = cipher.mod << params.logq;

	ZZX bxbx, axbx, axax, bxmult, axmult;

	Ring2Utils::square(bxbx, cipher.bx, cipher.mod, params.N);
	Ring2Utils::mult(axbx, cipher.bx, cipher.ax, cipher.mod, params.N);
	Ring2Utils::addAndEqual(axbx, axbx, cipher.mod, params.N);
	Ring2Utils::square(axax, cipher.ax, cipher.mod, params.N);

	RLWE key = publicKey.keyMap.at(MULTIPLICATION);
	Ring2Utils::mult(axmult, axax, key.ax, Pmod, params.N);
	Ring2Utils::mult(bxmult, axax, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(axmult, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bxmult, params.logq, params.N);

	Ring2Utils::addAndEqual(axmult, axbx, cipher.mod, params.N);
	Ring2Utils::addAndEqual(bxmult, bxbx, cipher.mod, params.N);

	cipher.bx = bxmult;
	cipher.ax = axmult;
}

//-----------------------------------------

Cipher Scheme::multByConst(Cipher& cipher, ZZ& cnst) {
	ZZX ax, bx;
	Ring2Utils::multByConst(ax, cipher.ax, cnst, cipher.mod, params.N);
	Ring2Utils::multByConst(bx, cipher.bx, cnst, cipher.mod, params.N);

	return Cipher(ax, bx, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::multByConstAndEqual(Cipher& cipher, ZZ& cnst) {
	Ring2Utils::multByConstAndEqual(cipher.ax, cnst, cipher.mod, params.N);
	Ring2Utils::multByConstAndEqual(cipher.bx, cnst, cipher.mod, params.N);
}

Cipher Scheme::multByPoly(Cipher& cipher, ZZX& poly) {
	ZZX axres, bxres;
	Ring2Utils::mult(axres, cipher.ax, poly, cipher.mod, params.N);
	Ring2Utils::mult(bxres, cipher.bx, poly, cipher.mod, params.N);
	return Cipher(axres, bxres, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::multByPolyAndEqual(Cipher& cipher, ZZX& poly) {
	Ring2Utils::multAndEqual(cipher.ax, poly, cipher.mod, params.N);
	Ring2Utils::multAndEqual(cipher.bx, poly, cipher.mod, params.N);
}

//-----------------------------------------

Cipher Scheme::multByMonomial(Cipher& cipher, const long degree) {
	ZZX ax, bx;

	Ring2Utils::multByMonomial(ax, cipher.ax, degree, params.N);
	Ring2Utils::multByMonomial(bx, cipher.bx, degree, params.N);

	return Cipher(ax, bx, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::multByMonomialAndEqual(Cipher& cipher, const long degree) {
	Ring2Utils::multByMonomialAndEqual(cipher.ax, degree, params.N);
	Ring2Utils::multByMonomialAndEqual(cipher.bx, degree, params.N);
}

//-----------------------------------------

Cipher Scheme::leftShift(Cipher& cipher, long bits) {
	ZZX ax, bx;

	Ring2Utils::leftShift(ax, cipher.ax, bits, cipher.mod, params.N);
	Ring2Utils::leftShift(bx, cipher.bx, bits, cipher.mod, params.N);

	return Cipher(ax, bx, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::leftShiftAndEqual(Cipher& cipher, long bits) {
	Ring2Utils::leftShiftAndEqual(cipher.ax, bits, cipher.mod, params.N);
	Ring2Utils::leftShiftAndEqual(cipher.bx, bits, cipher.mod, params.N);
}

void Scheme::doubleAndEqual(Cipher& cipher) {
	Ring2Utils::doubleAndEqual(cipher.ax, cipher.mod, params.N);
	Ring2Utils::doubleAndEqual(cipher.bx, cipher.mod, params.N);
}

//-----------------------------------------

Cipher Scheme::reScaleBy(Cipher& cipher, long bitsDown) {
	ZZX ax, bx;

	Ring2Utils::rightShift(ax, cipher.ax, bitsDown, params.N);
	Ring2Utils::rightShift(bx, cipher.bx, bitsDown, params.N);

	long newcbits = cipher.cbits - bitsDown;
	ZZ newmod = cipher.mod >> bitsDown;
	return Cipher(ax, bx, newmod, newcbits, cipher.slots);
}

Cipher Scheme::reScaleTo(Cipher& cipher, long newcbits) {
	ZZX ax, bx;

	long bitsDown = cipher.cbits - newcbits;
	Ring2Utils::rightShift(ax, cipher.ax, bitsDown, params.N);
	Ring2Utils::rightShift(bx, cipher.bx, bitsDown, params.N);

	ZZ newmod = power2_ZZ(newcbits);
	return Cipher(ax, bx, newmod, newcbits, cipher.slots);
}

void Scheme::reScaleByAndEqual(Cipher& cipher, long bitsDown) {
	Ring2Utils::rightShiftAndEqual(cipher.ax, bitsDown, params.N);
	Ring2Utils::rightShiftAndEqual(cipher.bx, bitsDown, params.N);
	cipher.cbits -= bitsDown;
	cipher.mod >>= bitsDown;
}

void Scheme::reScaleToAndEqual(Cipher& cipher, long newcbits) {
	long bitsDown = cipher.cbits - newcbits;
	Ring2Utils::rightShiftAndEqual(cipher.ax, bitsDown, params.N);
	Ring2Utils::rightShiftAndEqual(cipher.bx, bitsDown, params.N);
	cipher.cbits = newcbits;
	cipher.mod = power2_ZZ(newcbits);
}

Cipher Scheme::modDownBy(Cipher& cipher, long bitsDown) {
	ZZ newmod = cipher.mod >> bitsDown;
	long newcbits = cipher.cbits - bitsDown;
	ZZX bx, ax;
	Ring2Utils::mod(ax, cipher.ax, newmod, params.N);
	Ring2Utils::mod(bx, cipher.bx, newmod, params.N);
	return Cipher(ax, bx, newmod, newcbits, cipher.slots);
}

Cipher Scheme::modDownTo(Cipher& cipher, long newcbits) {
	ZZ newmod = power2_ZZ(newcbits);
	ZZX bx, ax;
	Ring2Utils::mod(ax, cipher.ax, newmod, params.N);
	Ring2Utils::mod(bx, cipher.bx, newmod, params.N);
	return Cipher(ax, bx, newmod, newcbits, cipher.slots);
}

void Scheme::modDownByAndEqual(Cipher& cipher, long bitsDown) {
	cipher.mod >>= bitsDown;
	cipher.cbits -= bitsDown;
	Ring2Utils::modAndEqual(cipher.ax, cipher.mod, params.N);
	Ring2Utils::modAndEqual(cipher.bx, cipher.mod, params.N);
}

void Scheme::modDownToAndEqual(Cipher& cipher, long newcbits) {
	if(cipher.cbits < newcbits) {
		invalid_argument("cbits of cipher should be larger than newcbits");
	}
	if(cipher.cbits == newcbits) {
		return;
	}
	cipher.mod = power2_ZZ(newcbits);
	cipher.cbits = newcbits;
	Ring2Utils::modAndEqual(cipher.ax, cipher.mod, params.N);
	Ring2Utils::modAndEqual(cipher.bx, cipher.mod, params.N);
}

Cipher Scheme::leftRotateFast(Cipher& cipher, long rotSlots) {
	ZZ Pmod = cipher.mod << params.logq;

	ZZX bxrot, bxres, axres;

	Ring2Utils::inpower(bxrot, cipher.bx, params.rotGroup[rotSlots], params.q, params.N);
	Ring2Utils::inpower(bxres, cipher.ax, params.rotGroup[rotSlots], params.q, params.N);

	RLWE key = publicKey.leftRotKeyMap.at(rotSlots);

	Ring2Utils::mult(axres, bxres, key.ax, Pmod, params.N);
	Ring2Utils::multAndEqual(bxres, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(axres, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bxres, params.logq, params.N);

	Ring2Utils::addAndEqual(bxres, bxrot, cipher.mod, params.N);
	return Cipher(axres, bxres, cipher.mod, cipher.cbits, cipher.slots);
}

void Scheme::leftRotateAndEqualFast(Cipher& cipher, long rotSlots) {
	ZZ Pmod = cipher.mod << params.logq;

	ZZX bxrot, bxres, axres;

	Ring2Utils::inpower(bxrot, cipher.bx, params.rotGroup[rotSlots], params.q, params.N);
	Ring2Utils::inpower(bxres, cipher.ax, params.rotGroup[rotSlots], params.q, params.N);

	RLWE key = publicKey.leftRotKeyMap.at(rotSlots);

	Ring2Utils::mult(axres, bxres, key.ax, Pmod, params.N);
	Ring2Utils::multAndEqual(bxres, key.bx, Pmod, params.N);

	Ring2Utils::rightShiftAndEqual(axres, params.logq, params.N);
	Ring2Utils::rightShiftAndEqual(bxres, params.logq, params.N);

	Ring2Utils::addAndEqual(bxres, bxrot, cipher.mod, params.N);

	cipher.ax = axres;
	cipher.bx = bxres;
}

Cipher Scheme::leftRotateByPo2(Cipher& cipher, long logrotSlots) {
	long rotSlots = (1 << logrotSlots);
	return leftRotateFast(cipher, rotSlots);
}

void Scheme::leftRotateByPo2AndEqual(Cipher& cipher, long logrotSlots) {
	long rotSlots = (1 << logrotSlots);
	leftRotateAndEqualFast(cipher, rotSlots);
}

Cipher Scheme::rightRotateByPo2(Cipher& cipher, long logrotSlots) {
	long rotSlots = params.N/2 - (1 << logrotSlots);
	return leftRotateFast(cipher, rotSlots);
}

void Scheme::rightRotateByPo2AndEqual(Cipher& cipher, long logrotSlots) {
	long rotSlots = params.N/2 - (1 << logrotSlots);
	leftRotateAndEqualFast(cipher, rotSlots);
}

Cipher Scheme::leftRotate(Cipher& cipher, long rotSlots) {
	Cipher res = cipher;
	leftRotateAndEqual(res, rotSlots);
	return res;
}

void Scheme::leftRotateAndEqual(Cipher& cipher, long rotSlots) {
	long remrotSlots = rotSlots % cipher.slots;
	long logrotSlots = log2((double)remrotSlots) + 1;
	for (long i = 0; i < logrotSlots; ++i) {
		if(bit(remrotSlots, i)) {
			leftRotateByPo2AndEqual(cipher, i);
		}
	}
}

Cipher Scheme::rightRotate(Cipher& cipher, long rotSlots) {
	Cipher res = cipher;
	rightRotateAndEqual(res, rotSlots);
	return res;
}

void Scheme::rightRotateAndEqual(Cipher& cipher, long rotSlots) {
	long remrotSlots = rotSlots % cipher.slots;
	long logrotSlots = log2((double)remrotSlots) + 1;
	for (long i = 0; i < logrotSlots; ++i) {
		if(bit(remrotSlots, i)) {
			rightRotateByPo2AndEqual(cipher, i);
		}
	}
}

Cipher Scheme::linearTransform(Cipher& cipher, long size) {
	long logSize = log2(size);
	long logSizeh = logSize / 2;
	long k = 1 << logSizeh;
	long m = 1 << (logSize - logSizeh);

	Cipher* encxrotvec = new Cipher[k];
	encxrotvec[0] = cipher;

	for (long i = 1; i < k; ++i) {
		encxrotvec[i] = leftRotateFast(encxrotvec[0], i);
	}

	BootKey bootKey = publicKey.bootKeyMap.at(logSize);
	Cipher res = multByPoly(encxrotvec[0], bootKey.pvec[0]);
	for (long j = 1; j < k; ++j) {
		Cipher cij = multByPoly(encxrotvec[j], bootKey.pvec[j]);
		addAndEqual(res, cij);
	}

	for (long i = 1; i < m; ++i) {
		long ki = k * i;
		Cipher ci0 = multByPoly(encxrotvec[0],bootKey.pvec[ki]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = multByPoly(encxrotvec[j], bootKey.pvec[j + ki]);
			addAndEqual(ci0, cij);
		}
		leftRotateAndEqualFast(ci0, ki);
		addAndEqual(res, ci0);
	}
	delete[] encxrotvec;
	return res;
}

Cipher Scheme::linearTransformInv(Cipher& cipher, long size) {
	long logSize = log2(size);
	long logSizeh = logSize / 2;
	long k = 1 << logSizeh;
	long m = 1 << (logSize - logSizeh);

	Cipher* cipherRotVec = new Cipher[k];
	cipherRotVec[0] = cipher;
	for (long i = 1; i < k; ++i) {
		cipherRotVec[i] = leftRotateFast(cipherRotVec[0], i);
	}
	BootKey bootKey = publicKey.bootKeyMap.at(logSize);

	Cipher res = multByPoly(cipherRotVec[0], bootKey.pvecInv[0]);

	for (long j = 1; j < k; ++j) {
		Cipher c0j = multByPoly(cipherRotVec[j], bootKey.pvecInv[j]);
		addAndEqual(res, c0j);
	}
	for (long i = 1; i < m; ++i) {
		long ki = k * i;
		Cipher ci0 = multByPoly(cipherRotVec[0], bootKey.pvecInv[ki]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = multByPoly(cipherRotVec[j], bootKey.pvecInv[j + ki]);
			addAndEqual(ci0, cij);
		}
		leftRotateAndEqualFast(ci0, ki);
		addAndEqual(res, ci0);
	}
	delete[] cipherRotVec;
	return res;
}

void Scheme::linearTransformAndEqual(Cipher& cipher, long size) {
	long logSize = log2(size);
	long logSizeh = logSize / 2;
	long k = 1 << logSizeh;
	long m = 1 << (logSize - logSizeh);

	Cipher* encxrotvec = new Cipher[k];
	encxrotvec[0] = cipher;

	for (long i = 1; i < k; ++i) {
		encxrotvec[i] = leftRotateFast(encxrotvec[0], i);
	}

	BootKey bootKey = publicKey.bootKeyMap.at(logSize);

	cipher = multByPoly(encxrotvec[0], bootKey.pvec[0]);
	for (long j = 1; j < k; ++j) {
		Cipher cij = multByPoly(encxrotvec[j], bootKey.pvec[j]);
		addAndEqual(cipher, cij);
	}

	for (long i = 1; i < m; ++i) {
		long ki = k * i;
		Cipher ci0 = multByPoly(encxrotvec[0], bootKey.pvec[ki]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = multByPoly(encxrotvec[j], bootKey.pvec[j + ki]);
			addAndEqual(ci0, cij);
		}
		leftRotateAndEqualFast(ci0, ki);
		addAndEqual(cipher, ci0);
	}
	delete[] encxrotvec;
}

void Scheme::linearTransformInvAndEqual(Cipher& cipher, long size) {
	long logSize = log2(size);
	long logSizeh = logSize / 2;
	long k = 1 << logSizeh;
	long m = 1 << (logSize - logSizeh);

	Cipher* cipherRotVec = new Cipher[k];
	cipherRotVec[0] = cipher;
	for (long i = 1; i < k; ++i) {
		cipherRotVec[i] = leftRotateFast(cipherRotVec[0], i);
	}
	BootKey bootKey = publicKey.bootKeyMap.at(logSize);

	cipher = multByPoly(cipherRotVec[0], bootKey.pvecInv[0]);

	for (long j = 1; j < k; ++j) {
		Cipher c0j = multByPoly(cipherRotVec[j], bootKey.pvecInv[j]);
		addAndEqual(cipher, c0j);
	}
	for (long i = 1; i < m; ++i) {
		long ki = k * i;
		Cipher ci0 = multByPoly(cipherRotVec[0], bootKey.pvecInv[ki]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = multByPoly(cipherRotVec[j], bootKey.pvecInv[j + ki]);
			addAndEqual(ci0, cij);
		}
		leftRotateAndEqualFast(ci0, ki);
		addAndEqual(cipher, ci0);
	}
	delete[] cipherRotVec;
}

Cipher Scheme::evaluateSin2pix7(Cipher& cipher, long pBits) {
	Cipher cipher2 = square(cipher); //depth 1
	reScaleByAndEqual(cipher2, pBits);
	Cipher cipher4 = square(cipher2); //depth 2
	reScaleByAndEqual(cipher4, pBits);
	RR c = -4*Pi*Pi*Pi/3;
	ZZ pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher tmp = multByConst(cipher, pc);
	reScaleByAndEqual(tmp, pBits); // depth 1

	c = -3/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher cipher13 = addConst(cipher2, pc);
	multAndEqual(cipher13, tmp);
	reScaleByAndEqual(cipher13, pBits); // depth 2

	c = -8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	tmp = multByConst(cipher, pc);
	reScaleByAndEqual(tmp, pBits); // depth 1

	c = -21/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher cipher57 = addConst(cipher2, pc);
	multAndEqual(cipher57, tmp);
	reScaleByAndEqual(cipher57, pBits); // depth 2
	multAndEqual(cipher57, cipher4);
	reScaleByAndEqual(cipher57, pBits); // depth 3

	modDownByAndEqual(cipher13, pBits); // depth 3
	addAndEqual(cipher57, cipher13); // depth 3
	return cipher57;
}

void Scheme::evaluateSin2pix7AndEqual(Cipher& cipher, long pBits) {
	Cipher cipher2 = square(cipher); //depth 1
	reScaleByAndEqual(cipher2, pBits);
	Cipher cipher4 = square(cipher2); //depth 2
	reScaleByAndEqual(cipher4, pBits);
	RR c = -4*Pi*Pi*Pi/3;
	ZZ pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher tmp = multByConst(cipher, pc);
	reScaleByAndEqual(tmp, pBits); // depth 1

	c = -3/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher cipher13 = addConst(cipher2, pc);
	multAndEqual(cipher13, tmp);
	reScaleByAndEqual(cipher13, pBits); // depth 2

	c = -8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	tmp = multByConst(cipher, pc);
	reScaleByAndEqual(tmp, pBits); // depth 1

	c = -21/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	cipher = addConst(cipher2, pc);
	multAndEqual(cipher, tmp);
	reScaleByAndEqual(cipher, pBits); // depth 2
	multAndEqual(cipher, cipher4);
	reScaleByAndEqual(cipher, pBits); // depth 3

	modDownByAndEqual(cipher13, pBits); // depth 3
	addAndEqual(cipher, cipher13); // depth 3
}

Cipher Scheme::evaluateCos2pix6(Cipher& cipher, long pBits) {
	Cipher cipher2 = square(cipher); //depth 1
	reScaleByAndEqual(cipher2, pBits);

	Cipher cipher4 = square(cipher2); //depth 2
	reScaleByAndEqual(cipher4, pBits);

	RR c = -1/(2*Pi*Pi);
	ZZ pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher cipher02 = addConst(cipher2, pc);

	c = -2*Pi*Pi;
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	multByConstAndEqual(cipher02, pc);
	reScaleByAndEqual(cipher02, pBits); // depth 2

	c = -15/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher cipher46 = addConst(cipher2, pc);

	c = -4*Pi*Pi*Pi*Pi*Pi*Pi/45;
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	multByConstAndEqual(cipher46, pc);
	reScaleByAndEqual(cipher46, pBits); // depth 2

	multAndEqual(cipher46, cipher4);
	reScaleByAndEqual(cipher46, pBits); // depth 3

	modDownByAndEqual(cipher02, pBits); // depth 3
	addAndEqual(cipher46, cipher02); // depth 3
	return cipher46;
}

void Scheme::evaluateCos2pix6AndEqual(Cipher& cipher, long pBits) {
	Cipher cipher2 = square(cipher); //depth 1
	reScaleByAndEqual(cipher2, pBits);

	Cipher cipher4 = square(cipher2); //depth 2
	reScaleByAndEqual(cipher4, pBits);

	RR c = -1/(2*Pi*Pi);
	ZZ pc = EvaluatorUtils::evaluateVal(c, pBits);
	Cipher cipher02 = addConst(cipher2, pc);

	c = -2*Pi*Pi;
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	multByConstAndEqual(cipher02, pc);
	reScaleByAndEqual(cipher02, pBits); // depth 2

	c = -15/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	cipher = addConst(cipher2, pc);

	c = -4*Pi*Pi*Pi*Pi*Pi*Pi/45;
	pc = EvaluatorUtils::evaluateVal(c, pBits);
	multByConstAndEqual(cipher, pc);
	reScaleByAndEqual(cipher, pBits); // depth 2

	multAndEqual(cipher, cipher4);
	reScaleByAndEqual(cipher, pBits); // depth 3

	modDownByAndEqual(cipher02, pBits); // depth 3
	addAndEqual(cipher, cipher02); // depth 3
}

Cipher Scheme::evaluateSin2x(Cipher& cSinx, Cipher& cCosx, long precisionBits) {
	Cipher res = mult(cSinx, cCosx);
	doubleAndEqual(res);
	reScaleByAndEqual(res, precisionBits);
	return res;
}

Cipher Scheme::evaluateCos2x(Cipher& cSinx, Cipher& cCosx, long precisionBits) {
	Cipher cSub = sub(cCosx, cSinx);
	Cipher cAdd = add(cCosx, cSinx);
	multAndEqual(cAdd, cSub);
	reScaleByAndEqual(cAdd, precisionBits);
	return cAdd;
}

Cipher Scheme::removeIpart(Cipher& cipher, long logq0, long logT, long logI) {
	Cipher cms = reScaleBy(cipher, logT);

	Cipher cipherSinx = evaluateSin2pix7(cms, logq0 + logI);
	Cipher cipherCosx = evaluateCos2pix6(cms, logq0 + logI);

	Cipher cipherSin2x, cipherCos2x;
	for (long i = 0; i < logI + logT - 1; ++i) {
		cipherSin2x = evaluateSin2x(cipherSinx, cipherCosx, logq0 + logI);
		cipherCos2x = evaluateCos2x(cipherSinx, cipherCosx, logq0 + logI);
		cipherSinx = cipherSin2x;
		cipherCosx = cipherCos2x;
	}
	cipherSinx = evaluateSin2x(cipherSinx, cipherCosx, logq0 + logI);
	ZZ temp = EvaluatorUtils::evaluateVal(1/(2*Pi), logq0 + logI);
	multByConstAndEqual(cipherSinx, temp);
	reScaleByAndEqual(cipherSinx, logq0 + 2 * logI);
	return cipherSinx;
}

void Scheme::removeIpartAndEqual(Cipher& cipher, long logq0, long logT, long logI) {
	Cipher cms = reScaleBy(cipher, logT);

	Cipher cipherSinx = evaluateSin2pix7(cms, logq0 + logI);
	Cipher cipherCosx = evaluateCos2pix6(cms, logq0 + logI);

	Cipher cipherSin2x, cipherCos2x;
	for (long i = 0; i < logI + logT - 1; ++i) {
		cipherSin2x = evaluateSin2x(cipherSinx, cipherCosx, logq0 + logI);
		cipherCos2x = evaluateCos2x(cipherSinx, cipherCosx, logq0 + logI);
		cipherSinx = cipherSin2x;
		cipherCosx = cipherCos2x;
	}
	cipher = evaluateSin2x(cipherSinx, cipherCosx, logq0 + logI);

	ZZ temp = EvaluatorUtils::evaluateVal(1/(2*Pi), logq0 + logI);
	multByConstAndEqual(cipher, temp);
	reScaleByAndEqual(cipher, logq0 + 2 * logI);
}

Cipher Scheme::bootstrap(Cipher& cipher, long logq0, long logq, long logT, long logI) {
	Cipher tmp = cipher;
	bootstrapAndEqual(tmp, logq0, logq, logT, logI);
	return tmp;
}

void Scheme::bootstrapAndEqual(Cipher& cipher, long logq0, long logq, long logT, long logI) {
	long logSlots = log2(cipher.slots);
	modDownToAndEqual(cipher, logq0);
	normalizeAndEqual(cipher);
	cipher.cbits = logq;
	cipher.mod = power2_ZZ(logq);

	if(logSlots == params.logN - 1) {
		Cipher cshift1 = multByMonomial(cipher, 2 * params.N - 1);
		linearTransformAndEqual(cipher, params.N / 2);
		linearTransformAndEqual(cshift1, params.N / 2);

		Cipher clinEvenConj = conjugate(cipher);
		addAndEqual(cipher, clinEvenConj);
		reScaleByAndEqual(cipher, logq0 + logI + logSlots);

		Cipher clinOddConj = conjugate(cshift1);
		addAndEqual(cshift1, clinOddConj);
		reScaleByAndEqual(cshift1, logq0 + logI + logSlots);

		removeIpartAndEqual(cipher, logq0, logT, logI);
		removeIpartAndEqual(cshift1, logq0, logT, logI);

		linearTransformInvAndEqual(cipher, params.N / 2);
		linearTransformInvAndEqual(cshift1, params.N / 2);

		multByMonomialAndEqual(cshift1, 1);
		addAndEqual(cipher, cshift1);
		reScaleByAndEqual(cipher, logq0 + logI);
	} else {
		for (long i = logSlots; i < params.logN - 1; ++i) {
			Cipher rot = leftRotateByPo2(cipher, i);
			addAndEqual(cipher, rot);
		}
		reScaleByAndEqual(cipher, params.logN - 1 - logSlots);
		linearTransformAndEqual(cipher, cipher.slots * 2);

		Cipher cipherConj = conjugate(cipher);
		addAndEqual(cipher, cipherConj);
		reScaleByAndEqual(cipher, logq0 + logI + logSlots + 2);

		removeIpartAndEqual(cipher, logq0, logT, logI);
		linearTransformInvAndEqual(cipher, cipher.slots * 2);
		reScaleByAndEqual(cipher, logq0 + logI);
	}
}
