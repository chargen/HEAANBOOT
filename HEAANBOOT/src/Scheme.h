#ifndef HEAAN_SCHEME_H_
#define HEAAN_SCHEME_H_

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>


#include "Common.h"
#include "CZZ.h"
#include "SecretKey.h"
#include "Ciphertext.h"
#include "Plaintext.h"
#include "Key.h"
#include "Context.h"
#include "BootKey.h"
#include "EvaluatorUtils.h"
#include "NumUtils.h"
#include "Params.h"
#include "Ring2Utils.h"

using namespace std;
using namespace NTL;

static long ENCRYPTION = 0;
static long MULTIPLICATION  = 1;
static long CONJUGATION = 2;

class Scheme {
private:
public:
	Context& context;
	map<long, Key> keyMap;
	map<long, Key> leftRotKeyMap;
	map<long, BootKey> bootKeyMap;

	Scheme(SecretKey& secretKey, Context& context);

	void addEncKey(SecretKey& secretKey);
	void addConjKey(SecretKey& secretKey);
	void addMultKey(SecretKey& secretKey);

	void addLeftRotKey(SecretKey& secretKey, long rot);
	void addLeftRotKeys(SecretKey& secretKey);
	void addRightRotKeys(SecretKey& secretKey);

	void addBootKeys(SecretKey& secretKey, long logsize, long pBits);
	void addSortKeys(SecretKey& secretKey, long size);

	//-----------------------------------------

	/**
	 * encodes vals into ZZX using fft inverse
	 * @param[in] vals
	 * @param[in] slots
	 * @return Message ZZX slots and level
	 */
	Plaintext encode(CZZ*& gvals, long slots, long cbits, bool isComplex = true);

	/**
	 * decodes ZZX into vals using fft
	 * @param[in] message
	 * @return array of CZZ values
	 */
	CZZ* decode(Plaintext& msg);

	Plaintext encodeSingle(CZZ& val, long cbits, bool isComplex = true);

	CZZ decodeSingle(Plaintext& msg);
	//-----------------------------------------

	/**
	 * encrypts message in RLWE instance bx = ax * sx + ex + mx mod 2^bits
	 * @param[in] message
	 * @param[in] bits
	 * @return cipher
	 */
	Ciphertext encryptMsg(Plaintext& msg);

	/**
	 * decrypts cipher
	 * @param[in] secret key
	 * @param[in] cipher
	 * @return message
	 */
	Plaintext decryptMsg(SecretKey& secretKey, Ciphertext& cipher);

	//-----------------------------------------

	/**
	 * All encryption process: regroup vals with adding conjugates, encode and encrypts
	 * @param[in] vals
	 * @param[in] slots
	 * @return cipher
	 */
	Ciphertext encrypt(CZZ*& vals, long slots, long cbits, bool isComplex = true);

	/**
	 * All decryption process: decrypt, decode, and degroup vals with removing conjugates
	 * @param[in] cipher
	 * @return vals
	 */
	CZZ* decrypt(SecretKey& secretKey, Ciphertext& cipher);

	//-----------------------------------------

	/**
	 * encrypt with only one value
	 * @param[in] val
	 * @return cipher
	 */
	Ciphertext encryptSingle(CZZ& val, long cbits, bool isComplex = true);

	/**
	 * decrypt with advance knowledge that slots = 1;
	 * @param[in] cipher
	 * @return val
	 */
	CZZ decryptSingle(SecretKey& secretKey, Ciphertext& cipher);

	//-----------------------------------------

	void normalizeAndEqual(Ciphertext& cipher);

	/**
	 * addition of ciphers
	 * @param[in] cipher(m1)
	 * @param[in] cipher(m2)
	 * @return cipher(m1 + m2)
	 */
	Ciphertext add(Ciphertext& cipher1, Ciphertext& cipher2);

	/**
	 * addition of ciphers
	 * @param[in, out] cipher(m1) -> cipher(m1 + m2)
	 * @param[in] cipher(m2)
	 */
	void addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	//-----------------------------------------

	/**
	 * constant addition
	 * @param[in] cipher(m)
	 * @param[in] constant
	 * @return cipher(m + constant)
	 */
	Ciphertext addConst(Ciphertext& cipher, ZZ& cnst);

	/**
	 * constant addition
	 * @param[in, out] cipher(m) -> cipher(m + constant)
	 * @param[in] constant
	 */
	void addConstAndEqual(Ciphertext& cipher, ZZ& cnst);

	//-----------------------------------------

	/**
	 * substraction of ciphers
	 * @param[in] cipher(m1)
	 * @param[in] cipher(m2)
	 * @return cipher(m1 - m2)
	 */
	Ciphertext sub(Ciphertext& cipher1, Ciphertext& cipher2);

	/**
	 * substraction of ciphers
	 * @param[in, out] cipher(m1) -> cipher(m1 - m2)
	 * @param[in] cipher(m2)
	 */
	void subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	/**
	 * substraction of ciphers
	 * @param[in] cipher(m1)
	 * @param[in, out] cipher(m2) -> cipher(m1 - m2)
	 */
	void subAndEqual2(Ciphertext& cipher1, Ciphertext& cipher2);

	/**
	 * conjugation in cipher
	 * @param[in] cipher(m = x + iy)
	 * @return cipher(x - iy)
	 */
	Ciphertext conjugate(Ciphertext& cipher);

	/**
	 * conjugation in cipher
	 * @param[in, out] cipher(m = x + iy) -> cipher(x - iy)
	 */
	void conjugateAndEqual(Ciphertext& cipher);

	/**
	 * multiplication by i (imaginary unit) in cipher
	 * @param[in] cipher(m)
	 * @return cipher(i * m)
	 */
	Ciphertext imult(Ciphertext& cipher, const long precisionBits);

	/**
	 * multiplication by i (imaginary unit) in cipher
	 * @param[in, out] cipher(m) -> cipher(i * m)
	 */
	void imultAndEqual(Ciphertext& cipher, const long precisionBits);

	//-----------------------------------------

	/**
	 * multiplication of ciphers. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after mult
	 * @param[in] cipher(m1)
	 * @param[in] cipher(m2)
	 * @return cipher(m1 * m2)
	 */
	Ciphertext mult(Ciphertext& cipher1, Ciphertext& cipher2);

	/**
	 * multiplication of ciphers. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after mult
	 * @param[in, out] cipher(m1) -> cipher(m1 * m2)
	 * @param[in] cipher(m2)
	 */
	void multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	//-----------------------------------------

	/**
	 * square of cipher. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after square
	 * @param[in] cipher(m)
	 * @return cipher(m^2)
	 */
	Ciphertext square(Ciphertext& cipher);

	/**
	 * square of cipher. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after square
	 * @param[in, out] cipher(m) -> cipher(m^2)
	 */
	void squareAndEqual(Ciphertext& cipher);

	//-----------------------------------------

	/**
	 * constant multiplication
	 * @param[in] cipher(m)
	 * @param[in] constant
	 * @return cipher(m * constant)
	 */
	Ciphertext multByConst(Ciphertext& cipher, ZZ& cnst);

	/**
	 * constant multiplication
	 * @param[in, out] cipher(m) -> cipher(m * constant)
	 * @param[in] constant
	 */
	void multByConstAndEqual(Ciphertext& cipher, ZZ& cnst);

	/**
	 * polynomial multiplication
	 * @param[in] cipher(m)
	 * @param[in] polynomial
	 * @return cipher(m * poly)
	 */
	Ciphertext multByPoly(Ciphertext& cipher, ZZX& poly);

	/**
	 * polynomial multiplication
	 * @param[in, out] cipher(m) -> cipher(m * poly)
	 * @param[in] polynomial
	 */
	void multByPolyAndEqual(Ciphertext& cipher, ZZX& poly);

	/**
	 * X^degree multiplication
	 * @param[in] cipher(m)
	 * @param[in] degree
	 * @return cipher(m * X^degree)
	 */
	Ciphertext multByMonomial(Ciphertext& cipher, const long degree);

	/**
	 * X^degree multiplication
	 * @param[in, out] cipher(m) -> cipher(m * X^degree)
	 * @param[in] degree
	 */
	void multByMonomialAndEqual(Ciphertext& cipher, const long degree);

	/**
	 * 2^bits multiplication
	 * @param[in] cipher(m)
	 * @param[in] bits
	 * @return cipher(m * 2^bits)
	 */
	Ciphertext leftShift(Ciphertext& cipher, long bits);

	/**
	 * 2^bits multiplication
	 * @param[in, out] cipher(m) -> cipher(m * 2^bits)
	 * @param[in] bits
	 */
	void leftShiftAndEqual(Ciphertext& cipher, long bits);

	/**
	 * doubles
	 * @param[in, out] cipher(m) -> cipher(2m)
	 */
	void doubleAndEqual(Ciphertext& cipher);

	/**
	 * rescaling procedure
	 * @param[in] cipher(m)
	 * @param[in] bitsDown
	 * @return cipher(m/2^bitsDown) with new cbits
	 */
	Ciphertext reScaleBy(Ciphertext& cipher, long bitsDown);

	Ciphertext reScaleTo(Ciphertext& cipher, long newcbits);
	/**
	 * rescaling procedure
	 * @param[in, out] cipher(m) -> cipher(m/2^bitsDown) with new cbits
	 * @param[in] bitsDown
	 */
	void reScaleByAndEqual(Ciphertext& cipher, long bitsDown);

	void reScaleToAndEqual(Ciphertext& cipher, long newcbits);
	/**
	 * modulus down procedure
	 * @param[in] cipher(m)
	 * @param[in] new cbits
	 * @return cipher(m) with new cbits
	 */
	Ciphertext modDownBy(Ciphertext& cipher, long bitsDown);

	Ciphertext modDownTo(Ciphertext& cipher, long cbits);
	/**
	 * modulus down procedure
	 * @param[in, out] cipher(m) -> cipher(m) with new cbits
	 * @param[in] new cbits
	 */
	void modDownByAndEqual(Ciphertext& cipher, long bitsDown);

	void modDownToAndEqual(Ciphertext& cipher, long cbits);
	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] log of rotation slots
	 * @return cipher(m(v_{1+2^logsteps}, v_{2+2^logsteps}, ..., v_{slots+2^logsteps})
	 */
	Ciphertext leftRotateByPo2(Ciphertext& cipher, long logrotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in, out] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1+2^logsteps}, v_{2+2^logsteps}, ..., v_{slots+2^logsteps})
	 * @param[in] log of rotation slots
	 * @return
	 */
	void leftRotateByPo2AndEqual(Ciphertext& cipher, long logrotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] log of rotation slots
	 * @return cipher(m(v_{1-2^logsteps}, v_{2-2^logsteps}, ..., v_{slots-2^logsteps})
	 */
	Ciphertext rightRotateByPo2(Ciphertext& cipher, long logrotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in, out] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1-2^logsteps}, v_{2-2^logsteps}, ..., v_{slots-2^logsteps})
	 * @param[in] log of rotation slots
	 * @return
	 */
	void rightRotateByPo2AndEqual(Ciphertext& cipher, long logrotSlots);

	Ciphertext leftRotateFast(Ciphertext& cipher, long rotSlots);

	void leftRotateAndEqualFast(Ciphertext& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] rotation slots
	 * @return cipher(m(v_{1+steps}, v_{2+steps}, ..., v_{slots+steps})
	 */
	Ciphertext leftRotate(Ciphertext& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1+steps}, v_{2+steps}, ..., v_{slots+steps})
	 * @param[in] rotation slots
	 */
	void leftRotateAndEqual(Ciphertext& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] rotation slots
	 * @return cipher(m(v_{1-steps}, v_{2-steps}, ..., v_{slots-steps})
	 */
	Ciphertext rightRotate(Ciphertext& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1-steps}, v_{2-steps}, ..., v_{slots-steps})
	 * @param[in] rotation slots
	 */
	void rightRotateAndEqual(Ciphertext& cipher, long rotSlots);

	//-----------------------------------------

	Ciphertext linearTransform(Ciphertext& cipher, long size);

	void linearTransformAndEqual(Ciphertext& cipher, long size);

	Ciphertext linearTransformInv(Ciphertext& cipher, long size);

	void linearTransformInvAndEqual(Ciphertext& cipher, long size);

	Ciphertext evaluateSin2pix7(Ciphertext& cipher, long pBits);

	void evaluateSin2pix7AndEqual(Ciphertext& cipher, long pBits);

	Ciphertext evaluateCos2pix6(Ciphertext& cipher, long pBits);

	void evaluateCos2pix6AndEqual(Ciphertext& cipher, long pBits);

	Ciphertext evaluateSin2x(Ciphertext& cSinx, Ciphertext& cCosx, long pBits);

	Ciphertext evaluateCos2x(Ciphertext& cSinx, Ciphertext& cCosx, long pBits);

	Ciphertext removeIpart(Ciphertext& cipher, long logq0, long logT, long logI = 4);

	void removeIpartAndEqual(Ciphertext& cipher, long logq0, long logT, long logI = 4);

	Ciphertext bootstrap(Ciphertext& cipher, long logq0, long logq, long logT, long logI = 4);

	void bootstrapAndEqual(Ciphertext& cipher, long logq0, long logq, long logT, long logI = 4);

	Ciphertext bootstrapOneReal(Ciphertext& cipher, long logq0, long logq, long logT, long logI = 4);

	void bootstrapOneRealAndEqual(Ciphertext& cipher, long logq0, long logq, long logT, long logI = 4);

};

#endif
