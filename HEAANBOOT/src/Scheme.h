#ifndef SCHEME_SCHEME_H_
#define SCHEME_SCHEME_H_

#include "CZZ.h"
#include "SecKey.h"
#include "Cipher.h"
#include "Message.h"
#include "PubKey.h"
#include "SchemeAux.h"

#include <stdexcept>

using namespace std;
using namespace NTL;

class Scheme {
private:
public:
	Params& params;
	PubKey& publicKey;
	SchemeAux& aux;


	Scheme(Params& params, PubKey& publicKey, SchemeAux& schemeaux) : params(params), publicKey(publicKey), aux(schemeaux) {};

	//-----------------------------------------

	/**
	 * encodes vals into ZZX using fft inverse
	 * @param[in] vals
	 * @param[in] slots
	 * @return Message ZZX slots and level
	 */
	Message encode(CZZ*& gvals, long slots, long cbits, bool isComplex = true);

	/**
	 * decodes ZZX into vals using fft
	 * @param[in] message
	 * @return array of CZZ values
	 */
	CZZ* decode(Message& msg);

	Message encodeSingle(CZZ& val, long cbits, bool isComplex = true);

	CZZ decodeSingle(Message& msg);
	//-----------------------------------------

	/**
	 * encrypts message in RLWE instance bx = ax * sx + ex + mx mod 2^bits
	 * @param[in] message
	 * @param[in] bits
	 * @return cipher
	 */
	Cipher encryptMsg(Message& msg);

	/**
	 * decrypts cipher
	 * @param[in] secret key
	 * @param[in] cipher
	 * @return message
	 */
	Message decryptMsg(SecKey& secretKey, Cipher& cipher);

	//-----------------------------------------

	/**
	 * All encryption process: regroup vals with adding conjugates, encode and encrypts
	 * @param[in] vals
	 * @param[in] slots
	 * @return cipher
	 */
	Cipher encrypt(CZZ*& vals, long slots, long cbits, bool isComplex = true);

	/**
	 * All decryption process: decrypt, decode, and degroup vals with removing conjugates
	 * @param[in] cipher
	 * @return vals
	 */
	CZZ* decrypt(SecKey& secretKey, Cipher& cipher);

	//-----------------------------------------

	/**
	 * encrypt with only one value
	 * @param[in] val
	 * @return cipher
	 */
	Cipher encryptSingle(CZZ& val, long cbits, bool isComplex = true);

	/**
	 * decrypt with advance knowledge that slots = 1;
	 * @param[in] cipher
	 * @return val
	 */
	CZZ decryptSingle(SecKey& secretKey, Cipher& cipher);

	//-----------------------------------------

	void normalizeAndEqual(Cipher& cipher);

	/**
	 * addition of ciphers
	 * @param[in] cipher(m1)
	 * @param[in] cipher(m2)
	 * @return cipher(m1 + m2)
	 */
	Cipher add(Cipher& cipher1, Cipher& cipher2);

	/**
	 * addition of ciphers
	 * @param[in, out] cipher(m1) -> cipher(m1 + m2)
	 * @param[in] cipher(m2)
	 */
	void addAndEqual(Cipher& cipher1, Cipher& cipher2);

	//-----------------------------------------

	/**
	 * constant addition
	 * @param[in] cipher(m)
	 * @param[in] constant
	 * @return cipher(m + constant)
	 */
	Cipher addConst(Cipher& cipher, ZZ& cnst);

	/**
	 * constant addition
	 * @param[in, out] cipher(m) -> cipher(m + constant)
	 * @param[in] constant
	 */
	void addConstAndEqual(Cipher& cipher, ZZ& cnst);

	//-----------------------------------------

	/**
	 * substraction of ciphers
	 * @param[in] cipher(m1)
	 * @param[in] cipher(m2)
	 * @return cipher(m1 - m2)
	 */
	Cipher sub(Cipher& cipher1, Cipher& cipher2);

	/**
	 * substraction of ciphers
	 * @param[in, out] cipher(m1) -> cipher(m1 - m2)
	 * @param[in] cipher(m2)
	 */
	void subAndEqual(Cipher& cipher1, Cipher& cipher2);

	/**
	 * substraction of ciphers
	 * @param[in] cipher(m1)
	 * @param[in, out] cipher(m2) -> cipher(m1 - m2)
	 */
	void subAndEqual2(Cipher& cipher1, Cipher& cipher2);

	/**
	 * conjugation in cipher
	 * @param[in] cipher(m = x + iy)
	 * @return cipher(x - iy)
	 */
	Cipher conjugate(Cipher& cipher);

	/**
	 * conjugation in cipher
	 * @param[in, out] cipher(m = x + iy) -> cipher(x - iy)
	 */
	void conjugateAndEqual(Cipher& cipher);

	/**
	 * multiplication by i (imaginary unit) in cipher
	 * @param[in] cipher(m)
	 * @return cipher(i * m)
	 */
	Cipher imult(Cipher& cipher, const long precisionBits);

	/**
	 * multiplication by i (imaginary unit) in cipher
	 * @param[in, out] cipher(m) -> cipher(i * m)
	 */
	void imultAndEqual(Cipher& cipher, const long precisionBits);

	//-----------------------------------------

	/**
	 * multiplication of ciphers. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after mult
	 * @param[in] cipher(m1)
	 * @param[in] cipher(m2)
	 * @return cipher(m1 * m2)
	 */
	Cipher mult(Cipher& cipher1, Cipher& cipher2);

	/**
	 * multiplication of ciphers. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after mult
	 * @param[in, out] cipher(m1) -> cipher(m1 * m2)
	 * @param[in] cipher(m2)
	 */
	void multAndEqual(Cipher& cipher1, Cipher& cipher2);

	//-----------------------------------------

	/**
	 * square of cipher. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after square
	 * @param[in] cipher(m)
	 * @return cipher(m^2)
	 */
	Cipher square(Cipher& cipher);

	/**
	 * square of cipher. This algorithm contain linearization.
	 * To manage the noise we usually need modular switching method after square
	 * @param[in, out] cipher(m) -> cipher(m^2)
	 */
	void squareAndEqual(Cipher& cipher);

	//-----------------------------------------

	/**
	 * constant multiplication
	 * @param[in] cipher(m)
	 * @param[in] constant
	 * @return cipher(m * constant)
	 */
	Cipher multByConst(Cipher& cipher, ZZ& cnst);

	/**
	 * constant multiplication
	 * @param[in, out] cipher(m) -> cipher(m * constant)
	 * @param[in] constant
	 */
	void multByConstAndEqual(Cipher& cipher, ZZ& cnst);

	/**
	 * polynomial multiplication
	 * @param[in] cipher(m)
	 * @param[in] polynomial
	 * @return cipher(m * poly)
	 */
	Cipher multByPoly(Cipher& cipher, ZZX& poly);

	/**
	 * polynomial multiplication
	 * @param[in, out] cipher(m) -> cipher(m * poly)
	 * @param[in] polynomial
	 */
	void multByPolyAndEqual(Cipher& cipher, ZZX& poly);

	/**
	 * X^degree multiplication
	 * @param[in] cipher(m)
	 * @param[in] degree
	 * @return cipher(m * X^degree)
	 */
	Cipher multByMonomial(Cipher& cipher, const long degree);

	/**
	 * X^degree multiplication
	 * @param[in, out] cipher(m) -> cipher(m * X^degree)
	 * @param[in] degree
	 */
	void multByMonomialAndEqual(Cipher& cipher, const long degree);

	/**
	 * 2^bits multiplication
	 * @param[in] cipher(m)
	 * @param[in] bits
	 * @return cipher(m * 2^bits)
	 */
	Cipher leftShift(Cipher& cipher, long bits);

	/**
	 * 2^bits multiplication
	 * @param[in, out] cipher(m) -> cipher(m * 2^bits)
	 * @param[in] bits
	 */
	void leftShiftAndEqual(Cipher& cipher, long bits);

	/**
	 * doubles
	 * @param[in, out] cipher(m) -> cipher(2m)
	 */
	void doubleAndEqual(Cipher& cipher);

	/**
	 * rescaling procedure
	 * @param[in] cipher(m)
	 * @param[in] bitsDown
	 * @return cipher(m/2^bitsDown) with new cbits
	 */
	Cipher reScaleBy(Cipher& cipher, long bitsDown);

	Cipher reScaleTo(Cipher& cipher, long newcbits);
	/**
	 * rescaling procedure
	 * @param[in, out] cipher(m) -> cipher(m/2^bitsDown) with new cbits
	 * @param[in] bitsDown
	 */
	void reScaleByAndEqual(Cipher& cipher, long bitsDown);

	void reScaleToAndEqual(Cipher& cipher, long newcbits);
	/**
	 * modulus down procedure
	 * @param[in] cipher(m)
	 * @param[in] new cbits
	 * @return cipher(m) with new cbits
	 */
	Cipher modDownBy(Cipher& cipher, long bitsDown);

	Cipher modDownTo(Cipher& cipher, long cbits);
	/**
	 * modulus down procedure
	 * @param[in, out] cipher(m) -> cipher(m) with new cbits
	 * @param[in] new cbits
	 */
	void modDownByAndEqual(Cipher& cipher, long bitsDown);

	void modDownToAndEqual(Cipher& cipher, long cbits);
	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] log of rotation slots
	 * @return cipher(m(v_{1+2^logsteps}, v_{2+2^logsteps}, ..., v_{slots+2^logsteps})
	 */
	Cipher leftRotateByPo2(Cipher& cipher, long logrotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in, out] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1+2^logsteps}, v_{2+2^logsteps}, ..., v_{slots+2^logsteps})
	 * @param[in] log of rotation slots
	 * @return
	 */
	void leftRotateByPo2AndEqual(Cipher& cipher, long logrotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] log of rotation slots
	 * @return cipher(m(v_{1-2^logsteps}, v_{2-2^logsteps}, ..., v_{slots-2^logsteps})
	 */
	Cipher rightRotateByPo2(Cipher& cipher, long logrotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in, out] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1-2^logsteps}, v_{2-2^logsteps}, ..., v_{slots-2^logsteps})
	 * @param[in] log of rotation slots
	 * @return
	 */
	void rightRotateByPo2AndEqual(Cipher& cipher, long logrotSlots);

	Cipher leftRotateFast(Cipher& cipher, long rotSlots);

	void leftRotateAndEqualFast(Cipher& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] rotation slots
	 * @return cipher(m(v_{1+steps}, v_{2+steps}, ..., v_{slots+steps})
	 */
	Cipher leftRotate(Cipher& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1+steps}, v_{2+steps}, ..., v_{slots+steps})
	 * @param[in] rotation slots
	 */
	void leftRotateAndEqual(Cipher& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots))
	 * @param[in] rotation slots
	 * @return cipher(m(v_{1-steps}, v_{2-steps}, ..., v_{slots-steps})
	 */
	Cipher rightRotate(Cipher& cipher, long rotSlots);

	/**
	 * calculates cipher of array with rotated indexes
	 * @param[in] cipher(m(v_1, v_2, ..., v_slots)) -> cipher(m(v_{1-steps}, v_{2-steps}, ..., v_{slots-steps})
	 * @param[in] rotation slots
	 */
	void rightRotateAndEqual(Cipher& cipher, long rotSlots);

	//-----------------------------------------

	Cipher linearTransform(Cipher& cipher, long size);

	void linearTransformAndEqual(Cipher& cipher, long size);

	Cipher linearTransformInv(Cipher& cipher, long size);

	void linearTransformInvAndEqual(Cipher& cipher, long size);

	Cipher evaluateSin2pix7(Cipher& cipher, long pBits);

	void evaluateSin2pix7AndEqual(Cipher& cipher, long pBits);

	Cipher evaluateCos2pix6(Cipher& cipher, long pBits);

	void evaluateCos2pix6AndEqual(Cipher& cipher, long pBits);

	Cipher evaluateSin2x(Cipher& cSinx, Cipher& cCosx, long pBits);

	Cipher evaluateCos2x(Cipher& cSinx, Cipher& cCosx, long pBits);

	Cipher removeIpart(Cipher& cipher, long logq0, long logT, long logI = 4);

	void removeIpartAndEqual(Cipher& cipher, long logq0, long logT, long logI = 4);

	Cipher bootstrap(Cipher& cipher, long logq0, long logq, long logT, long logI = 4);

	void bootstrapAndEqual(Cipher& cipher, long logq0, long logq, long logT, long logI = 4);

	Cipher bootstrapOneReal(Cipher& cipher, long logq0, long logq, long logT, long logI = 4);

	void bootstrapOneRealAndEqual(Cipher& cipher, long logq0, long logq, long logT, long logI = 4);

};

#endif /* SCHEME_SCHEME_H_ */
