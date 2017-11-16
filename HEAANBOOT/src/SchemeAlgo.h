#ifndef HEAAN_SCHEMEALGO_H_
#define HEAAN_SCHEMEALGO_H_

#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>

#include "Common.h"
#include "CZZ.h"
#include "EvaluatorUtils.h"
#include "Params.h"
#include "Plaintext.h"
#include "SecretKey.h"
#include "Ciphertext.h"
#include "Scheme.h"

class SchemeAlgo {
public:
	Scheme& scheme;

	//-----------------------------------------

	SchemeAlgo(Scheme& scheme) : scheme(scheme) {};

	//-----------------------------------------

	/**
	 * encrypting array of vals, each to one cipher
	 * @param[in] [m_1, m_2,...,m_size]
	 * @param[in] size
	 * @return [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 */
	Ciphertext* encryptSingleArray(CZZ*& vals, long size);

	/**
	 * decrypting array of ciphers with single val encrypted in each
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] size
	 * @return [m_1, m_2,...,m_size]
	 */
	CZZ* decryptSingleArray(SecretKey& secretKey, Ciphertext*& ciphers, long size);

	/**
	 * Calculating power of 2 cipher
	 * @param[in] cipher(m)
	 * @param[in] precision of initial m
	 * @param[in] logdeg
	 * @return cipher(m^2^logdeg)
	 */
	Ciphertext powerOf2(Ciphertext& cipher, const long precisionBits, const long logDegree);

	/**
	 * Calculating and storing all power of 2 ciphers up to 2^logdeg
	 * @param[in] cipher(m)
	 * @param[in] precision of initial m
	 * @param[in] logdeg
	 * @return [cipher(m), cipher(m^2), cipher(m^4), ... , cipher(m^2^logdeg)]
	 */
	Ciphertext* powerOf2Extended(Ciphertext& cipher, const long precisionBits, const long logDegree);

	//-----------------------------------------

	/**
	 * Calculating power of cipher
	 * @param[in] cipher(m)
	 * @param[in] precision of initial m
	 * @param[in] deg
	 * @return cipher(m^deg)
	 */
	Ciphertext power(Ciphertext& cipher, const long precisionBits, const long degree);

	/**
	 * Calculating and storing powers of cipher up to deg
	 * @param[in] cipher(m)
	 * @param[in] precision of initial m
	 * @param[in] deg
	 * @return [cipher(m), cipher(m^2), ... , cipher(m^deg)]
	 */
	Ciphertext* powerExtended(Ciphertext& cipher, const long precisionBits, const long degree);

	//-----------------------------------------

	/**
	 * Calculating product of ciphers, number of ciphers here is power of 2
	 * @param[in] [cipher(m_1), cipher(m_2), ... ,cipher(m_{2^logdeg})]
	 * @param[in] precision of initial m_i
	 * @param[in] logdeg
	 * @return cipher(m_1 * m_2 *...*m_{2^logdeg})
	 */
	Ciphertext prodOfPo2(Ciphertext*& ciphers, const long precisionBits, const long logDegree);

	/**
	 * Calculating product of ciphers
	 * @param[in] [cipher(m_1), cipher(m_2), ... ,cipher(m_{degree})]
	 * @param[in] precision of initial m_i
	 * @param[in] degree
	 * @return cipher(m_1 * m_2 *...*m_{degree})
	 */
	Ciphertext prod(Ciphertext*& ciphers, const long precisionBits, const long degree);

	//-----------------------------------------

	/**
	 * Calculating sum of ciphers
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] size
	 * @return cipher(m_1 + m_2 + ... + m_size)
	 */
	Ciphertext sum(Ciphertext*& ciphers, const long size);

	/**
	 * Calculating distance of ciphers
	 * @param[in] cipher(m_1, m_2, ..., m_slots)
	 * @param[in] cipher(m'_1, m'_2, ..., m'_slots)
	 * @return cipher(sum((m_i-m'_i)^2) >> precisionBits)
	 */
	Ciphertext distance(Ciphertext& cipher1, Ciphertext& cipher2, const long precisionBits);
	//-----------------------------------------

	/**
	 * Pairwise ciphers multiplication
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] [cipher(n_1), cipher(n_2),...,cipher(n_size)]
	 * @param[in] size
	 * @return [cipher(m_1 * n_1), cipher(m_2 * n_2),...,cipher(m_size * n_size)]
	 */
	Ciphertext* multVec(Ciphertext*& ciphers1, Ciphertext*& ciphers2, const long size);

	/**
	 * Pairwise ciphers multiplication
	 * @param[in, out] [cipher(m_1), cipher(m_2),...,cipher(m_size)] -> [cipher(m_1 * n_1), cipher(m_2 * n_2),...,cipher(m_size * n_size)]
	 * @param[in] [cipher(n_1), cipher(n_2),...,cipher(n_size)]
	 * @param[in] size
	 */
	void multAndEqualVec(Ciphertext*& ciphers1, Ciphertext*& ciphers2, const long size);

	/**
	 * Pairwise ciphers multiplication and modulus switching
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] [cipher(n_1), cipher(n_2),...,cipher(n_size)]
	 * @param[in] precision of initial m_i, n_i
	 * @param[in] size
	 * @return [cipher(m_1 * n_1 / p), cipher(m_2 * n_2 / p),...,cipher(m_size * n_size / p)] one level higher
	 */
	Ciphertext* multAndModSwitchVec(Ciphertext*& ciphers1, Ciphertext*& ciphers2, const long precisionBits, const long size);

	/**
	 * Pairwise ciphers multiplication and modulus switching
	 * @param[in, out] [cipher(m_1), cipher(m_2),...,cipher(m_size)] -> [cipher(m_1 * n_1 / p), cipher(m_2 * n_2 / p),...,cipher(m_size * n_size / p)]
	 * @param[in] [cipher(n_1), cipher(n_2),...,cipher(n_size)]
	 * @param[in] precision of initial m_i, n_i
	 * @param[in] size
	 */
	void multModSwitchAndEqualVec(Ciphertext*& ciphers1, Ciphertext*& ciphers2, const long precisionBits, const long size);

	//-----------------------------------------

	/**
	 * Calculating inner product of ciphers
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] [cipher(n_1), cipher(n_2),...,cipher(n_size)]
	 * @param[in] precision of initial m_i, n_i
	 * @param[in] size
	 * @return cipher(m_1 * n_1 + m_2 * n_2 + ... + m_size * n_size)
	 */
	Ciphertext innerProd(Ciphertext*& ciphers1, Ciphertext*& ciphers2, const long precision, const long size);

	/**
	 * Calculating cipher of partial sums
	 * @param[in] cipher(m_1, m_2,..., m_size)
	 * @param[in] slots summed in partial sums
	 * @return cipher(m_1 + ... + m_slots, m_2 + ... + m_{slots + 1},...,m_size + m_1 + ... + m_{slots - 1})
	 *
	 */
	Ciphertext partialSlotsSum(Ciphertext& cipher, const long slots);

	/**
	 * Calculating cipher of partial sums
	 * @param[in, out] cipher(m_1, m_2,..., m_size) -> cipher(m_1 + ... + m_slots, m_2 + ... + m_{slots + 1},...,m_size + m_1 + ... + m_{slots - 1})
	 * @param[in] slots summed in partial sums
	 */
	void partialSlotsSumAndEqual(Ciphertext& cipher, const long slots);

	//-----------------------------------------

	/**
	 * Calculating inverse of a cipher, using steps number of approximations
	 * @param[in] cipher(m)
	 * @param[in] precision of initial m
	 * @param[in] steps
	 * @return cipher(1 / m << (2 * precisionBits))
	 */
	Ciphertext inverse(Ciphertext& cipher, const long precisionBits, const long steps);

	/**
	 * Calculating and storing inverse of a cipher at each step up to stepsNum
	 * @param[in] cipher(m)
	 * @param[in] precision of initial m
	 * @param[in] steps
	 * @return [cipher(1 / m << (2 * precisionBits)), ... ,cipher(1 / m << (2 * precisionBits))]
	 */
	Ciphertext* inverseExtended(Ciphertext& cipher, const long precisionBits, const long steps);

	//-----------------------------------------

	/**
	 * Calculating function using Taylor Series approximation, more information in SchemeAux
	 * @param[in] cipher(m)
	 * @param[in] funcName
	 * @param[in] precision of initial m
	 * @param[in] degree
	 * @return cipher(funcName(m))
	 */
	Ciphertext function(Ciphertext& cipher, string& funcName, const long precisionBits, const long degree);

	/**
	 * Calculating function using Taylor Series approximation, more information in SchemeAux
	 * @param[in] cipher(m)
	 * @param[in] funcName
	 * @param[in] precision of initial m
	 * @param[in] degee
	 * @return cipher(funcName(m) * p), but saves one level
	 */
	Ciphertext functionLazy(Ciphertext& cipher, string& funcName, const long precisionBits, const long degree);

	/**
	 * Calculating function using Taylor Series approximation, and storing intermediate results
	 * @param[in] cipher(m)
	 * @param[in] funcName
	 * @param[in] precision of initial m
	 * @param[in] degee
	 * @return [cipher(funcName(m)), ... ,cipher(funcName(m))]
	 */
	Ciphertext* functionExtended(Ciphertext& cipher, string& funcName, const long precisionBits, const long degree);

	//-----------------------------------------

	/**
	 * Calculating fft of ciphers
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] precision of initial m_i
	 * @param[in] size is a power of 2
	 * @param[in] boolean is forward?
	 * @return [cipher(fft_1), ... ,cipher(fft_size)]
	 */
	void fftRaw(Ciphertext*& ciphers, const long size, const bool isForward);

	/**
	 * Calculating fft of ciphers
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] precision of initial m_i
	 * @param[in] size is a power of 2
	 * @return [cipher(fft_1), ... ,cipher(fft_size)]
	 */
	void fft(Ciphertext*& ciphers, const long size);

	/**
	 * Calculating fft inverse of ciphers
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] precision of initial m_i
	 * @param[in] size is a power of 2
	 * @return [cipher(fftinv_1), ... ,cipher(fftinv_size)]
	 */
	void fftInv(Ciphertext*& ciphers, const long size);

	/**
	 * Calculating fft inverse of ciphers
	 * @param[in] [cipher(m_1), cipher(m_2),...,cipher(m_size)]
	 * @param[in] precision of initial m_i
	 * @param[in] size is a power of 2
	 * @return [cipher(fftinv_1 * size), ... ,cipher(fftinv_size * size)] but saves level
	 */
	void fftInvLazy(Ciphertext*& ciphers, const long size);

	//-----------------------------------------

};

#endif
