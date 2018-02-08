#ifndef HEAAN_CONTEXT_H_
#define HEAAN_CONTEXT_H_

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <complex>

#include "Common.h"
#include "Params.h"
#include "BootContext.h"

using namespace std;
using namespace NTL;

static RR const Pi = ComputePi_RR();

static string LOGARITHM = "Logarithm"; ///< log(x)
static string EXPONENT  = "Exponent"; ///< exp(x)
static string SIGMOID   = "Sigmoid"; ///< sigmoid(x) = exp(x) / (1 + exp(x))

class Context {
public:

	long logN; ///< log of N
	long logQ; ///< log of Q
	double sigma; ///< standard deviation for Gaussian distribution
	long h; ///< parameter for HWT distribution

	long N; ///< N is a power-of-two that corresponds to the ring Z[X]/(X^N + 1)
	long Nh; ///< Nh = N/2
	long logNh; ///< logNh = logN - 1
	long M; ///< M = 2N
	long logPQ; ///< log of PQ

	ZZ Q; ///< Q corresponds to the highest modulus
	ZZ PQ; ///< PQ = Q * Q

	long* rotGroup; ///< auxiliary information about rotation group indexes for batch encoding
	RR* ksiPowsr; ///< storing ksi pows for fft calculation
	RR* ksiPowsi; ///< storing ksi pows for fft calculation

	ZZ* qpowvec;

	map<string, double*> taylorCoeffsMap; ///< storing taylor coefficients for function calculation
	map<long, BootContext> bootContextMap; ///< storing bootstrapping information, if generated

	Context(Params& params);

	/**
	 * adding information for Bootstrapping
	 * @param[in] logl: log of slots
	 * @param[in] logp: log of precision
	 */
	void addBootContext(long logSlots, long logp);

	/**
	 * encoding of values to polynomial, with bigger coefficients (used in Scheme.encode method)
	 * @param[in] vals: array of values
	 * @param[in] slots: size of array
	 */
	ZZX encode(complex<double>* vals, long slots, long logp);
	ZZX encode(double* vals, long slots, long logp);


	/**
	 * encoding of values to polynomial
	 * @param[in] vals: array of values
	 * @param[in] slots: size of array
	 */
	ZZX encodeSmall(complex<double>* vals, long slots, long logp);

	/**
	 * reverse bits for fft calculations
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 */
	void bitReverse(complex<double>* vals, const long size);

	/**
	 * calculates fft in Z_q[X] / (X^N + 1)
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 */
	void fft(complex<double>* vals, const long size);

	/**
	 * calculates fft inverse lazy in Z_q[X] / (X^N + 1)
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 */
	void fftInvLazy(complex<double>* vals, const long size);

	/**
	 * calculates fft inverse in Z_q[X] / (X^N + 1)
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 */
	void fftInv(complex<double>* vals, const long size);

	/**
	 * calculates special fft in Z_q[X] / (X^N + 1) for encoding/decoding
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 */
	void fftSpecial(complex<double>* vals, const long size);

	/**
	 * calculates special fft inverse lazy in Z_q[X] / (X^N + 1) for encoding/decoding
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 */
	void fftSpecialInvLazy(complex<double>* vals, const long size);

	/**
	 * calculates special fft inverse in Z_q[X] / (X^N + 1) for encoding/decoding
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 */
	void fftSpecialInv(complex<double>* vals, const long size);

	virtual ~Context();
};

#endif /* CONTEXT_H_ */
