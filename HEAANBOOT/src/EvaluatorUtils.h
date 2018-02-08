#ifndef HEAAN_EVALUATORUTILS_H_
#define HEAAN_EVALUATORUTILS_H_

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include "Context.h"

using namespace NTL;

class EvaluatorUtils {
public:


	//----------------------------------------------------------------------------------
	//   ZZ & CZZ TO DOUBLE
	//----------------------------------------------------------------------------------

	static double randomReal(double bound = 1.0);

	static complex<double> randomComplex(double bound = 1.0);

	static complex<double> randomCircle(double anglebound = 1.0);

	static double* randomRealArray(long size, double bound = 1.0);

	static complex<double>* randomComplexArray(long size, double bound = 1.0);

	static complex<double>* randomCircleArray(long size, double bound);

	/**
	 * evaluates double value (x >> logp)
	 * @param[in] x: ZZ scaled up value
	 * @param[in] logp: log of precision
	 * @return x >> logp
	 */
	static double evalReal(const ZZ& x, const long logp);

	/**
	 * evaluates double value (x >> logp)
	 * @param[out] res: double value (x >> logp)
	 * @param[in] x: ZZ scaled up value
	 * @param[in] logp: log of precision
	 */
	static void evalReal(double& res, const ZZ& x, const long logp);

	//----------------------------------------------------------------------------------
	//   DOUBLE & RR TO ZZ & CZZ
	//----------------------------------------------------------------------------------


	/**
	 * evaluates value x << logp
	 * @param[in] x: double value
	 * @param[in] logp: log of precision
	 * @return x << logp
	 */
	static ZZ evalZZ(const double x, const long logp);

	/**
	 * evaluates value x << logp
	 * @param[in] x: double value
	 * @param[in] logp: log of precision
	 * @return x << logp
	 */
	static ZZ evalZZ(const RR& x, const long logp);

	//----------------------------------------------------------------------------------
	//   ROTATIONS
	//----------------------------------------------------------------------------------


	/**
	 * left indexes rotation of values
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 * @param[in] rotSize: rotation size
	 */
	static void leftRotateAndEqual(complex<double>* vals, const long size, const long rotSize);

	/**
	 * right indexes rotation of values
	 * @param[in, out] vals: array of values
	 * @param[in] size: array size
	 * @param[in] rotSize: rotation size
	 */
	static void rightRotateAndEqual(complex<double>* vals, const long size, const long rotSize);

};

#endif
