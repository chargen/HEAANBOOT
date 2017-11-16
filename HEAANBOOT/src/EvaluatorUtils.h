#ifndef HEAAN_EVALUATORUTILS_H_
#define HEAAN_EVALUATORUTILS_H_

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include "Context.h"
#include "CZZ.h"

using namespace NTL;

class EvaluatorUtils {
public:

	/**
	 * evaluates value xr << bits
	 * @param[in] x
	 * @param[in] bits
	 * @return x << bits
	 */
	static ZZ evalZZ(const double& x, const long& bits);

	/**
	 * evaluates value xr << bits
	 * @param[in] x
	 * @param[in] bits
	 * @return x << bits
	 */
	static ZZ evalZZ(const RR& x, const long& bits);

	/**
	 * evaluates Z[i] value (xr + i * xi) << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value (xr + i * xi) << bits
	 */
	static CZZ evalCZZ(const double& xr, const double& xi, const long& bits);

	/**
	 * evaluates Z[i] value (xr + i * xi) << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value (xr + i * xi) << bits
	 */
	static CZZ evalCZZ(const RR& xr, const RR& xi, const long& bits);

	/**
	 * evaluates random bits in Z[i]
	 * @param[in] bits
	 * @return random bits in Z[i]
	 */
	static CZZ evalRandCZZ(const long& bits);

	static CZZ evalRandCZZ0(const long& bits);

	/**
	 * evaluates random bits in Z[i]
	 * @param[in] bits
	 * @return random bits in Z[i]
	 */
	static CZZ evalRandCZZCircle(const long& bits);

	static ZZ* evalRandZZArray(const long& size, const long& bits);

	/**
	 * evaluates array of random bits in Z[i]
	 * @param[in] size of array
	 * @param[in] bits
	 * @return array of random bits in Z[i]
	 */
	static CZZ* evalRandCZZArray(const long& size, const long& bits);

	static CZZ* evalRandCZZ0Array(const long& size, const long& bits);

	//-----------------------------------------

	/**
	 * evaluates Z[i] value (xr + i * xi)^d << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] power degree
	 * @param[in] bits
	 * @return Z[i] value (xr + i * xi)^d << bits
	 */
	static CZZ evalCZZPow(const double& xr, const double& xi, const long& degree, const long& bits);

	/**
	 * evaluates Z[i] value (xr + i * xi)^d << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] power degree
	 * @param[in] bits
	 * @return Z[i] value (xr + i * xi)^d << bits
	 */
	static CZZ evalCZZPow(const RR& xr, const RR& xi, const long& degree, const long& bits);

	/**
	 * evaluates Z[i] value (xr + i * xi)^d << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] log(degree)
	 * @param[in] bits
	 * @return Z[i] value (xr + i * xi)^d << bits
	 */
	static CZZ evalCZZPow2(const double& xr, const double& xi, const long& logDegree, const long& bits);

	/**
	 * evaluates Z[i] value (xr + i * xi)^d << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] log(degree)
	 * @param[in] bits
	 * @return Z[i] value (xr + i * xi)^d << bits
	 */
	static CZZ evalCZZPow2(const RR& xr, const RR& xi, const long& logDegree, const long& bits);

	//-----------------------------------------

	/**
	 * evaluates array of Z[i] values (xr + i * xi)^j << bits for j=1,...,d
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] power degree
	 * @param[in] bits
	 * @return array of Z[i] values (xr + i * xi)^j << bits for j=1,...,d
	 */
	static CZZ* evalCZZPowArray(const double& xr, const double& xi, const long& degree, const long& bits);

	/**
	 * evaluates array of Z[i] values (xr + i * xi)^j << bits for j=1,...,d
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] power degree
	 * @param[in] bits
	 * @return array of Z[i] values (xr + i * xi)^j << bits for j=1,...,d
	 */
	static CZZ* evalCZZPowArray(const RR& xr, const RR& xi, const long& degree, const long& bits);

	/**
	 * evaluates array of Z[i] values (xr + i * xi)^j << bits for j=1,2,2^2,...,d
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] log(d)
	 * @param[in] bits
	 * @return array of Z[i] values (xr + i * xi)^j << bits for j=1,2,2^2,...,d
	 */
	static CZZ* evalCZZPow2Array(const double& xr, const double& xi, const long& logDegree, const long& bits);

	/**
	 * evaluates array of Z[i] values (xr + i * xi)^j << bits for j=1,2,2^2,...,d
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] log(d)
	 * @param[in] bits
	 * @return array of Z[i] values (xr + i * xi)^j << bits for j=1,2,2^2,...,d
	 */
	static CZZ* evalCZZPow2Array(const RR& xr, const RR& xi, const long& logDegree, const long& bits);

	//-----------------------------------------

	/**
	 * evaluates Z[i] value 1 / (xr + i * xi) << (2 * bits)
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value 1 / (xr + i * xi) << (2 * bits)
	 */
	static CZZ evalCZZInv(const double& xr, const double& xi, const long& bits);

	/**
	 * evaluates Z[i] value 1 / (xr + i * xi) << (2 * bits)
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] log(p)
	 * @return Z[i] value 1 / (xr + i * xi) << (2 * bits)
	 */
	static CZZ evalCZZInv(const RR& xr, const RR& xi, const long& bits);

	/**
	 * evaluates Z[i] value log(xr + i * xi) << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value log(xr + i * xi) << bits
	 */
	static CZZ evalCZZLog(const double& xr, const double& xi, const long& bits);

	/**
	 * evaluates Z[i] value exp(xr + i * xi) << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value exp(xr + i * xi) << bits
	 */
	static CZZ evalCZZExp(const double& xr, const double& xi, const long& bits);

	/**
	 * evaluates Z[i] value exp(xr + i * xi) << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value exp(xr + i * xi) << bits
	 */
	static CZZ evalCZZExp(const RR& xr, const RR& xi, const long& bits);

	/**
	 * evaluates Z[i] value sigmoid(xr + i * xi) << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value sigmoid(xr + i * xi) << bits
	 */
	static CZZ evalCZZSigmoid(const double& xr, const double& xi, const long& bits);

	/**
	 * evaluates Z[i] value sigmoid(xr + i * xi) << bits
	 * @param[in] real part
	 * @param[in] imaginary part
	 * @param[in] bits
	 * @return Z[i] value sigmoid(xr + i * xi) << bits
	 */
	static CZZ evalCZZSigmoid(const RR& xr, const RR& xi, const long& bits);

	//-----------------------------------------

	/**
	 * left shift array of values by bits
	 * @param[in, out] array of values
	 * @param[in] array size
	 * @param[in] bits
	 */
	static void leftShiftAndEqual(CZZ*& vals, const long& size, const long& bits);

	//-----------------------------------------

	/**
	 * left indexes rotation of values
	 * @param[in, out] array of values
	 * @param[in] array size
	 * @param[in] rotation size
	 */
	static void leftRotateAndEqual(CZZ*& vals, const long& size, const long& rotSize);

	/**
	 * right indexes rotation of values
	 * @param[in, out] array of values
	 * @param[in] array size
	 * @param[in] rotation size
	 */
	static void rightRotateAndEqual(CZZ*& vals, const long& size, const long& rotSize);
	//-----------------------------------------
};

#endif
