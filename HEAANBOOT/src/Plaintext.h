#ifndef HEAAN_PLAINTEXT_H_
#define HEAAN_PLAINTEXT_H_

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

using namespace std;
using namespace NTL;

class Plaintext {
public:

	ZZX mx; ///< message mod X^N + 1

	ZZ mod;
	long cbits; ///< bits in cipher
	long slots; ///< number of slots

	bool isComplex;
	//-----------------------------------------

	/**
	 * Plaintext: mx
	 * @param[in] polynomial mx
	 * @param[in] bits: bits in cipher
	 * @param[in] slots: number of slots
	 */
	Plaintext(ZZX mx = ZZX::zero(), ZZ mod = ZZ::zero(), long cbits = 0, long slots = 1, bool isComplex = true) : mx(mx), mod(mod), cbits(cbits), slots(slots), isComplex(isComplex) {}

	Plaintext(const Plaintext& o) : mx(o.mx), mod(o.mod), cbits(o.cbits), slots(o.slots), isComplex(o.isComplex) {}
};

#endif
