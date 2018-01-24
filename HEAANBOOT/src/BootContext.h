#ifndef HEAAN_BOOTCONTEXT_H_
#define HEAAN_BOOTCONTEXT_H_

#include <NTL/ZZX.h>

using namespace NTL;

class BootContext {

public:

	ZZX* pvec; ///< encodings of "diagonal" values of linearization matrix
	ZZX* pvecInv; ///< encodings of "diagonal" values of inverse linearization matrix

	ZZX p1; ///< auxiliary encoding for removeIpart
	ZZX p2; ///< auxiliary encoding for removeIpart

	ZZ c; ///< encoding of 1/4pi

	long logp; ///< number of preicision bits used in evaluation polynomials

	BootContext(ZZX* pvec = NULL, ZZX* pvecInv = NULL, ZZX p1 = ZZX::zero(), ZZX p2 = ZZX::zero(), ZZ c = ZZ::zero(), long logp = 0);

};

#endif
