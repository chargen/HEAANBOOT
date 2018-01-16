#ifndef HEAAN_BOOTCONTEXT_H_
#define HEAAN_BOOTCONTEXT_H_

#include <NTL/ZZX.h>

using namespace NTL;

class BootContext {

public:

	ZZX* pvec; ///< encodings of "diagonal" values of linearization matrix
	ZZX* pvecInv; ///< encodings of "diagonal" values of inverse linearization matrix

	long logp; ///< number of preicision bits used in evaluation polynomials

	BootContext(ZZX* pvec = NULL, ZZX* pvecInv = NULL, long logp = 0);

};

#endif
