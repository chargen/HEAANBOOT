#ifndef HEAAN_BOOTCONTEXT_H_
#define HEAAN_BOOTCONTEXT_H_

#include "Context.h"

class BootContext {
public:

	long logp; ///< number of preicision bits used in evaluation polynomials

	ZZX* pvec; ///< encodings of "diagonal" values of linearization matrix
	ZZX* pvecInv; ///< encodings of "diagonal" values of inverse linearization matrix

	BootContext(Context& context, long logp, long logl);
};

#endif
