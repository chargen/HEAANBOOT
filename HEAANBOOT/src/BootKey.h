#ifndef HEAAN_BOOTKEY_H_
#define HEAAN_BOOTKEY_H_

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include "CZZ.h"
#include "EvaluatorUtils.h"
#include "NumUtils.h"
#include "Ring2Utils.h"
#include "Context.h"

class BootKey {
public:

	long logp; ///< number of preicision bits used in evaluation polynomials

	ZZX* pvec; ///< encodings of "diagonal" values of linearization matrix
	ZZX* pvecInv; ///< encodings of "diagonal" values of inverse linearization matrix

	BootKey(Context& context, long logp, long logl);
};

#endif
