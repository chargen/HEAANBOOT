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
	long pBits;

	ZZX* pvec;
	ZZX* pvecInv;

	BootKey(Context& context, long pBits, long l);
};

#endif
