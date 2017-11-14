#ifndef HEAAN_BOOTKEY_H_
#define HEAAN_BOOTKEY_H_

#include "SchemeAux.h"

class BootKey {
public:
	long pBits;

	ZZX* pvec;
	ZZX* pvecInv;

	BootKey(Params& params, SchemeAux& aux, long pBits, long l);
};

#endif
