#ifndef BOOTKEY_H_
#define BOOTKEY_H_

#include <SchemeAux.h>
#include <SecKey.h>

class BootKey {
public:
	long pBits;

	ZZX* pvec;
	ZZX* pvecInv;

	BootKey(Params& params, SchemeAux& aux, SecKey& secretKey, long pBits, long l);
};

#endif /* BOOTKEY_H_ */
