#ifndef SCHEME_PUBKEY_H_
#define SCHEME_PUBKEY_H_

#include <map>

#include "BootKey.h"
#include "RLWE.h"
#include "SecKey.h"

using namespace std;
using namespace NTL;

class PubKey {
public:

	map<long, RLWE> keyMap;
	map<long, RLWE> leftRotKeyMap;
	map<long, BootKey> bootKeyMap;
	//-----------------------------------------

	PubKey(Params& params, SecKey& secretKey);

	//-----------------------------------------

	void addEncKey(Params& params, SecKey& secretKey);
	void addConjKey(Params& params, SecKey& secretKey);
	void addMultKey(Params& params, SecKey& secretKey);

	void addLeftRotKey(Params& params, SecKey& secretKey, long rot);
	void addLeftRotKeys(Params& params, SecKey& secretKey);
	void addRightRotKeys(Params& params, SecKey& secretKey);

	void addBootKeys(Params& params, SecKey& secretKey, SchemeAux& aux, long size, long pBits);
	void addSortKeys(Params& params, SecKey& secretKey, long size);
};

static long ENCRYPTION = 0;
static long MULTIPLICATION  = 1;
static long CONJUGATION = 2;

#endif /* SCHEME_PUBKEY_H_ */
