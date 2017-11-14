#ifndef HEAAN_PUBLICKEY_H_
#define HEAAN_PUBLICKEY_H_

#include <map>

#include "BootKey.h"
#include "RLWE.h"
#include "SecretKey.h"

using namespace std;
using namespace NTL;

class PublicKey {
public:

	map<long, RLWE> keyMap;
	map<long, RLWE> leftRotKeyMap;
	map<long, BootKey> bootKeyMap;
	//-----------------------------------------

	PublicKey(Params& params, SecretKey& secretKey);

	//-----------------------------------------

	void addEncKey(Params& params, SecretKey& secretKey);
	void addConjKey(Params& params, SecretKey& secretKey);
	void addMultKey(Params& params, SecretKey& secretKey);

	void addLeftRotKey(Params& params, SecretKey& secretKey, long rot);
	void addLeftRotKeys(Params& params, SecretKey& secretKey);
	void addRightRotKeys(Params& params, SecretKey& secretKey);

	void addBootKeys(Params& params, SecretKey& secretKey, SchemeAux& aux, long logsize, long pBits);
	void addSortKeys(Params& params, SecretKey& secretKey, long size);
};

static long ENCRYPTION = 0;
static long MULTIPLICATION  = 1;
static long CONJUGATION = 2;

#endif
