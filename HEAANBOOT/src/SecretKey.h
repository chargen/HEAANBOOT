#ifndef HEAAN_SECRETKEY_H_
#define HEAAN_SECRETKEY_H_

#include <NTL/ZZX.h>

#include <fstream>

#include "Params.h"
#include "NumUtils.h"


using namespace std;
using namespace NTL;

class SecretKey {
public:

	ZZX sx; ///< secret key

	SecretKey(Params& params);

	void Write(long SecretKeyID);

	void Read(long SecretKeyID);
};

#endif
