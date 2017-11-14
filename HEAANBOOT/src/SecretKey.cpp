#include "SecretKey.h"

#include "NumUtils.h"

SecretKey::SecretKey(Params& params) {
	NumUtils::sampleHWT(sx, params.N, params.h);
}
