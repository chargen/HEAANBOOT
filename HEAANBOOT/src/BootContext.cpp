#include "BootContext.h"

BootContext::BootContext(ZZX* pvec, ZZX* pvecInv, long logp) : pvec(pvec), pvecInv(pvecInv), logp(logp) {}

//BootContext::~BootContext() {
//	delete[] pvec;
//	delete[] pvecInv;
//}
