#include <vector>

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include "../src/Common.h"
#include "../src/Ciphertext.h"
#include "../src/CZZ.h"
#include "../src/EvaluatorUtils.h"
#include "../src/NumUtils.h"
#include "../src/Params.h"
#include "../src/Scheme.h"
#include "../src/SchemeAlgo.h"
#include "../src/SecretKey.h"
#include "../src/StringUtils.h"
#include "../src/TimeUtils.h"
#include "../src/Context.h"

using namespace std;
using namespace NTL;

void TestMult(long logN, long logq, long precisionBits, long logSlots) {
	// Parameter Setting //
	TimeUtils timeutils;
	Params params(logN, logq);
	// SecretKey and PublicKey generate //
	timeutils.start("KeyGen");
	Context context(params);
	SecretKey secretKey(params);
	Scheme scheme(secretKey, context);
	SchemeAlgo algo(scheme);
	timeutils.stop("KeyGen");
	//-----------------------------------------
	long slots = (1 << logSlots);
	CZZ* mvec = EvaluatorUtils::evalRandCZZArray(slots, precisionBits);
	//-----------------------------------------
	timeutils.start("Encrypt batch");
	vector<Ciphertext> ciphers;
	for(long i = 0; i < 8; i++) {
		ciphers.push_back(scheme.encrypt(mvec, slots, logq));
	}
	timeutils.stop("Encrypt batch");
	//-----------------------------------------
	for(long i = 0; i < 4; i++) {
		timeutils.start("Homomorphic Multiplication");
		scheme.multAndEqual(ciphers[2*i], ciphers[2*i + 1]);
		scheme.reScaleByAndEqual(ciphers[2*i], precisionBits);
		timeutils.stop("Homomorphic Multiplication");
	}
	for(long i = 0; i < 2; i++) {
		timeutils.start("Homomorphic Multiplication");
		scheme.multAndEqual(ciphers[4*i], ciphers[4*i + 2]);
		scheme.reScaleByAndEqual(ciphers[4*i], precisionBits);
		timeutils.stop("Homomorphic Multiplication");
	}
	for(long i = 0; i < 1; i++) {
		timeutils.start("Homomorphic Multiplication");
		scheme.multAndEqual(ciphers[8*i], ciphers[8*i + 4]);
		scheme.reScaleByAndEqual(ciphers[8*i], precisionBits);
		timeutils.stop("Homomorphic Multiplication");
	}
	cout << "Ciphertext modulus bits = " << ciphers[0].cbits << endl;
	timeutils.start("Decrypt batch");
	CZZ* dvec = scheme.decrypt(secretKey, ciphers[0]);		

	timeutils.stop("Decrypt batch");
}

int main() {
	long depth = 3;
	long pBits = 16;
	long lBits = 10;
	long logq = pBits * (depth + 1) + lBits;
	long logN = 12;
	cout << "RLWE parameter (N, logq) = (" << (1 << logN) << "," << logq << ")" << endl; 
	TestMult(logN, logq, pBits, 10);

	return 0;
}

