#include "TestScheme.h"

#include <iostream>

using namespace std;

int main() {

	long logN;
	long logQ;
	long pBits;
	long logSlots;
	long testNum = 1;

	cout << "logN: ";
	cin >> logN;

	cout << "logQ: ";
	cin >> logQ;

	cout << "pBits: ";
	cin >> pBits;

	cout << "logSlots: ";
	cin >> logSlots;

	while(testNum != 0) {

		cout << "Lists of test" << endl;
		cout << "0. EXIT" << endl;
		cout << "1. testEncodeBatch" << endl;
		cout << "2. testBasic" << endl;

		cout << "\nTest Num: ";
		cin >> testNum;

		switch(testNum) {
			case 0:
				break;
			case 1:
				TestScheme::testEncodeBatch(logN, logQ, pBits, logSlots);
			case 2:
				TestScheme::testBasic(logN, logQ, pBits, logSlots);
		}

	}

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, logSlots
	 * Suggested: 13, 65, 30, 3
	 */
//	TestScheme::testEncodeBatch(13, 150, 30, 3);

	/*
	 * Params: logN, logQ, logp, isComplex
	 * Suggested: 13, 65, 30, 3
	 */
//	TestScheme::testEncodeSingle(13, 150, 30, false);

	/*
	 * Params: logN, logQ, logp, logSlots
	 * Suggested: 13, 65, 30, 3
	 */

//	TestScheme::testConjugateBatch(13, 65, 30, 3);

	/*
	 * Params: logN, logQ, logp, logSlots
	 * Suggested: 13, 65, 30, 3
	 */

//	TestScheme::testimultBatch(13, 65, 30, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, rotlogSlots, logSlots, isLeft
	 * Suggested: 13, 65, 30, 2, 5, true
	 */

//	TestScheme::testRotateByPo2Batch(13, 65, 30, 2, 5, true);
//	TestScheme::testRotateBatch(13, 65, 30, 17, 5, true);

	/*
	 * Params: logN, logQ, logp, logSlots
	 * Suggested: 13, 65, 30, 3
	 */

//	TestScheme::testSlotsSum(13, 65, 30, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, logDegree, logSlots
	 * Suggested: 13, 155, 30, 4, 3
	 * Suggested: 15, 618, 56, 10, 3
	 */

//	TestScheme::testPowerOf2Batch(13, 155, 30, 4, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, degree, logSlots
	 * Suggested: 13, 155, 30, 13, 3
	 * Suggested: 15, 618, 56, 903, 3
	 */

//	TestScheme::testPowerBatch(13, 155, 30, 13, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, logDegree, logSlots
	 * Suggested: 13, 155, 30, 4, 3
	 * Suggested: 15, 618, 56, 10, 3
	 */

//	TestScheme::testProdOfPo2Batch(13, 155, 30, 4, 3);

	/*
	 * Params: logN, logQ, logp, degree, logSlots
	 * Suggested: 13, 155, 30, 13, 3
	 * Suggested: 15, 618, 56, 903, 3
	 */

//	TestScheme::testProdBatch(13, 155, 30, 13, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, invSteps, logSlots
	 * Suggested: 14, 255, 25, 8, 3
	 * Suggested: 15, 325, 32, 8, 3
	 */

//	TestScheme::testInverseBatch(14, 255, 25, 8, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, degree, logSlots
	 * Suggested: 13, 155, 30, 7, 3
	 */

//	TestScheme::testLogarithmBatch(13, 155, 30, 7, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, degree, logSlots
	 * Suggested: 13, 155, 30, 7, 3
	 */

//	TestScheme::testExponentBatch(13, 155, 30, 7, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, degree, logSlots
	 * Suggested: 13, 155, 30, 7, 3
	 */

//	TestScheme::testSigmoidBatch(13, 155, 30, 7, 3);
//	TestScheme::testSigmoidBatchLazy(13, 155, 30, 7, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logQ, logp, logSlots, logFFTdim
	 * Suggested: 13, 100, 42, 3, 4;
	 */

//	TestScheme::testFFTBatch(13, 100, 42, 3, 4);
//	TestScheme::testFFTBatchLazy(13, 100, 42, 3, 4);

	/*
	 * Params: logN, logQ, logp, logSlots, logFFTdim, logHdim
	 * Suggested: 13, 140, 42, 3, 3, 2;
	 */

//	TestScheme::testFFTBatchLazyMultipleHadamard(13, 140, 42, 3, 3, 2);

	/*
	 * Params: logN, logQ, logp, logSlots
	 * Suggested: 13, 65, 30, 3
	 */

//	TestScheme::testCiphertextWriteAndRead(15, 620, 30, 3);

	//-----------------------------------------

	/*
	 * Params: logN, logq, logQ, logSlots, nu, logT
	 * Suggested: 15, 29, 620, 3, 6, 2
	 * Suggested: 15, 37, 620, 3, 10, 3
	 * Suggested: 16, 41, 1240, 3, 10, 3
	 * Suggested: 16, 54, 1240, 3, 15, 5
	 */
	TestScheme::testBootstrap(15, 29, 620, 3, 6, 2);

	/*
	 * Params: logN, logq, logQ, nu, logT
	 * Suggested: 15, 29, 620, 6, 2
	 * Suggested: 15, 37, 620, 10, 3
	 * Suggested: 16, 41, 1240, 10, 3
	 * Suggested: 16, 54, 1240, 15, 5
	 */
//	TestScheme::testBootstrapSingleReal(15, 29, 620, 6, 2);

//	TestScheme::test();

	return 0;
}
