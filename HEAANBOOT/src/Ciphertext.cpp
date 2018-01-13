#include "Ciphertext.h"

#include <NTL/tools.h>
#include "Common.h"

void Ciphertext::Write(long CiphertextID) {
	ofstream myfile;
	myfile.open("Ciphertext" + to_string(CiphertextID) + ".txt");
	myfile << "Ciphertext Information" << endl;
	myfile << "CiphertextID = " << CiphertextID << endl;
	myfile << deg(ax) << endl;
	myfile << deg(bx) << endl;
	myfile << logq << endl;
	myfile << slots << endl;
	myfile << isComplex << endl;
	for(long i = 0; i < deg(ax) + 1; i++) {
		myfile << ax[i] << endl;
	}
	for(long i = 0; i < deg(bx) + 1; i++) {
		myfile << bx[i] << endl;
	}
	myfile.close();
}

void Ciphertext::Read(long CiphertextID) {
	ifstream myfile("Ciphertext" + to_string(CiphertextID) + ".txt");
	if(myfile.is_open()) {
		// kill previous memory
		ax.kill();
		bx.kill();
		// start reading
		long temp;
		string line;
		// pass first two lines
		getline(myfile, line);
		getline(myfile, line);

		// read 3rd line and get degree of ax
		getline(myfile, line);
		temp = atol(line.c_str());
		ax.SetLength(temp + 1);

		// read 4th line and get degree of bx
		getline(myfile, line);
		temp = atol(line.c_str());
		bx.SetLength(temp + 1);

		// read 5th line and get q
		getline(myfile, line);
		logq = atol(line.c_str());

		q = power2_ZZ(logq);

		// read 6th line and get slots
		getline(myfile, line);
		slots = atol(line.c_str());

		// read 7th line and get isComplex
		getline(myfile, line);
		isComplex = atoi(line.c_str());

		// read other lines and get ax and bx
		for(long i = 0; i < deg(ax) + 1; i++) {
			getline(myfile, line);
			ax[i] = conv<ZZ>(line.c_str());
		}
		for(long i = 0; i < deg(bx) + 1; i++) {
			getline(myfile, line);
			bx[i] = conv<ZZ>(line.c_str());
		}
		myfile.close();
	} else {
		throw std::invalid_argument("Unable to open file");
	}
}
