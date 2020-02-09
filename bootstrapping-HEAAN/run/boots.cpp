/*
* Test bootstrapping's effect on precision
*/
#include <cstddef>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <stdlib.h>
#include "../src/HEAAN.h"

using namespace std;
using namespace NTL;

void TestScheme::testBootstrap(long logq, long logp, long logSlots, long logT) {
	cout << "!!! START TEST BOOTSTRAP !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	timeutils.start("Key generating");
	scheme.addBootKey(secretKey, logSlots, logq + 4);
	timeutils.stop("Key generated");

	long slots = (1 << logSlots);
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(slots);

	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, slots, logp, logq);

	cout << "cipher logq before: " << cipher.logq << endl;

	scheme.modDownToAndEqual(cipher, logq);
	scheme.normalizeAndEqual(cipher);
	cipher.logq = logQ;
	cipher.logp = logq + 4;

	Ciphertext rot;
	timeutils.start("SubSum");
	for (long i = logSlots; i < logNh; ++i) {
		scheme.leftRotateFast(rot, cipher, (1 << i));
		scheme.addAndEqual(cipher, rot);
	}
	scheme.divByPo2AndEqual(cipher, logNh);
	timeutils.stop("SubSum");

	timeutils.start("CoeffToSlot");
	scheme.coeffToSlotAndEqual(cipher);
	timeutils.stop("CoeffToSlot");

	timeutils.start("EvalExp");
	scheme.evalExpAndEqual(cipher, logT);
	timeutils.stop("EvalExp");

	timeutils.start("SlotToCoeff");
	scheme.slotToCoeffAndEqual(cipher);
	timeutils.stop("SlotToCoeff");

	cipher.logp = logp;
	cout << "cipher logq after: " << cipher.logq << endl;

	complex<double>* dvec = scheme.decrypt(secretKey, cipher);

	StringUtils::compare(mvec, dvec, slots, "boot");

	cout << "!!! END TEST BOOTSRTAP !!!" << endl;
}

void TestScheme::testBootstrapSingleReal(long logq, long logp, long logT) {
	cout << "!!! START TEST BOOTSTRAP SINGLE REAL !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	timeutils.start("Key generating");
	scheme.addBootKey(secretKey, 0, logq + 4);
	timeutils.stop("Key generated");

	double mval = EvaluatorUtils::randomReal();

	Ciphertext cipher;
	scheme.encryptSingle(cipher, mval, logp, logq);

	cout << "cipher logq before: " << cipher.logq << endl;
	scheme.modDownToAndEqual(cipher, logq);
	scheme.normalizeAndEqual(cipher);
	cipher.logq = logQ;

	Ciphertext rot, cconj;
	timeutils.start("SubSum");
	for (long i = 0; i < logNh; ++i) {
		scheme.leftRotateFast(rot, cipher, 1 << i);
		scheme.addAndEqual(cipher, rot);
	}
	scheme.conjugate(cconj, cipher);
	scheme.addAndEqual(cipher, cconj);
	scheme.divByPo2AndEqual(cipher, logN);
	timeutils.stop("SubSum");

	timeutils.start("EvalExp");
	scheme.evalExpAndEqual(cipher, logT);
	timeutils.stop("EvalExp");

	cout << "cipher logq after: " << cipher.logq << endl;

	cipher.logp = logp;
	complex<double> dval = scheme.decryptSingle(secretKey, cipher);

	StringUtils::compare(mval, dval.real(), "boot");

	cout << "!!! END TEST BOOTSRTAP SINGLE REAL !!!" << endl;
}

//********************************************************************************
// Main function
//********************************************************************************
int main() {

	long logq = 600; // Ciphertext polynomial modulus
	long logp = 30; // Scale

    // long logq = 1200; ///< Ciphertext Modulus
	// long logp = 30; ///< Real message will be quantized by multiplying 2^40
	// long logn = 5; ///< log2(The number of slots)

    logq = logp + 10; //< suppose the input ciphertext of bootstrapping has logq = logp + 10
    long logn = 3; //< larger logn will make bootstrapping tech much slower
    long logT = 4; //< this means that we use Taylor approximation in [-1/T,1/T] with double angle fomula

    TestScheme::testBootstrap(logq, logp, logn, logT);

	return 0;
}
