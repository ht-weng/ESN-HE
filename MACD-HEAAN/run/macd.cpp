/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
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

inline vector<double> csv2vec(string inputFileName) {
    vector<double> data;
    ifstream inputFile(inputFileName);
    int l = 0;

    while (inputFile) {
        l++;
        string line;
        if (!getline(inputFile, line)) {
            break;
        }
        try {
            data.push_back(stof(line));
        }
        catch (const std::invalid_argument e) {
            cout << "NaN found in file " << inputFileName << " line " << l
                 << endl;
            e.what();
        }
    }

    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }

    return data;
}

inline void printVector(vector<double>& v) {
    for(int i=0; i < v.size(); i++) {
        cout << v[i] << ' ';
    }
    cout << "\n";
}

inline Ciphertext sum(vector<Ciphertext>& v, Scheme &scheme) {
    Ciphertext result_encrypted = v[0];

    for(int i = 1; i < v.size(); i++) {
        Ciphertext val = v[i];
        scheme.addAndEqual(result_encrypted, val);
    }
    return result_encrypted;
}

inline vector<Ciphertext> slice(vector<Ciphertext>& v, int start , int end) {
    int oldlen = v.size();
    int newlen;

    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }

    vector<Ciphertext> nv;

    for (int i = 0; i < newlen; i++) {
        nv.push_back(v[start+i]);
    }
    return nv;
}

inline Ciphertext getSample(int time, long logq, long logp, long logn, Scheme &scheme) {
	vector<double> prices;
	prices = csv2vec("apple_prices.csv");

    complex<double> sample;
	sample.real(prices[time]);

    Ciphertext sample_encrypted;
	scheme.encryptSingle(sample_encrypted, sample, logp, logq);
    return sample_encrypted;
}

inline void assembleSample(Ciphertext sample, vector<Ciphertext>& past_prices) {
    past_prices.push_back(sample);
}

inline vector<Ciphertext> wma(vector<Ciphertext>& s, int m, Scheme &scheme, long logq, long logp, long logn) {
    vector<Ciphertext> wma;
    vector<double> weights;
    // generate a list of weights of the window size
    for (int i = 0; i < m; i++) {
        double w = (2.0*(i+1.0)/(m*(m+1.0)));
        weights.push_back(w);
    }
    // multiply corresponding data points and weights to get WMA
    for (int i = 0; i < s.size()-m; i++) {
        vector<Ciphertext> s_sliced;
        vector<Ciphertext> window;
        s_sliced = slice(s, i, i+m);
        for (int j = 0; j < m; j++) {
			Ciphertext tmp;
            tmp.copy(s_sliced[j]);
            complex<double> tmp_weight;
			tmp_weight.real(weights[j]);
            // scheme.modDownToAndEqual();
            scheme.multByConstAndEqual(tmp, tmp_weight, logp);
			scheme.reScaleByAndEqual(tmp, logp);
			
            window.push_back(tmp);
        }
        Ciphertext sumup = sum(window, scheme);
        wma.push_back(sumup);
    }
    return wma;
}

inline void decryptCipherVec(vector<Ciphertext>& v, Scheme &scheme, SecretKey &secretKey) {
    // vector<double> res;
    for (int i = 0; i < v.size(); i++) {
        complex<double> val;
        Ciphertext val_encrypted;
        val_encrypted = v[i];
        val = scheme.decryptSingle(secretKey, val_encrypted);
        cout << val << endl;
        // res.push_back(real(val));
    }
    // return res;
}

int main() {

	long logq = 600; ///< Ciphertext Modulus
	long logp = 30; ///< Real message will be quantized by multiplying 2^40
	long logn = 10; ///< log2(The number of slots)
	int time_max = 100;	

	SetNumThreads(8);
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = pow(2, logn);

	cout << "Data Import Starts" << endl;
	cout << "\n";

    vector<Ciphertext> prices_encrypted;
    for (int i = 0; i < time_max; i++) {
        Ciphertext sample = getSample(i, logq, logp, logn, scheme);
        assembleSample(sample, prices_encrypted);
    }

	cout << "Data Import Complete" << endl;
	cout << "\n";
	cout << "MACD Analysis Starts" << endl;
	cout << "\n";

	vector<Ciphertext> wma12_encrypted = wma(prices_encrypted, 12, scheme, logq, logp, logn);
	// vector<Ciphertext> wma12_encrypted_sliced = slice(wma12_encrypted, 14, wma12_encrypted.size());

	cout << "wma12 done" << endl;

	// vector<Ciphertext> wma26_encrypted;
    // wma26_encrypted = wma(prices_encrypted, 26, scheme, logq, logp, logn);

	// cout << "wma26 done" << endl;

	// vector<Ciphertext> wma_diff_encrypted;
    // for (int i = 0; i < wma26_encrypted.size(); i++) {
    //     Ciphertext tmp_diff;
    //     scheme.sub(tmp_diff, wma12_encrypted_sliced[i], wma26_encrypted[i]);
    //     wma_diff_encrypted.push_back(tmp_diff);
    // }

	// cout << "wma diff done" << endl;

    // vector<Ciphertext> wma9_encrypted;
    // wma9_encrypted = wma(wma_diff_encrypted, 9, scheme, logq, logp, logn);

	// cout << "wma9 done" << endl;

    // vector<Ciphertext> wma_diff_sliced = slice(wma_diff_encrypted, 9, wma_diff_encrypted.size());

    // vector<Ciphertext> macd_encrypted;
    // for (int i = 0; i < wma9_encrypted.size(); i++) {
    //     Ciphertext tmp_macd;
    //     scheme.sub(tmp_macd, wma_diff_sliced[i], wma9_encrypted[i]);
    //     macd_encrypted.push_back(tmp_macd);
    // }

	// cout << "MACD Analysis Complete" << endl;
	// cout << "\n";
	// cout << "Data Export Starts" << endl;
	// cout << "\n";

	vector<double> wma12;
    decryptCipherVec(wma12_encrypted, scheme, secretKey);
    
    // for (int i = 0; i < wma12_encrypted.size(); i++) {
    //     complex<double> val;
    //     Ciphertext val_encrypted;
    //     val_encrypted = wma12_encrypted[i];
    //     val = scheme.decryptSingle(secretKey, val_encrypted);
    //     wma12.push_back(real(val));
    // }

	// cout << "wma12 exported" << endl;

    // vector<double> wma26;
    // for (int i = 0; i < wma26_encrypted.size(); i++) {
    //     complex<double> val;
    //     Ciphertext val_encrypted;
    //     val_encrypted = wma26_encrypted[i];
    //     val = scheme.decryptSingle(secretKey, val_encrypted);
    //     wma26.push_back(real(val));
    // }

	// cout << "wma26 exported" << endl;

    // vector<double> wma_diff;
    // for (int i = 0; i < wma_diff_sliced.size(); i++) {
    //     complex<double> val;
    //     Ciphertext val_encrypted;
    //     val_encrypted = wma_diff_sliced[i];
    //     val = scheme.decryptSingle(secretKey, val_encrypted);
    //     wma_diff.push_back(real(val));
    // }

	// cout << "wma diff exported" << endl;

    // vector<double> wma9;
    // for (int i = 0; i < wma9_encrypted.size(); i++) {
	// 	cout << i << endl;
    //     complex<double> val;
    //     Ciphertext val_encrypted;
    //     val_encrypted = wma9_encrypted[i];
    //     val = scheme.decryptSingle(secretKey, val_encrypted);
    //     wma9.push_back(real(val));
    // }

	// cout << "wma9 exported" << endl;

    // vector<double> macd;
    // for (int i = 0; i < macd_encrypted.size(); i++) {
    //     complex<double> val;
    //     Ciphertext val_encrypted;
    //     val_encrypted = macd_encrypted[i];
    //     val = scheme.decryptSingle(secretKey, val_encrypted);
    //     macd.push_back(real(val));
    // }

	// cout << "macd exported" << endl;

	// printVector(wma12);
	// printVector(wma26);
	// printVector(wma_diff);
	// printVector(macd);

	// vector<Ciphertext> decisions_encrypted;
    // for (int i = 1; i < macd_encrypted.size(); i++) {

    //     // First store the latest two macd signals m(t) and m(t-1)
    //     Ciphertext mt = macd_encrypted[i];
    //     Ciphertext mt_1 = macd_encrypted[i-1];

    //     complex<double> coeff_norm;
	// 	coeff_norm.real(0.25);

    //     scheme.multByConstAndEqual(mt, coeff_norm, logp);
    //     scheme.reScaleByAndEqual(mt, logq);
	// 	scheme.multByConstAndEqual(mt_1, coeff_norm, logp);
    //     scheme.reScaleByAndEqual(mt_1, logq);

    //     Ciphertext recent_diff;
    //     scheme.sub(recent_diff, mt, mt_1);
    //     Ciphertext recent_mult;
	// 	scheme.mult(recent_mult, mt, mt_1);
	// 	scheme.reScaleByAndEqual(recent_mult, logq);

    //     // Set up coefficients
    //     complex<double> coeff1, coeff2, coeff3, offset;
	// 	coeff1.real(0.375);
	// 	coeff2.real(-1.25);
	// 	coeff3.real(1.875);
	// 	offset.real(1.0);

    //     Ciphertext x2_encrypted;
    //     scheme.square(x2_encrypted, recent_mult);
	// 	scheme.reScaleByAndEqual(x2_encrypted, logq);

    //     // Calculate x^3, which is at level 6
    //     Ciphertext x3_encrypted;
	// 	scheme.mult(x3_encrypted, recent_mult, x2_encrypted);
	// 	scheme.reScaleByAndEqual(x3_encrypted, logq);

    //     // Calculate x^5, which is at level 7
    //     Ciphertext x5_encrypted;
	// 	scheme.mult(x5_encrypted, x2_encrypted, x3_encrypted);
	// 	scheme.reScaleByAndEqual(x5_encrypted, logq);

	// 	scheme.multByConstAndEqual(x5_encrypted, coeff1, logp);
	// 	scheme.reScaleByAndEqual(x5_encrypted, logq);

	// 	scheme.multByConstAndEqual(x3_encrypted, coeff2, logp);
	// 	scheme.reScaleByAndEqual(x3_encrypted, logq);

	// 	scheme.multByConstAndEqual(x5_encrypted, coeff1, logp);
	// 	scheme.reScaleByAndEqual(x5_encrypted, logq);

	// 	scheme.multByConstAndEqual(recent_mult, coeff3, logp);
	// 	scheme.reScaleByAndEqual(recent_mult, logq);

	// 	Ciphertext sign_encrypted;
	// 	scheme.add(sign_encrypted, x5_encrypted, x3_encrypted);
	// 	scheme.addAndEqual(sign_encrypted, recent_mult);
	// 	scheme.addConstAndEqual(sign_encrypted, offset, logp);

    //     Ciphertext result_encrypted;
	// 	scheme.mult(result_encrypted, sign_encrypted, recent_diff);
	// 	scheme.reScaleByAndEqual(result_encrypted, logq);

    //     decisions_encrypted.push_back(result_encrypted);
    // }

    // vector<double> decisions;
    // for (int i = 0; i < decisions_encrypted.size(); i++) {
    //     complex<double> result;
    //     Plaintext result_plain; 
    //     result = scheme.decryptSingle(secretKey, decisions_encrypted[i]);
    //     decisions.push_back(real(result));
    // }

    // cout << "Decisions calculation done" << endl;

    // // write data to csv for plotting
    // ofstream output_file1("macd_heaan.csv");
    // ostream_iterator<double> output_iterator1(output_file1, "\n");
    // copy(macd.begin(), macd.end(), output_iterator1);

    // ofstream output_file2("wma12_heaan.csv");
    // ostream_iterator<double> output_iterator2(output_file2, "\n");
    // copy(wma12.begin(), wma12.end(), output_iterator2);

    // ofstream output_file3("wma26_heaan.csv");
    // ostream_iterator<double> output_iterator3(output_file3, "\n");
    // copy(wma26.begin(), wma26.end(), output_iterator3);

    // ofstream output_file4("wma_diff_heaan.csv");
    // ostream_iterator<double> output_iterator4(output_file4, "\n");
    // copy(wma_diff.begin(), wma_diff.end(), output_iterator4);

    // ofstream output_file5("wma9_heaan.csv");
    // ostream_iterator<double> output_iterator5(output_file5, "\n");
    // copy(wma9.begin(), wma9.end(), output_iterator5);

	// ofstream output_file6("decisions_heaan.csv");
    // ostream_iterator<double> output_iterator6(output_file6, "\n");
    // copy(decisions.begin(), decisions.end(), output_iterator6);

	// cout << "Data Output Complete" << endl;

	return 0;
}
