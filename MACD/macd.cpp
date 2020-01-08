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

using namespace std;

inline double sum(const vector<double>& v) {
    double result;

    for(int i=0; i < v.size(); i++) {
        result += v[i];
    }
    return result;
}

inline vector<double> slice(const vector<double>& v, int start=0, int end=-1) {
    int oldlen = v.size();
    int newlen;

    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }

    vector<double> nv(newlen);

    for (int i=0; i<newlen; i++) {
        nv[i] = v[start+i];
    }
    return nv;
}

inline void print_vector(const vector<double>& v) {
    for(int i=0; i < v.size(); i++) {
        cout << v[i] << ' ';
    }
    cout << "\n";
}

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

inline vector<double> ema(const vector<double>& s, int n) {
    vector<double> ema;
    int j = 1;

    // get n sma first and calculate the next n period ema
    vector<double> s_head = slice(s, 0, n);
    double sma = sum(s_head)/n;
    double multiplier = 2.0/(1.0+n);
    ema.push_back(sma);

    // EMA(current) = ( (Price(current) - EMA(prev) ) x Multiplier) + EMA(prev)
    ema.push_back(((s[n] - sma) * multiplier) + sma);

    // now calculate the rest of the values
    for (int i = n+1; i < s.size(); i++) {
        double tmp;
        tmp = (s[i] - ema[j]) * multiplier + ema[j];
        j++;
        ema.push_back(tmp);
    }
    return ema;
}

inline vector<double> decision(const vector<double>& s) {
    vector<double> d;
    d.push_back(0);
    for (int i = 1; i < s.size(); i++) {
        if (s[i]*s[i-1] < 0) {
            if (s[i] < 0) {
                // sell
                d.push_back(-1);
            } else {
                // buy
                d.push_back(1);
            }
        } else if (s[i]*s[i-1] == 0) {
            if (s[i-1] < 0 || s[i] > 0) {
                // buy
                d.push_back(1);
            } else if (s[i-1] > 0 || s[i] < 0) {
                //sell
                d.push_back(-1);
            } else {
                // do nothing
                d.push_back(0);
            }
        } else {
            // do nothing
            d.push_back(0);
        }
    }
    return d;
}

inline vector<double> tanh(const vector<double>& s) {
    vector<double> result;
    for (int i = 0; i < s.size(); i++) {
        double x = s[i];
        // f1(x) = -0.5*x^3 + 1.5*x
        // f2(x) = 0.375*x^5 - 1.25*x^3 + 1.875*x
        // f3(x) = −0.3125*x^7 + 1.3125*x^5 - 2.1875*x^3 + 2.1875*x
        // f4(x)= 0.2734375*x^9 - 1.40625*x^7 + 2.953125*x^5 - 3.28125*x^3 + 2.46093758*x
        // result.push_back(0.5*pow(x, 3) + 1.5*x);
        result.push_back(0.375*pow(x, 5) - 1.25*pow(x, 3) + 1.875*x);
    }
    return result;
}

int main() {

    vector<double> prices = csv2vec("apple_prices.csv");

    vector<double> ema12 = ema(prices, 12);
    vector<double> ema26 = ema(prices, 26);

    // slice ema12 to make sure its size matches the size of ema26
    vector<double> ema12_sliced = slice(ema12, 14, ema12.size());

    // calculate macd
    vector<double> ema_diff;
    for (int i = 0; i < ema26.size(); i++) {
        ema_diff.push_back(ema12_sliced[i]-ema26[i]);
    }
    
    vector<double> ema9 = ema(ema_diff, 9);
    vector<double> ema_diff_sliced = slice(ema_diff, 8, ema_diff.size());

    vector<double> macd;
    for (int i = 0; i < ema9.size(); i++) {
        macd.push_back(ema_diff_sliced[i]-ema9[i]);
    }

    vector<double> decisions = decision(macd);
    vector<double> tanh_decisions = tanh(macd);

    // write data to csv for plotting
    ofstream output_file1("macd.csv");
    ostream_iterator<double> output_iterator1(output_file1, "\n");
    copy(macd.begin(), macd.end(), output_iterator1);

    ofstream output_file2("decisions.csv");
    ostream_iterator<double> output_iterator2(output_file2, "\n");
    copy(decisions.begin(), decisions.end(), output_iterator2);

    // write decisions to csv
    ofstream output_file3("tanh_decisions.csv");
    ostream_iterator<double> output_iterator3(output_file3, "\n");
    copy(tanh_decisions.begin(), tanh_decisions.end(), output_iterator3);

    return 0;
}