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

inline vector<double> ema(const vector<double>& s, int n) {
    vector<double> ema;
    int j = 1;

    // get n sma first and calculate the next n period ema
    vector<double> s_head = slice(s, 0, n);
    double sma = sum(s_head)/n;
    double multiplier = 2/(1+n);
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

int main() {

    vector<double> prices;

    for (int i = 0; i < 100; i++ ) {
        prices.push_back(i);
    }

    print_vector(ema(prices, 5));

    return 0;
}