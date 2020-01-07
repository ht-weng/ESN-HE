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

template<typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size = 4, int prec = 3)
{
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if(slot_count <= 2 * print_size)
    {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        vec.resize(std::max(vec.size(), 2 * print_size));
        std::cout << "    [";
        for (std::size_t i = 0; i < print_size; i++)
        {
            std::cout << " " << vec[i] << ",";
        }
        if(vec.size() > 2 * print_size)
        {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    std::cout.copyfmt(old_fmt);
}

template<typename T>
inline std::vector<T> rot(std::vector<T> ct, int l)
{
    vector<T> ct_cp = ct;
    if (l >= 0){
        rotate(ct_cp.begin(), ct_cp.begin() + l, ct_cp.end());
    } else {
        rotate(ct_cp.rbegin(), ct_cp.rbegin() - l, ct_cp.rend());
    }
    return ct_cp;
}

template<typename T>
inline std::vector<T> add(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template<typename T>
inline std::vector<T> mult(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::multiplies<T>());
    return result;
}
  
template<typename T>
inline std::vector<T> cmult(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::multiplies<T>());
    return result;
}

inline vector<double> genUs(int d, int k)
{
    vector<double> result;
    result.reserve(d*d);
    for (int i = 0; i < d*d; i++) {
        result.push_back(0);
    }

    for (int l = 0; l < d*d; l++){
        if (k >= 0){
            if ((l-d*k >= 0) && (l-d*k < d-k)){
                result[l] = 1;
            }
        } else {
            if ((l-(d+k)*d >= -k) && (l-(d+k)*d < d)){
                result[l] = 1;
            }
        }   
    }
    return result;
}

inline vector<double> genUt(int d, int k)
{
    vector<double> result;
    result.reserve(d*d);
    for (int i = 0; i < d*d; i++) {
        result.push_back(0);
    }
    for (int i = 0; i < d; i++){
        for (int l = 0; l < d*d; l++){
            if (l == k+d*i){
                result[l] = 1;
            }
        }
    }
    return result;
}

inline vector<double> genVk(int d, int k)
{
    vector<double> result;
    result.reserve(d*d);
    for (int i = 0; i < d*d; i++) {
        result.push_back(0);
    }

    for (int l = 0; l < d*d; l++){
        if ((l % d >= 0) && (l % d < d-k)){
            result[l] = 1;
        }
    }
    return result;
}

inline vector<double> genVkd(int d, int k)
{
    vector<double> result;
    result.reserve(d*d);
    for (int i = 0; i < d*d; i++) {
        result.push_back(0);
    }

    for (int l = 0; l < d*d; l++){
        if ((l % d >= d-k) && (l % d < d)){
            result[l] = 1;
        }
    }
    return result;
}

inline vector<double> linear_tran(vector<double> ct, int d, int U_type, int k_m)
{
    vector<double> result;

    switch (U_type)
    {
        case 0:{
            int i = 0;
            vector<double> ct_ls[2*d-1];
            ct_ls[i] = cmult(ct, genUs(d, 0));
            i++;
            for (int k = 1-d; k < 0; k++){
                ct_ls[i] = add(ct_ls[i-1], cmult(rot(ct, k), genUs(d, k)));
                i++;
            }
            for (int k = 1; k < d; k++){
                ct_ls[i] = add(ct_ls[i-1], cmult(rot(ct, k), genUs(d, k)));
                i++;
            }
            result = ct_ls[i-1];
            break;
        }
        case 1:{
            int i = 0;
            vector<double> ct_ls [2*d-1];
            ct_ls[i] = cmult(ct, genUt(d, 0));
            i++;
            for (int k = 1; k < d; k++){
                ct_ls[i] = add(ct_ls[i-1], cmult(rot(ct, d*k), genUt(d, k)));
                i++;
            }
            result = ct_ls[i-1];
            break;
        }
        case 2:{
            result = add(cmult(rot(ct, k_m), genVk(d, k_m)), cmult(rot(ct, k_m-d), genVkd(d, k_m)));
            break;
        }
        case 3:{
            result = rot(ct, d*k_m);
            break;
        }
        default:{
            cout << "Error: Wrong U_type value!" << endl;
            break;
        }
    }
    return result;
}

inline vector<double> mat_mult(vector<double> ct_a, vector<double> ct_b, int d)
{
    vector<double> result;
    vector<double> ct_a0;
    vector<double> ct_b0;

    ct_a0 = linear_tran(ct_a, d, 0, 0);
    ct_b0 = linear_tran(ct_b, d, 1, 0);

    int i = 0;
    vector<double> ab_ls[d];
    ab_ls[i] = mult(ct_a0, ct_b0);
    i++;

    for (int k = 1; k < d; k++){
        vector<double> ct_ak;
        vector<double> ct_bk;

        ct_ak = linear_tran(ct_a0, d, 2, k);
        ct_bk = linear_tran(ct_b0, d, 3, k);

        ab_ls[i] = add(ab_ls[i-1], mult(ct_ak, ct_bk));
        i++;
    }
    result = ab_ls[i-1];
    return result;
}

// rmat_mult can only be applied when matrix dimension equals to 64
inline vector<double> rmat_mult(vector<double> ct_a, vector<double> ct_b, int d, int l)
{
    vector<double> result;
    vector<double> ct_a0;
    vector<double> ct_b0;

    ct_a0 = linear_tran(ct_a, d, 0, 0);
    ct_b0 = linear_tran(ct_b, d, 1, 0);

    int i = 0;
    vector<double> ab_ls[d];
    ab_ls[i] = mult(ct_a0, ct_b0);
    i++;

    for (int k = 1; k < l; k++){
        vector<double> ct_ak;
        vector<double> ct_bk;

        ct_ak = linear_tran(ct_a0, d, 2, k);
        ct_bk = linear_tran(ct_b0, d, 3, k);

        ab_ls[i] = add(ab_ls[i-1], mult(ct_ak, ct_bk));
        i++;
    }

    vector<double> ab;
    ab = ab_ls[i-1];
    for (int k = 0; k < ceil(log2(d/l)); k++){
        ab_ls[i] = add(ab_ls[i-1], rot(ab, l*d*pow(2, k)));
        i++;
    }
    result = ab_ls[i-1];
    return result;
}

int main() {
    int max_d = 64;
    int mat_d = 3;
    int rmat_d = 1;
    vector<double> ct_a;
    double elem_a = 0.0;
    for (int i = 0; i < max_d*max_d; i++) {
        if ((i%max_d < mat_d) && i < max_d*mat_d) {
            ct_a.push_back(elem_a);
            elem_a++;
        }
        else{
            ct_a.push_back(0.0);
        }
    }

    vector<double> ct_ar;
    for (int i = 0; i < max_d*max_d; i++) {
        if ((i%max_d < mat_d) && i < max_d*mat_d) {
            ct_ar.push_back(i%max_d);
        }
        else{
            ct_ar.push_back(0.0);
        }
    }

    vector<double> ct_b;
    double elem_b = 0.0;
    for (int i = 0; i < max_d*max_d; i++) {
        if (i%max_d < mat_d && i < max_d*mat_d) {
            ct_b.push_back(elem_b);
            elem_b++;
        }
        else{
            ct_b.push_back(0.0);
        }
    }

    cout<<"Square Matrix Multiplication"<<endl;
    cout<<"A: "<<endl;
    print_vector(ct_a, mat_d*max_d);
    cout<<"B: "<<endl;
    print_vector(ct_b, mat_d*max_d);
    cout<<"A * B: "<<endl;
    vector<double> result_1 = mat_mult(ct_a, ct_b, max_d);
    print_vector(result_1, mat_d*max_d);

    cout<<"Rectangular Matrix Multiplication"<<endl;
    cout<<"AR: "<<endl;
    print_vector(ct_ar, mat_d*max_d);
    cout<<"B: "<<endl;
    print_vector(ct_b, mat_d*max_d);
    cout<<"AR * B: "<<endl;
    // rmat_mult is not used here because the matrix dimension is only 3
    vector<double> result_2 = mat_mult(ct_ar, ct_b, max_d);
    print_vector(result_2, mat_d*max_d);

    return 0;  
    }  
