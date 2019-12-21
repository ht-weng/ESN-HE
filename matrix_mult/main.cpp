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
#include "seal/seal.h"

using namespace std;
using namespace seal;

//********************************************************************************
// SEAL helper functions
//********************************************************************************
/*
Helper function: Prints the name of the example in a fancy banner.
*/
inline void print_example_banner(std::string title)
{
    if (!title.empty())
    {
        std::size_t title_length = title.length();
        std::size_t banner_length = title_length + 2 * 10;
        std::string banner_top = "+" + std::string(banner_length - 2, '-') + "+";
        std::string banner_middle =
            "|" + std::string(9, ' ') + title + std::string(9, ' ') + "|";

        std::cout << std::endl
            << banner_top << std::endl
            << banner_middle << std::endl
            << banner_top << std::endl;
    }
}

/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(std::shared_ptr<seal::SEALContext> context)
{
    // Verify parameters
    if (!context)
    {
        throw std::invalid_argument("context is not set");
    }
    auto &context_data = *context->key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::BFV:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::CKKS:
        scheme_name = "CKKS";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " <<
        context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_mod_count = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_mod_count - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::BFV)
    {
        std::cout << "|   plain_modulus: " << context_data.
            parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

/*
Helper function: Prints the `parms_id' to std::ostream.
*/
inline std::ostream &operator <<(std::ostream &stream, seal::parms_id_type parms_id)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    stream << std::hex << std::setfill('0')
        << std::setw(16) << parms_id[0] << " "
        << std::setw(16) << parms_id[1] << " "
        << std::setw(16) << parms_id[2] << " "
        << std::setw(16) << parms_id[3] << " ";

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);

    return stream;
}

/*
Helper function: Prints a vector of floating-point values.
*/
template<typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size = 4, int prec = 3)
{
    /*
    Save the formatting information for std::cout.
    */
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

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}


/*
Helper function: Prints a matrix of values.
*/
template<typename T>
inline void print_matrix(std::vector<T> matrix, std::size_t row_size)
{
    /*
    We're not going to print every column of the matrix (there are 2048). Instead
    print this many slots from beginning and end of the matrix.
    */
    std::size_t print_size = 5;

    std::cout << std::endl;
    std::cout << "    [";
    for (std::size_t i = 0; i < print_size; i++)
    {
        std::cout << std::setw(3) << std::right << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = row_size - print_size; i < row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i]
            << ((i != row_size - 1) ? "," : " ]\n");
    }
    std::cout << "    [";
    for (std::size_t i = row_size; i < row_size + print_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = 2 * row_size - print_size; i < 2 * row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i]
            << ((i != 2 * row_size - 1) ? "," : " ]\n");
    }
    std::cout << std::endl;
};

/*
Helper function: Print line number.
*/
inline void print_line(int line_number)
{
    std::cout << "Line " << std::setw(3) << line_number << " --> ";
}

//********************************************************************************
// Functions for linear transformation
//********************************************************************************
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
            vector<double> ct_ls [2*d-1];
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

//********************************************************************************
// Matrix multiplication functions
//********************************************************************************
inline vector<double> mat_mult(vector<double> ct_a, vector<double> ct_b, int d)
{
    vector<double> result;
    vector<double> ct_a0;
    vector<double> ct_b0;

    ct_a0 = linear_tran(ct_a, d, 0, 0);
    ct_b0 = linear_tran(ct_b, d, 1, 0);

    int i = 0;
    vector<double> ab_ls [d];
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

inline vector<double> rmat_mult(vector<double> ct_a, vector<double> ct_b, int d, int l)
{
    vector<double> result;
    vector<double> ct_a0;
    vector<double> ct_b0;

    ct_a0 = linear_tran(ct_a, d, 0, 0);
    ct_b0 = linear_tran(ct_b, d, 1, 0);

    int i = 0;
    vector<double> ab_ls [d];
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
//********************************************************************************
// Main function
//********************************************************************************
int main()
{
    //////////////////////////////////////////////////////////////////////////////
    // SEAL settings
    //////////////////////////////////////////////////////////////////////////////
    // Set up the CKKS scheme.
    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 4096;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    // TODO: Modify poly_modulus_degree, coeff_modulus and scale
    parms.set_coeff_modulus(CoeffModulus::Create(
        poly_modulus_degree, { 60, 40, 40, 60 }));
    double scale = pow(2.0, 40);

    // Set up Context
    auto context = SEALContext::Create(parms);
    print_parameters(context);
    cout << endl;

    // Generate public and private keys
    KeyGenerator keygen(context);
    auto public_key = keygen.public_key();
    auto secret_key = keygen.secret_key();
    auto relin_keys = keygen.relin_keys();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << endl;
    cout << endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // SEAL Example
    //////////////////////////////////////////////////////////////////////////////
    // Average of mackey glass data
    // Initialise the sum variable
    Plaintext x_sum;
    encoder.encode(0, scale, x_sum);
    Ciphertext x_sum_encrypted;
    encryptor.encrypt(x_sum, x_sum_encrypted);

    // Loop on all time series values
    for (size_t i = 0; i < sample_n; i++) {
        // Define input plain text
        vector<double> input;
        input.reserve(slot_count);
        input.push_back(X[i]);

        // Encode and encrypt input vector
        Plaintext x_plain;
        encoder.encode(input, scale, x_plain);
        Ciphertext x_encrypted;
        encryptor.encrypt(x_plain, x_encrypted);

        // Sum
        evaluator.add_inplace(x_sum_encrypted, x_encrypted);
    }

    // Divide to get average
    Plaintext denom;
    encoder.encode(0.001, scale, denom);
    evaluator.multiply_plain_inplace(x_sum_encrypted, denom);

    // Decrypt, decode, and print the result
    Plaintext plain_result;
    decryptor.decrypt(x_sum_encrypted, plain_result);
    vector<double> result;
    encoder.decode(plain_result, result);
    cout << "Result vector: " << endl;
    print_vector(result, 10, 7);
    cout << endl;
    cout << "Correct average value: " << X_avg << endl;
    cout << endl;

    return 0;
}
