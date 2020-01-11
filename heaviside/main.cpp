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
// main function for homomorphic encryption of mackey glass equation
int main()
{
    //****************************************************************************
    // SEAL settings
    // Set up the CKKS scheme.
    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(
        poly_modulus_degree, { 40, 30, 30, 30, 30, 40 }));
    double scale = pow(2.0, 30);

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

    // Construct an array to store parameters ids of level 1 to level 5
    int p_i = 0;
    parms_id_type parms_ids[5];
    auto context_data = context->first_context_data();
    while (context_data)
    {
        parms_ids[p_i] = context_data->parms_id();
        context_data = context_data->next_context_data();
        p_i++;
    }

    //****************************************************************************
    // Approximation of sigmoid function
    //****************************************************************************
    // Set up input data
    vector<double> input;
    input.reserve(slot_count);
    double x = 0.99;
    input.push_back(x);
    cout << "Input data x: " << endl;
    print_vector(input, 1, 6);

    cout << "Evaluating polynomial f(x) = 0.375*x^5 - 1.25*x^3 + 1.875*x to approximate Heaviside function" << endl;

    // Set up coefficients
    Plaintext coeff1_plain, coeff2_plain, coeff3_plain;
    encoder.encode(0.375, scale, coeff1_plain);
    encoder.encode((-1.25), scale, coeff2_plain);
    encoder.encode(1.875, scale, coeff3_plain);

    // Encode and encrypt input data
    Plaintext x_plain;
    encoder.encode(input, scale, x_plain);
    // x is at level 0
    Ciphertext x_encrypted;
    encryptor.encrypt(x_plain, x_encrypted);

    // Calculate x^2, which is at level 1
    Ciphertext x2_encrypted;
    evaluator.square(x_encrypted, x2_encrypted);
    evaluator.relinearize_inplace(x2_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x2_encrypted);
    x2_encrypted.scale() = scale;

    // Calculate x^3, which is at level 2
    Ciphertext x3_encrypted;
    // Switch mod to ensure parameters id match
    evaluator.mod_switch_to_inplace(x_encrypted, parms_ids[1]);
    evaluator.multiply(x_encrypted, x2_encrypted, x3_encrypted);
    evaluator.relinearize_inplace(x3_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x3_encrypted);
    x3_encrypted.scale() = scale;

    // Calculate x^5, which is at level 3
    Ciphertext x5_encrypted;
    // Switch mod to ensure parameters id match
    evaluator.mod_switch_to_inplace(x2_encrypted, parms_ids[2]);
    evaluator.multiply(x2_encrypted, x3_encrypted, x5_encrypted);
    evaluator.relinearize_inplace(x5_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x5_encrypted);
    x5_encrypted.scale() = scale;

    // Calculate 0.375*x^5, which is at level 4
    // Switch mod to ensure parameters id match
    evaluator.mod_switch_to_inplace(coeff1_plain, parms_ids[3]);
    evaluator.multiply_plain_inplace(x5_encrypted, coeff1_plain);
    evaluator.rescale_to_next_inplace(x5_encrypted);
    x5_encrypted.scale() = scale;

    // Calculate -1.25*x^3, which is at level 3
    // Switch mod to ensure parameters id match
    evaluator.mod_switch_to_inplace(coeff2_plain, parms_ids[2]);
    evaluator.multiply_plain_inplace(x3_encrypted, coeff2_plain);
    evaluator.rescale_to_next_inplace(x3_encrypted);
    x3_encrypted.scale() = scale;

    // Calculate 1.875*x, which is at level 2
    // Switch mod to ensure parameters id match
    // Note that previously we haved switched x_encrypted' level to 1
    evaluator.mod_switch_to_inplace(coeff3_plain, parms_ids[1]);
    evaluator.multiply_plain_inplace(x_encrypted, coeff3_plain);
    evaluator.rescale_to_next_inplace(x_encrypted);
    x_encrypted.scale() = scale;

    // To do the addition, we have to ensure that the terms have same parms_id and scale
    evaluator.mod_switch_to_inplace(x3_encrypted, parms_ids[4]);
    evaluator.mod_switch_to_inplace(x_encrypted, parms_ids[4]);

    Ciphertext result_encrypted;
    evaluator.add(x5_encrypted, x3_encrypted, result_encrypted);
    evaluator.add_inplace(result_encrypted, x_encrypted);
    Plaintext result_plain;
    decryptor.decrypt(result_encrypted, result_plain);
    vector<double> result;
    encoder.decode(result_plain, result);
    cout << "Evaluated Result: " << endl;
    print_vector(result, 1, 6);

    cout << "True Result Without Homomorphic Encryption: " << endl;
    cout << 0.375*pow(x, 5) - 1.25*pow(x, 3) + 1.875*x << endl;

    return 0;
}
