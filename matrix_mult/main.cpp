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

// ********************************************************************************
// Functions for linear transformation
// ********************************************************************************
inline Plaintext genUs(int d, int k, double scale, CKKSEncoder &encoder)
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
    Plaintext result_plain;
    encoder.encode(result, scale, result_plain);
    return result_plain;
}

inline Plaintext genUt(int d, int k, double scale, CKKSEncoder &encoder)
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
    Plaintext result_plain;
    encoder.encode(result, scale, result_plain);
    return result_plain;
}

inline Plaintext genVk(int d, int k, double scale, CKKSEncoder &encoder)
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
    Plaintext result_plain;
    encoder.encode(result, scale, result_plain);
    return result_plain;
}

inline Plaintext genVkd(int d, int k, double scale, CKKSEncoder &encoder)
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
    Plaintext result_plain;
    encoder.encode(result, scale, result_plain);
    return result_plain;
}

inline Plaintext genOnes(int d, double scale, CKKSEncoder &encoder)
{
    vector<double> result;
    result.reserve(d*d);
    for (int i = 0; i < d*d; i++) {
        result.push_back(1);
    }

    Plaintext result_plain;
    encoder.encode(result, scale, result_plain);
    return result_plain;
}

inline Ciphertext linear_tran(Ciphertext &ct, int d, int U_type, int k_m, double scale, 
    CKKSEncoder &encoder, Evaluator &evaluator, const GaloisKeys &gal_keys, 
    const parms_id_type& parms_id)
{
    Ciphertext result;
    
    switch (U_type)
    {
        case 0:{
            int i = 0;
            Ciphertext ct_ls[2*d-1];
            Plaintext mat_Us;
            ct_ls[i] = ct;
            mat_Us = genUs(d, 0, scale, encoder);
            evaluator.mod_switch_to_inplace(mat_Us, parms_id);

            evaluator.multiply_plain_inplace(ct_ls[i], mat_Us);
            // cout << "ct_ls[i]: " << ct_ls[i].parms_id() << endl;
            evaluator.rescale_to_next_inplace(ct_ls[i]);
            ct_ls[i].scale() = scale;
            // evaluator.mod_switch_to_inplace(ct_ls[i], last_parms_id);
            // cout << "ct_ls[i]: " << ct_ls[i].parms_id() << endl;
            i++;
            for (int k = 1-d; k < 0; k++){
                Ciphertext temp_rot;
                mat_Us = genUs(d, k, scale, encoder);
                evaluator.mod_switch_to_inplace(mat_Us, parms_id);
                
                evaluator.rotate_vector(ct, k, gal_keys, temp_rot);
                // cout << "ct: " << ct.parms_id() << endl;
                // cout << "temp_rot: " << temp_rot.parms_id() << endl;
                evaluator.multiply_plain_inplace(temp_rot, mat_Us);
                // cout << "temp_rot: " << temp_rot.parms_id() << endl;
                evaluator.rescale_to_next_inplace(temp_rot);
                temp_rot.scale() = scale;
                // evaluator.mod_switch_to_inplace(temp_rot, last_parms_id);
                // cout << "ct_ls[i]: " << ct_ls[i].parms_id() << endl;
                evaluator.add(ct_ls[i-1], temp_rot, ct_ls[i]);
                // cout << "ct_ls[i]: " << ct_ls[i].parms_id() << endl;
                i++;
            }
            for (int k = 1; k < d; k++){
                Ciphertext temp_rot;
                mat_Us = genUs(d, k, scale, encoder);
                evaluator.mod_switch_to_inplace(mat_Us, parms_id);
                // cout << "ct: " << ct.parms_id() << endl;
                evaluator.rotate_vector(ct, k, gal_keys, temp_rot);
                // cout << "ct: " << ct.parms_id() << endl;
                // cout << "temp_rot: " << temp_rot.parms_id() << endl;
                evaluator.multiply_plain_inplace(temp_rot, mat_Us);
                // cout << "temp_rot: " << temp_rot.parms_id() << endl;
                evaluator.rescale_to_next_inplace(temp_rot);
                temp_rot.scale() = scale;
                // evaluator.mod_switch_to_inplace(temp_rot, last_parms_id);
                // cout << "ct_ls[i]: " << ct_ls[i].parms_id() << endl;
                evaluator.add(ct_ls[i-1], temp_rot, ct_ls[i]);
                // cout << "ct_ls[i]: " << ct_ls[i].parms_id() << endl;
                i++;
            }
            result = ct_ls[i-1];
            break;
        }
        case 1:{
            int i = 0;
            Ciphertext ct_ls[2*d-1];
            Plaintext mat_Ut;
            ct_ls[i] = ct;
            mat_Ut = genUt(d, 0, scale, encoder);
            evaluator.mod_switch_to_inplace(mat_Ut, parms_id);
            evaluator.multiply_plain_inplace(ct_ls[i], mat_Ut);
            evaluator.rescale_to_next_inplace(ct_ls[i]);
            ct_ls[i].scale() = scale;
            // evaluator.mod_switch_to_inplace(ct_ls[i], last_parms_id);
            i++;
            for (int k = 1; k < d; k++){
                Ciphertext temp_rot;
                mat_Ut = genUt(d, k, scale, encoder);
                evaluator.mod_switch_to_inplace(mat_Ut, parms_id);

                evaluator.rotate_vector(ct, d*k, gal_keys, temp_rot);
                evaluator.multiply_plain_inplace(temp_rot, mat_Ut);
                evaluator.rescale_to_next_inplace(temp_rot);
                temp_rot.scale() = scale;
                // evaluator.mod_switch_to_inplace(temp_rot, last_parms_id);
                evaluator.add(ct_ls[i-1], temp_rot, ct_ls[i]);
                i++;
            }
            result = ct_ls[i-1];
            break;
        }
        case 2:{
            Ciphertext temp_rot1, temp_rot2;
            Plaintext mat_Vk, mat_Vkd;
            mat_Vk = genVk(d, k_m, scale, encoder);
            mat_Vkd = genVkd(d, k_m, scale, encoder);
            evaluator.mod_switch_to_inplace(mat_Vk, parms_id);
            evaluator.mod_switch_to_inplace(mat_Vkd, parms_id);

            evaluator.rotate_vector(ct, k_m, gal_keys, temp_rot1);
            evaluator.rotate_vector(ct, k_m-d, gal_keys, temp_rot2);
            // cout << "a" << endl;
            evaluator.multiply_plain_inplace(temp_rot1, mat_Vk);
            // cout << "b" << endl;
            evaluator.rescale_to_next_inplace(temp_rot1);
            temp_rot1.scale() = scale;
            // evaluator.mod_switch_to_inplace(temp_rot1, last_parms_id);
            // cout << "c" << endl;
            evaluator.multiply_plain_inplace(temp_rot2, mat_Vkd);
            // cout << "d" << endl;
            evaluator.rescale_to_next_inplace(temp_rot2);
            temp_rot2.scale() = scale;
            // evaluator.mod_switch_to_inplace(temp_rot1, last_parms_id);
            // cout << "e" << endl;
            evaluator.add(temp_rot1, temp_rot2, result);
            // cout << "f" << endl;
            break;
        }
        case 3:{
            Plaintext mat_ones;
            mat_ones = genOnes(d, scale, encoder);
            evaluator.mod_switch_to_inplace(mat_ones, parms_id);
            evaluator.rotate_vector(ct, d*k_m, gal_keys, result);
            evaluator.multiply_plain_inplace(result, mat_ones);
            evaluator.rescale_to_next_inplace(result);
            result.scale() = scale;
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
// template <typename T>
// inline Ciphertext mat_mult(Ciphertext ct_a, Ciphertext ct_b, int d, double scale, 
//     CKKSEncoder &encoder, Evaluator &evaluator, const GaloisKeys &gal_keys, 
//     T relin_keys)
// {
//     Ciphertext result, ct_a0, ct_b0;

//     ct_a0 = linear_tran(ct_a, d, 0, 0, scale, encoder, evaluator, gal_keys);
//     ct_b0 = linear_tran(ct_b, d, 1, 0, scale, encoder, evaluator, gal_keys);

//     int i = 0;
//     Ciphertext ab_ls[d];
//     evaluator.multiply(ct_a0, ct_b0, ab_ls[i]);
//     evaluator.relinearize_inplace(ab_ls[i], relin_keys);
//     evaluator.rescale_to_next_inplace(ab_ls[i]);
//     i++;

//     for (int k = 1; k < d; k++){
//         Ciphertext ct_ak, ct_bk, temp_mult;
//         ct_ak = linear_tran(ct_a0, d, 2, k, scale, encoder, evaluator, gal_keys);
//         ct_bk = linear_tran(ct_b0, d, 3, k, scale, encoder, evaluator, gal_keys);
//         evaluator.multiply(ct_ak, ct_bk, temp_mult);
//         evaluator.relinearize_inplace(temp_mult, relin_keys);
//         evaluator.rescale_to_next_inplace(temp_mult);
//         evaluator.add(ab_ls[i-1], temp_mult, ab_ls[i]);
//         i++;
//     }
//     result = ab_ls[i-1];
//     return result;
// }

// inline vector<double> rmat_mult(vector<double> ct_a, vector<double> ct_b, int d, int l)
// {
//     vector<double> result;
//     vector<double> ct_a0;
//     vector<double> ct_b0;

//     ct_a0 = linear_tran(ct_a, d, 0, 0);
//     ct_b0 = linear_tran(ct_b, d, 1, 0);

//     int i = 0;
//     vector<double> ab_ls [d];
//     ab_ls[i] = mult(ct_a0, ct_b0);
//     i++;

//     for (int k = 1; k < l; k++){
//         vector<double> ct_ak;
//         vector<double> ct_bk;

//         ct_ak = linear_tran(ct_a0, d, 2, k);
//         ct_bk = linear_tran(ct_b0, d, 3, k);

//         ab_ls[i] = add(ab_ls[i-1], mult(ct_ak, ct_bk));
//         i++;
//     }
//     vector<double> ab;
//     ab = ab_ls[i-1];
//     for (int k = 0; k < ceil(log2(d/l)); k++){
//         ab_ls[i] = add(ab_ls[i-1], rot(ab, l*d*pow(2, k)));
//         i++;
//     }
//     result = ab_ls[i-1];
//     return result;
// }

//********************************************************************************
// Mackey glass helper functions
//********************************************************************************
#define HARD_MACKEY_GLASS  	30
#define SOFT_MACKEY_GLASS	17

// Define difficulty of the problem
#define MACKEY_GLASS_DIFFICULTY			SOFT_MACKEY_GLASS

// Mackey glass equation
double mackeyglass_eq(double x_t, double x_t_minus_tau, double a, double b) {
	double x_dot = -b*x_t + a*x_t_minus_tau/(1 + pow(x_t_minus_tau,10));
	return x_dot;
}

double mackeyglass_rk4(double x_t, double x_t_minus_tau, double deltat, double a, double b) {
	double k1 = deltat*mackeyglass_eq(x_t,          x_t_minus_tau, a, b);
	double k2 = deltat*mackeyglass_eq(x_t+0.5*k1,   x_t_minus_tau, a, b);
	double k3 = deltat*mackeyglass_eq(x_t+0.5*k2,   x_t_minus_tau, a, b);
	double k4 = deltat*mackeyglass_eq(x_t+k3,       x_t_minus_tau, a, b);
	double x_t_plus_deltat = (x_t + k1/6 + k2/3 + k3/3 + k4/6);
	return x_t_plus_deltat;
}

// Generate a Mackey-Glass time series
void mackey(double *X, double *T, int sample_n) {
	double a        = 0.2;     // value for a in eq (1)
	double b        = 0.1;     // value for b in eq (1)
	int tau      	= MACKEY_GLASS_DIFFICULTY;		// delay constant in eq (1)
	double x0       = 1.2;		// initial condition: x(t=0)=x0
	double deltat   = 0.1;	    // time step size (which coincides with the integration step)
	int interval 	= 1;	    // output is printed at every 'interval' time steps

	double time = 0;
	int index = 1;
	int history_length = floor(tau/deltat);
	double x_history[history_length];
	for (int i = 0; i < x_history[i]; ++i) x_history[i] = 0.0;
	double x_t = x0;
	double x_t_minus_tau, x_t_plus_deltat;

	// Set every value to the default value
	for (int i = 0; i < sample_n; i++) {
		X[i] = x_t;

//		if ((i % interval == 0) && (i > 0)) {
//			printf("%f %f\n", T[i-1], X[i]);
//		}

		if (tau == 0)
			x_t_minus_tau = 0.0;
		else
			x_t_minus_tau = x_history[index];


		x_t_plus_deltat = mackeyglass_rk4(x_t, x_t_minus_tau, deltat, a, b);

		if (tau != 0) {
			x_history[index] = x_t_plus_deltat;
			index = (index % history_length)+1;
		}

		time = time + deltat;
		T[i] = time;
		x_t = x_t_plus_deltat;
	}
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
    size_t poly_modulus_degree = 8192;
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
    GaloisKeys gal_keys = keygen.galois_keys();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << endl;
    cout << endl;

    int p_i = 0;
    parms_id_type parms_ids[3];
    auto context_data = context->first_context_data();
    while (context_data)
    {
        parms_ids[p_i] = context_data->parms_id();
        context_data = context_data->next_context_data();
        p_i++;
    }


    //////////////////////////////////////////////////////////////////////////////
    // Mackey glass settings
    //////////////////////////////////////////////////////////////////////////////
    int sample_all = 10000;	// total no. of samples, excluding the given initial condition
	assert (sample_all >= 2000); // if sample_n < 2000 then the time series is incorrect!!
	double M[sample_all];
	double T[sample_all];
	for (int i = 0; i < sample_all; ++i) M[i] = 0.0;
	for (int i = 0; i < sample_all; ++i) T[i] = 0.0;

    // Generate mackey glass time series data
	mackey(M,T,sample_all);

	// Downsample
	int down_sample = 10;
	int sample_n = sample_all / down_sample;

	// Normalize mackey to -1 1 using hyperbolic tangent
	double X[sample_n];
	for (int i = 0; i < sample_n; i++) {
		X[i] = tanh(M[i*down_sample] - 1.0);
	}
	
    double X_sum;
    double X_avg;
    // Calculate the average value of and print the mackey glass time series data
    // cout << "Mackey glass: " << endl;
	for (int i = 0; i < sample_n; i++) {
        X_sum = X_sum + X[i];
    	// cout << X[i] << ' ';
	}
    X_avg = X_sum/sample_n;
    cout << endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // SEAL Mackey Glass Averaging Example
    //////////////////////////////////////////////////////////////////////////////
    // Initialise the sum variable
    // Plaintext x_sum;
    // encoder.encode(0, scale, x_sum);
    // Ciphertext x_sum_encrypted;
    // encryptor.encrypt(x_sum, x_sum_encrypted);

    // // Loop on all time series values
    // for (size_t i = 0; i < sample_n; i++) {
    //     // Define input plain text
    //     vector<double> input;
    //     input.reserve(slot_count);
    //     input.push_back(X[i]);

    //     // Encode and encrypt input vector
    //     Plaintext x_plain;
    //     encoder.encode(input, scale, x_plain);
    //     Ciphertext x_encrypted;
    //     encryptor.encrypt(x_plain, x_encrypted);

    //     // Sum
    //     evaluator.add_inplace(x_sum_encrypted, x_encrypted);
    // }

    // // Divide to get average
    // Plaintext denom;
    // encoder.encode(0.001, scale, denom);
    // evaluator.multiply_plain_inplace(x_sum_encrypted, denom);

    // // Decrypt, decode, and print the result
    // Plaintext plain_result;
    // decryptor.decrypt(x_sum_encrypted, plain_result);
    // vector<double> result;
    // encoder.decode(plain_result, result);
    // cout << "Result vector: " << endl;
    // print_vector(result, 10, 7);
    // cout << endl;
    // cout << "Correct average value: " << X_avg << endl;
    // cout << endl;

    //////////////////////////////////////////////////////////////////////////////
    // Matrix multiplication
    //////////////////////////////////////////////////////////////////////////////
    /*
    First print the key level parameter information.
    */
    // auto context_data = context->key_context_data();
    // cout << "----> Level (chain index): " << context_data->chain_index();
    // cout << " ...... key_context_data()" << endl;
    // cout << "      parms_id: " << context_data->parms_id() << endl;
    // cout << "      coeff_modulus primes: ";
    // cout << hex;
    // for(const auto &prime : context_data->parms().coeff_modulus())
    // {
    //     cout << prime.value() << " ";
    // }
    // cout << dec << endl;
    // cout << "\\" << endl;
    // cout << " \\-->";

    // /*
    // Next iterate over the remaining (data) levels.
    // */
    // context_data = context->first_context_data();
    // while (context_data)
    // {
    //     cout << " Level (chain index): " << context_data->chain_index();
    //     if (context_data->parms_id() == context->first_parms_id())
    //     {
    //         cout << " ...... first_context_data()" << endl;
    //     }
    //     else if (context_data->parms_id() == context->last_parms_id())
    //     {
    //         cout << " ...... last_context_data()" << endl;
    //     }
    //     else
    //     {
    //         cout << endl;
    //     }
    //     cout << "      parms_id: " << context_data->parms_id() << endl;
    //     cout << "      coeff_modulus primes: ";
    //     cout << hex;
    //     for(const auto &prime : context_data->parms().coeff_modulus())
    //     {
    //         cout << prime.value() << " ";
    //     }
    //     cout << dec << endl;
    //     cout << "\\" << endl;
    //     cout << " \\-->";

    //     /*
    //     Step forward in the chain.
    //     */
    //     context_data = context_data->next_context_data();
    // }
    // cout << " End of chain reached" << endl << endl;
    // /////////////////////////////////////////////////////////////////////

    vector<double> input_a;
    for (int i = 0; i < 4096; i++) {
        input_a.push_back(i);
    }
    vector<double> input_b;
    for (int i = 0; i < 4096; i++) {
        input_b.push_back(4095-i);
    }

    Plaintext a_plain, b_plain;
    encoder.encode(input_a, scale, a_plain);
    encoder.encode(input_b, scale, b_plain);
    Ciphertext a_encrypted, b_encrypted;
    encryptor.encrypt(a_plain, a_encrypted);
    encryptor.encrypt(b_plain, b_encrypted);
    

    Plaintext plain_mult_result;
    Ciphertext cipher_mult_result;

    // //////////////
    // Ciphertext ct_a0;
    // ct_a0 = linear_tran(a_encrypted, 64, 0, 0, scale, encoder, evaluator, gal_keys);
    // // cout << "ct_a0: " << ct_a0.parms_id() << endl;
    // // cout << "a_encrypted: " << a_encrypted.parms_id() << endl;
    // Ciphertext temp_rot1;
    // evaluator.rotate_vector(ct_a0, 1, gal_keys, temp_rot1);
    // Plaintext vk;
    // vk = genVk(64, 1, scale, encoder);
    // // cout << "vk: " << vk.parms_id() << endl;
    // // cout << "temp_rot1: " << temp_rot1.parms_id() << endl;
    // evaluator.multiply_plain_inplace(temp_rot1, vk);
    // // cout << "temp_rot1: " << temp_rot1.parms_id() << endl;
    // cipher_mult_result = temp_rot1;
    // //////////////////

    /////////////////
    int d = 64;
    Ciphertext ct_a0, ct_b0;

    ct_a0 = linear_tran(a_encrypted, d, 0, 0, scale, encoder, evaluator, gal_keys, parms_ids[0]);
    ct_b0 = linear_tran(b_encrypted, d, 1, 0, scale, encoder, evaluator, gal_keys, parms_ids[0]);

    int i = 0;
    Ciphertext ab_ls[d];
    cout << "ct_a0: " << ct_a0.parms_id() << endl;
    cout << "ct_b0: " << ct_b0.parms_id() << endl;
    evaluator.multiply(ct_a0, ct_b0, ab_ls[i]);
    evaluator.relinearize_inplace(ab_ls[i], relin_keys);
    evaluator.rescale_to_next_inplace(ab_ls[i]);
    ab_ls[i].scale() = scale;
    evaluator.mod_switch_to_inplace(ab_ls[i], parms_ids[2]);
    i++;

    for (int k = 1; k < d; k++){
        cout << i << endl;
        Ciphertext ct_ak, ct_bk, temp_mult;
        
        ct_ak = linear_tran(ct_a0, d, 2, k, scale, encoder, evaluator, gal_keys, parms_ids[1]);
        ct_bk = linear_tran(ct_b0, d, 3, k, scale, encoder, evaluator, gal_keys, parms_ids[1]);
        cout << "ct_ak: " << ct_ak.parms_id() << endl;
        cout << "ct_bk: " << ct_bk.parms_id() << endl;
        evaluator.multiply(ct_ak, ct_bk, temp_mult);
        evaluator.relinearize_inplace(temp_mult, relin_keys);
        evaluator.rescale_to_next_inplace(temp_mult);
        temp_mult.scale() = scale;
        evaluator.add(ab_ls[i-1], temp_mult, ab_ls[i]);
        i++;
    }

    cipher_mult_result = ab_ls[i-1];
    //////////////

    decryptor.decrypt(cipher_mult_result, plain_mult_result);
    vector<double> mult_result;
    encoder.decode(plain_mult_result, mult_result);
    cout << "Matrix multiplication result: " << endl;
    print_vector(mult_result, 4096, 3);



    return 0;
}
