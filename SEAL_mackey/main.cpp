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
// Mackey glass helper functions
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
    //****************************************************************************

    //****************************************************************************
    // Mackey glass settings
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
    cout << "Mackey glass data: " << endl;
	for (int i = 0; i < sample_n; i++) {
        X_sum = X_sum + X[i];
    	cout << X[i] << ' ';
	}
    X_avg = X_sum/sample_n;
    cout << endl;

    //****************************************************************************
    // SEAL
    //****************************************************************************
    
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
    print_vector(result, 5, 7);
    cout << endl;
    cout << "Correct average value: " << X_avg << endl;

    return 0;
}
