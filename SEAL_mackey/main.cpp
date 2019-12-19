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
    cout << endl;
    //****************************************************************************

    // //****************************************************************************
    // // Mackey glass settings
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
    cout << "Mackey glass: " << endl;
	for (int i = 0; i < sample_n; i++) {
        X_sum = X_sum + X[i];
    	// cout << X[i] << ' ';
	}
    X_avg = X_sum/sample_n;
    cout << endl;

    //****************************************************************************
    // SEAL Examples
    //****************************************************************************
    
    // //****************************************************************************
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
    // //****************************************************************************

    //****************************************************************************
    // Approximation of sigmoid function
    //****************************************************************************
    vector<double> input;
    input.reserve(slot_count);
    // Set the input data
    double input_data = 0.35;
    input.push_back(input_data);
    cout << "Approximation of sigmoid function: " << endl;
    cout << "Input data: " << endl;
    print_vector(input, 10, 10);

    cout << "Evaluating polynomial f(x) = 0.5+0.197x-0.004x^3 ..." << endl;

    vector<double> coeff1;
    coeff1.reserve(slot_count);
    coeff1.push_back(0.5);
    Plaintext coeff1_plain;
    encoder.encode(coeff1, scale, coeff1_plain);
    Ciphertext coeff1_encrypted;
    encryptor.encrypt(coeff1_plain, coeff1_encrypted);

    Plaintext plain_coeff2, plain_coeff3;
    encoder.encode(0.197, scale, plain_coeff2);
    encoder.encode(-0.004, scale, plain_coeff3);

    Plaintext x_plain;
    encoder.encode(input, scale, x_plain);
    Ciphertext x1_encrypted;
    encryptor.encrypt(x_plain, x1_encrypted);

    /*
    To compute x^3 we first compute x^2 and relinearize. However, the scale has
    now grown to 2^80.
    */
    print_line(__LINE__);
    cout << "Compute x^2 and relinearize:" << endl;
    Ciphertext x3_encrypted;
    evaluator.square(x1_encrypted, x3_encrypted);
    evaluator.relinearize_inplace(x3_encrypted, relin_keys);
    cout << "    + Scale of x^2 before rescale: " << log2(x3_encrypted.scale())
        << " bits" << endl;
    /*
    Now rescale; in addition to a modulus switch, the scale is reduced down by
    a factor equal to the prime that was switched away (40-bit prime). Hence, the
    new scale should be close to 2^40. Note, however, that the scale is not equal
    to 2^40: this is because the 40-bit prime is only close to 2^40.
    */
    print_line(__LINE__);
    cout << "Rescale x^2." << endl;
    evaluator.rescale_to_next_inplace(x3_encrypted);
    cout << "    + Scale of x^2 after rescale: " << log2(x3_encrypted.scale())
        << " bits" << endl;

    /*
    Now x3_encrypted is at a different level than x1_encrypted, which prevents us
    from multiplying them to compute x^3. We could simply switch x1_encrypted to
    the next parameters in the modulus switching chain. However, since we still
    need to multiply the x^3 term with -0.004 (plain_coeff3), we instead compute -0.004*x
    first and multiply that with x^2 to obtain -0.004*x^3. To this end, we compute
    -0.004*x and rescale it back from scale 2^80 to something close to 2^40.
    */
    print_line(__LINE__);
    cout << "Compute and rescale -0.004*x." << endl;
    Ciphertext x1_encrypted_coeff3;
    evaluator.multiply_plain(x1_encrypted, plain_coeff3, x1_encrypted_coeff3);
    cout << "    + Scale of -0.004*x before rescale: " << log2(x1_encrypted_coeff3.scale())
        << " bits" << endl;
    evaluator.rescale_to_next_inplace(x1_encrypted_coeff3);
    cout << "    + Scale of -0.004*x after rescale: " << log2(x1_encrypted_coeff3.scale())
        << " bits" << endl;


    /*
    Since x3_encrypted and x1_encrypted_coeff3 have the same exact scale and use
    the same encryption parameters, we can multiply them together. We write the
    result to x3_encrypted, relinearize, and rescale. Note that again the scale
    is something close to 2^40, but not exactly 2^40 due to yet another scaling
    by a prime. We are down to the last level in the modulus switching chain.
    */
    print_line(__LINE__);
    cout << "Compute, relinearize, and rescale (-0.004*x)*x^2." << endl;
    evaluator.multiply_inplace(x3_encrypted, x1_encrypted_coeff3);
    evaluator.relinearize_inplace(x3_encrypted, relin_keys);
    cout << "    + Scale of -0.004*x^3 before rescale: " << log2(x3_encrypted.scale())
        << " bits" << endl;
    evaluator.rescale_to_next_inplace(x3_encrypted);
    cout << "    + Scale of -0.004*x^3 after rescale: " << log2(x3_encrypted.scale())
        << " bits" << endl;

    /*
    Next we compute the degree one term. All this requires is one multiply_plain
    with plain_coeff1. We overwrite x1_encrypted with the result.
    */
    print_line(__LINE__);
    cout << "Compute and rescale 0.197*x." << endl;
    evaluator.multiply_plain_inplace(x1_encrypted, plain_coeff2);
    cout << "    + Scale of 0.197*x before rescale: " << log2(x1_encrypted.scale())
        << " bits" << endl;
    evaluator.rescale_to_next_inplace(x1_encrypted);
    cout << "    + Scale of 0.197*x after rescale: " << log2(x1_encrypted.scale())
        << " bits" << endl;

    /*
    Now we would hope to compute the sum of all three terms. However, there is
    a serious problem: the encryption parameters used by all three terms are
    different due to modulus switching from rescaling.

    Encrypted addition and subtraction require that the scales of the inputs are
    the same, and also that the encryption parameters (parms_id) match. If there
    is a mismatch, Evaluator will throw an exception.
    */
    cout << endl;
    print_line(__LINE__);
    cout << "Parameters used by all three terms are different." << endl;
    cout << "    + Modulus chain index for x3_encrypted: "
        << context->get_context_data(x3_encrypted.parms_id())->chain_index() << endl;
    cout << "    + Modulus chain index for x1_encrypted: "
        << context->get_context_data(x1_encrypted.parms_id())->chain_index() << endl;
    cout << "    + Modulus chain index for coeff1_encrypted: "
        << context->get_context_data(coeff1_encrypted.parms_id())->chain_index() << endl;
    cout << endl;

    /*
    Let us carefully consider what the scales are at this point. We denote the
    primes in coeff_modulus as P_0, P_1, P_2, P_3, in this order. P_3 is used as
    the special modulus and is not involved in rescalings. After the computations
    above the scales in ciphertexts are:

        - Product x^2 has scale 2^80 and is at level 2;
        - Product -0.004*x has scale 2^80 and is at level 2;
        - We rescaled both down to scale 2^80/P_2 and level 1;
        - Product -0.004*x^3 has scale (2^80/P_2)^2;
        - We rescaled it down to scale (2^80/P_2)^2/P_1 and level 0;
        - Product 0.197*x has scale 2^80;
        - We rescaled it down to scale 2^80/P_2 and level 1;
        - The constant term 0.5 has scale 2^40 and is at level 2.

    Although the scales of all three terms are approximately 2^40, their exact
    values are different, hence they cannot be added together.
    */
    print_line(__LINE__);
    cout << "The exact scales of all three terms are different:" << endl;
    ios old_fmt(nullptr);
    old_fmt.copyfmt(cout);
    cout << fixed << setprecision(10);
    cout << "    + Exact scale in -0.004*x^3: " << x3_encrypted.scale() << endl;
    cout << "    + Exact scale in  0.197*x: " << x1_encrypted.scale() << endl;
    cout << "    + Exact scale in      0.5: " << coeff1_encrypted.scale() << endl;
    cout << endl;
    cout.copyfmt(old_fmt);

    /*
    There are many ways to fix this problem. Since P_2 and P_1 are really close
    to 2^40, we can simply "lie" to Microsoft SEAL and set the scales to be the
    same. For example, changing the scale of -0.004*x^3 to 2^40 simply means that we
    scale the value of -0.004*x^3 by 2^120/(P_2^2*P_1), which is very close to 1.
    This should not result in any noticeable error.

    Another option would be to encode 1 with scale 2^80/P_2, do a multiply_plain
    with 0.197*x, and finally rescale. In this case we would need to additionally
    make sure to encode 1 with appropriate encryption parameters (parms_id).

    In this example we will use the first (simplest) approach and simply change
    the scale of -0.004*x^3 and 0.197*x to 2^40.
    */
    print_line(__LINE__);
    cout << "Normalize scales to 2^40." << endl;
    x3_encrypted.scale() = pow(2.0, 40);
    x1_encrypted.scale() = pow(2.0, 40);

    /*
    We still have a problem with mismatching encryption parameters. This is easy
    to fix by using traditional modulus switching (no rescaling). CKKS supports
    modulus switching just like the BFV scheme, allowing us to switch away parts
    of the coefficient modulus when it is simply not needed.
    */
    print_line(__LINE__);
    cout << "Normalize encryption parameters to the lowest level." << endl;
    parms_id_type last_parms_id = x3_encrypted.parms_id();
    evaluator.mod_switch_to_inplace(x1_encrypted, last_parms_id);
    evaluator.mod_switch_to_inplace(coeff1_encrypted, last_parms_id);

    /*
    All three ciphertexts are now compatible and can be added.
    */
    print_line(__LINE__);
    cout << "Compute -0.004*x^3 + 0.197*x + 0.5." << endl;
    Ciphertext encrypted_result;
    evaluator.add(x3_encrypted, x1_encrypted, encrypted_result);
    evaluator.add_inplace(encrypted_result, coeff1_encrypted);

    /*
    Decrypt, decode, and print the result.
    */
    print_line(__LINE__);
    cout << "Decrypt and decode -0.004*x^3 + 0.197x + 0.5." << endl;
    Plaintext plain_result_sig;
    decryptor.decrypt(encrypted_result, plain_result_sig);
    vector<double> sig_result;
    encoder.decode(plain_result_sig, sig_result);
    cout << "Computed result: " << endl;
    print_vector(sig_result, 10, 10);

    // Print true result
    vector<double> true_result;
    true_result.push_back(0.5+0.197*input_data-0.004*input_data*input_data*input_data);
    cout << "True result: " << endl;
    print_vector(true_result, 10, 10);
    //****************************************************************************

    return 0;
}
