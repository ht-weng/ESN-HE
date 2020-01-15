#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;

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

int main() {

	int sample_all = 10000;	// total no. of samples, excluding the given initial condition
	assert (sample_all >= 2000); // if sample_n < 2000 then the time series is incorrect!!
	double M[sample_all];
	double T[sample_all];
	for (int i = 0; i < sample_all; ++i) M[i] = 0.0;
	for (int i = 0; i < sample_all; ++i) T[i] = 0.0;

	mackey(M,T,sample_all);

	// Downsample (very important!)
	int down_sample = 10;
	int sample_n = sample_all / down_sample;

	// Normalize mackey to -1 1 using hyperbolic tangent
	double X[sample_n];
	for (int i = 0; i < sample_n; i++) {
		X[i] = tanh(M[i*down_sample] - 1.0);
	}
	
	cout << "Sample No: " << sample_n << endl;
	cout << "Sample Data: " << endl;
	for (int i = 0; i < sample_n; i++) {
    	cout << X[i] << ' ';
	}
	cout << endl;

	return 0;
}
