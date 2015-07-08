#ifndef SPLINE_H
#define SPLINE_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <istream>

class spline{

public:
	spline(const long double *x_values, const long double *y_values, int num_vals);
	long double splint(long double x);
//	void initialize(long double *x_values, long double *y_values, int num_vals); 
	long double addition_adjustment = 0.0;
	void adjust_zero(long double x_step, long double lower_bound, long double upper_bound);
	void print_spline(long double x_step, long double lower_bound, long double upper_bound);
/*
	~spline(){
		delete[] x_array;
		delete[] y_array;
		delete[] y2_array;
	}

*/
private:
	void initial_spline(int n, long double yp1, long double ypn);

	long double x_array[14];
	long double y_array[14];

	int num_data;

	long double y2_array[14];
};

#endif
