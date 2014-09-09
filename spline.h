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

	void initial_spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	float splint(float x);
	void initialize(std::istream& ifs);

	void print_spline(float x_step, float lower_bound, float upper_bound);

	void adjust_zero(float x_step, float lower_bound, float upper_bound);

	float *x_array;
	float *y_array;

	int num_data;

	float yp1, ypn;

	float addition_adjustment;

	float *y2_array;
};

#endif
