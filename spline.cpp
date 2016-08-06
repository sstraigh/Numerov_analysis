//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 8 Sep 2014

//implementation of the class "spline", which reads in an arbitrary potential
//specified by two columns from STDIN, and provides interpolation of the
//function within the specified domain

#include "spline.h"

//constructor
spline::spline(const long double *x_values, const long double *y_values, int num_vals){

        num_data=num_vals;

  //      x_array= new long double[num_data];
//        y_array= new long double[num_data];

        for (int i=1; i<=num_data; ++i){
            y_array[i] = y_values[i - 1];
            x_array[i] = x_values[i - 1];
        }

	long double yp1 = (y_array[1] - y_array[2]) / (x_array[1] - x_array[2]);
	
	long double ypn = (y_array[num_data - 1] - y_array[num_data - 0]) / (x_array[num_data - 1] - x_array[num_data - 0]);

   //     y2_array=new long double[num_data];
        std::fill(y2_array, y2_array+(1+num_data), 0.0);

//	std::cout<<y_array[7]<<std::endl;
//	std::cout<<y_array[7]<<std::endl;
        //arrays of doubles are now filled with relevant (and arbitrary)
        //y and x values

        initial_spline((num_data-0), yp1, ypn);

//	std::cout<<y2_array[7]<<std::endl;
//
//	std::cout<<splint(1.0)<<std::endl;

        //spline is initialized with the x and y data read in and interpolated
}


//print_spline is a routine to double check the spline function interpolates
//data correctly within upper and lower bound. x_step is the interval between
//successive domain values for which the function is evaluated
void spline::print_spline(long double x_step, long double lower_bound, long double upper_bound){

	long double x_pos=lower_bound;

	while (x_pos<=upper_bound){

		std::cerr
			<<std::setw(18)<<x_pos
			<<std::setw(18)<<splint(x_pos)
			<<'\n';

		x_pos+=x_step;

	}

}

//adjust_zero sets the lowest point in the interpolated function to zero-if
//the function is positive, the "addition_adjustment" is a negative value 
//and vice-versa
void spline::adjust_zero(long double x_step, long double lower_bound, long double upper_bound){

	long double x_pos=lower_bound;

	long double lowest_y_val=0.0;

		//loop to find the lowest y value in the specified domain

	while (x_pos<=upper_bound){

		if (lowest_y_val>splint(x_pos)) lowest_y_val=splint(x_pos);

		x_pos+=x_step;

	}

	addition_adjustment= -1.0* lowest_y_val;


}

//creates the y-vector which, for a particular function, holds the values of the function.
//Adapted from numerical recipes in C
void spline::initial_spline(int n, long double yp1, long double ypn){

	int i,k;
	long double p, qn, sig, un;

        long double * u = new long double[n];

	if (yp1 > 0.99e30) y2_array[1]=u[1]=0.0;

	else{
		y2_array[1]=-0.5;
	
		u[1]=(3.0/(x_array[2]-x_array[1]))*((y_array[2]-y_array[1])/(x_array[2]-x_array[1])-yp1);
	}

//	std::cout<<y2_array[1]<<std::endl;
//	std::cout<<x_array[n-1]<<'\t'<<y_array[n-1]<<std::endl;
//	std::cout<<x_array[n]<<'\t'<<y_array[n]<<std::endl;
//	std::cout<<x_array[n+1]<<'\t'<<y_array[n+1]<<std::endl;	

	for (i = 2; i <= n-1; ++i){

		sig=(x_array[i]-x_array[i-1])/(x_array[i+1]-x_array[i-1]);
		p=sig*y2_array[i-1]+2.0;
		y2_array[i]=(sig-1.0)/p;

		u[i]=(y_array[i+1]-y_array[i])/(x_array[i+1]-x_array[i]) - (y_array[i]-y_array[i-1])/(x_array[i]-x_array[i-1]);
		u[i]=(6.0*u[i]/(x_array[i+1]-x_array[i-1])-sig*u[i-1])/p;

//		std::cout<<x_array[i]<<'\t'<<y_array[i]<<'\t'<<y2_array[i]
//			 <<'\t'<<u[i]<<std::endl;

	}

	if(ypn>0.99e30) qn=un=0.0;

	else{
		qn=0.5;
		un=(3.0/(x_array[n]-x_array[n-1]))*(ypn-(y_array[n]-y_array[n-1])/(x_array[n]-x_array[n-1]));
	}

//	std::cout<<y2_array[n]<<'\t'<<un<<'\t'<<qn<<'\t'<<u[n-1]
//		 <<'\t'<<y2_array[n-1]<<std::endl;
	y2_array[n]=(un-qn*u[n-1])/(qn*y2_array[n-1]+1.0);
//	std::cout<<y2_array[n]<<std::endl;

	for (k=n-1; k>=1; --k){
	    y2_array[k]=y2_array[k]*y2_array[k+1]+u[k];
//	    std::cout<<y2_array[k]<<std::endl;
	}
//	std::cout<<y2_array[n-1]<<std::endl;

	delete[] u;

}
//splint (SPLine INTerpolation) finds the value of the interpolated function
//for a given value within the function's domain.
//From Numerical Recipes in C
long double spline::splint(long double x){

//	void nrerror(char error_text[]);

	int klo, khi, k;

	long double h,b,a;

	klo=1;
	khi=num_data;
	while (khi-klo>1){
		k=(khi+klo) >> 1;
		if (x_array[k] > x ) khi=k;
		else klo=k;
	}

	h=x_array[khi]-x_array[klo];

	if (h==0.0){
		std::cout<<"Bad xa input to routine splint\n";
		exit(1);
	}

	a=(x_array[khi]-x)/h;
	b=(x-x_array[klo])/h;
	return (a*y_array[klo]+b*y_array[khi]+((a*a*a-a)*y2_array[klo]+(b*b*b-b)*y2_array[khi])*(h*h)/6.0);
}

