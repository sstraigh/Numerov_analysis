//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 8 Sep 2014

//implementation of the class "spline", which reads in an arbitrary potential
//specified by two columns from STDIN, and provides interpolation of the
//function within the specified domain

#include "spline.h"

//print_spline is a routine to double check the spline function interpolates
//data correctly within upper and lower bound. x_step is the interval between
//successive domain values for which the function is evaluated
void spline::print_spline(float x_step, float lower_bound, float upper_bound){

	float x_pos=lower_bound;

	while (x_pos<=upper_bound){

		std::cout
			<<std::setw(18)<<x_pos
			<<std::setw(18)<<(splint(x_pos)+addition_adjustment)
			<<'\n';

		x_pos+=x_step;

	}

}

//adjust_zero sets the lowest point in the interpolated function to zero-if
//the function is positive, the "addition_adjustment" is a negative value 
//and vice-versa
void spline::adjust_zero(float x_step, float lower_bound, float upper_bound){

	float x_pos=lower_bound;

	float lowest_y_val=0.0;

		//loop to find the lowest y value in the specified domain

	while (x_pos<=upper_bound){

		if (lowest_y_val>splint(x_pos)) lowest_y_val=splint(x_pos);

		x_pos+=x_step;

	}

	addition_adjustment= -1.0* lowest_y_val;


}

//creates the y-vector which, for a particular function, holds the values of the function.
//Adapted from numerical recipes in C
void spline::initial_spline(float x[], float y[], int n, float yp1, float ypn, float y2[]){

	int i,k;
	float p, qn, sig, un;

	float * u=new float[n];

	if (yp1 > 0.99e30) y2[1]=u[1]=0.0;

	else{
		y2[1]=-0.5;
	
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	
	for (i=2; i<=n-1; ++i){

		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	if(ypn>0.99e30) qn=un=0.0;

	else{
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}

	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1; k>=1; --k) y2[k]=y2[k]*y2[k+1]+u[k];

	delete u;
}
//splint (SPLine INTerpolation) finds the value of the interpolated function
//for a given value within the function's domain.
//From Numerical Recipes in C
float spline::splint(float x){

//	void nrerror(char error_text[]);

	int klo, khi, k;

	float h,b,a;

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
	}

	a=(x_array[khi]-x)/h;
	b=(x-x_array[klo])/h;
	return (a*y_array[klo]+b*y_array[khi]+((a*a*a-a)*y2_array[klo]+(b*b*b-b)*y2_array[khi])*(h*h)/6.0);
}

//Initialization of a spline object-implementation as written here
//is specific for use in Numerov-Analysis.
void spline::initialize(std::istream& ifs){

	yp1=5.0e30;
	ypn=yp1;

	num_data=0;

	std::vector<float> x,y;
	
	while(ifs){

		std::string xholder;
		std::cin>>xholder;
		
		std::stringstream ssx(xholder);
		float xtemp;
		ssx>>xtemp;
		
		x.push_back(xtemp);

		std::string yholder;
		std::cin>>yholder;
			
		std::stringstream ssy(yholder);
		float ytemp;
		ssy>>ytemp;

		y.push_back(ytemp);

		++num_data;
	}

	x_array=new float[num_data];
	std::fill(x_array, x_array+num_data, 0.0);

	y_array=new float[num_data];
	std::fill(y_array, y_array+num_data, 0.0);
	
	int index=1;
	std::vector<float>::iterator yit=y.begin();
	for (std::vector<float>::iterator xit=x.begin(); (xit!=x.end() && index<num_data); ++xit){

		x_array[index]=(*xit);
		y_array[index]=(*yit);

		++yit;

		++index;
	}

	y2_array=new float[num_data];
	std::fill(y2_array, y2_array+num_data, 0.0);	

	initial_spline(x_array, y_array, num_data, yp1, ypn, y2_array);

}
