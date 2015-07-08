//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 17 Fed 2015

#include "wavefunction.h"
#include <cmath>

//initializer
void wavefunction::initialize_wavefunction(const long double x_start, const long double 
					   x_end, const long double interval, const 
					   long double input_mass){

	x_step = interval;
	x_min = x_start;
	x_max = x_end;
	system_mass = input_mass;
	
	num_intervals = (x_max - x_min)/x_step;
	
	long double hbar = 1.0;

	A_param = hbar/std::sqrt(system_mass);
	B_param = std::sqrt(hbar)/(std::sqrt(std::sqrt(system_mass)));
	
	reduced_x_step = x_step / B_param;
	reduced_x_max = x_max / B_param;
	reduced_x_min = x_min / B_param;
	
	wavefunction_holder = new long double[num_intervals];

	std::fill(wavefunction_holder, wavefunction_holder + num_intervals, 0.0);

	wavefunction_holder[0] = 0.0;
	wavefunction_holder[1] = x_step;
}
void wavefunction::set_energy(long double red_E_guess, long double add_adjustment){

	energy = red_E_guess - add_adjustment;
	energy *= A_param;

}

bool wavefunction::check_wavefunction(int energy_level){

	normalize_wavefunction();

	int node_counter = 0;

//	std::cout<<"here 1"<<std::endl;

	//check if there are the correct number of nodes in psi(x)
	for (int i =0; i < num_intervals-1; ++i){
		if (wavefunction_holder[i]*wavefunction_holder[i+1] < 0){
			++node_counter;
		}
	}

//	std::cout<<"here 2"<<std::endl;
	
	num_nodes = node_counter;

	if (node_counter != energy_level) return false;

	long double taper_value = wavefunction_holder[num_intervals-1] * 
			     wavefunction_holder[num_intervals-2];

	if (taper_value < x_step) return true;

//	std::cout<<"here 3"<<std::endl;

	return false;

}
void wavefunction::deinitialize(){

	for (int i=0; i<num_intervals; ++i) wavefunction_holder[i]=0.0;

	wavefunction_holder[1]=x_step;

}


void wavefunction::compute_next_psi(long double current_G_reduced, long double
                              prev_G_reduced, long double next_G_reduced,
			      int interval_index){

	long double next_psi = 2.0 * wavefunction_holder[interval_index] 
			  - wavefunction_holder[interval_index - 1];

	next_psi += 5.0 * current_G_reduced * 
		    wavefunction_holder[interval_index] * reduced_x_step
		    * reduced_x_step / 6.0;
	next_psi += prev_G_reduced * wavefunction_holder[interval_index - 1]
		    * reduced_x_step * reduced_x_step / 12.0;
	next_psi /= (1 - (next_G_reduced * reduced_x_step * reduced_x_step
		    / 12.0));

	wavefunction_holder[interval_index + 1] = next_psi;

//	std::cout<<next_psi<<'\t'<<wavefunction_holder[interval_index + 1]<<std::endl;
}


long double wavefunction::get_psi_at_index(int index){

	return wavefunction_holder[index];

}

long double wavefunction::get_energy(){

	return energy;
}
void wavefunction::normalize_wavefunction(){

	long double reduced_total_sum = 0.0;
	for (int i = 0; i < num_intervals; ++i){
		reduced_total_sum += wavefunction_holder[i]*
	   			     wavefunction_holder[i]*
				     reduced_x_step;
	}

	long double normalization_factor = 1.0 / std::sqrt(reduced_total_sum);

	for (int i = 0; i < num_intervals; ++i){
		wavefunction_holder[i]*=normalization_factor;
	}


}
