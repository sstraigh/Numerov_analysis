#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <istream>

class wavefunction{

public:

	void initialize_wavefunction(const long double x_start, const long double x_end, const long double 
		     interval, const long double input_mass);

	long double get_energy();
	long double get_psi_at_index(int index);

	long double get_node_num(){return num_nodes;}
	long double get_reduced_x_min(){return reduced_x_min;}
	long double get_reduced_x_step(){return reduced_x_step;}

	long double get_A_param(){return A_param;}
	long double get_B_param(){return B_param;}

	void compute_next_psi(long double current_G_reduced, long double 
			      prev_G_reduced, long double next_G_reduced,
			      int interval_index);

	bool check_wavefunction(int energy_level);
	void deinitialize();

	void set_energy(long double red_E_guess, long double add_adjustment);
	long double energy = 0.0;
	
	void free_holder(){delete[] wavefunction_holder;}

private:
	void normalize_wavefunction();

	long double x_step, reduced_x_step;
	long double reduced_x_min, reduced_x_max;
	long double x_min, x_max;
	long double A_param;
	long double B_param;
	long double system_mass;
	long double * wavefunction_holder;

	int num_nodes=0;
	int num_intervals;

};

#endif
