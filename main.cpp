#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include "spline.h"
#include <cstring>

long double compute_next_psi(long double current_psi, long double prev_psi, long double current_G, long double prev_G, long double next_G, long double s){

	long double next_psi= 2.0*current_psi - prev_psi;
	next_psi += (5.0 * current_G * current_psi * s * s / 6.0);
	next_psi += (prev_G * prev_psi * s * s /12.0);
	next_psi /= (1- (next_G * s * s / 12.0));

	return next_psi;

}

long double compute_reduced_G(long double x_reduced, long double red_E_guess, int pot_choice, spline func, long double A_param, long double B_param){

	if (pot_choice==3){
		long double x=x_reduced*B_param;
		long double potential=func.splint(x)+func.addition_adjustment;
		potential/=A_param;
		return (2*potential-2*red_E_guess);

	}

}

long double compute_G(long double x, long double E_guess, int pot_choice, spline func, long double mass){

	if (pot_choice==1){
		long double potential = 0.5 * x * x;
		return (2*potential - 2*E_guess)*mass;
	}

	else if (pot_choice==2){

		long double potential= 0.25 * x * x * x * x;
		return (2*potential-2*E_guess)*mass;
	}

	else if (pot_choice==3){ //custom potential, computed using the spline interpolation function adapted from
				//Numerical Recipes in C
		long double potential=(func.splint(x)+func.addition_adjustment);
		return (2*potential-2*E_guess)*mass;
	}

}

long double find_custom_classical_lower(long double E_guess, spline func, long double step_size, long double lower_bound, long double upper_bound){

	long double y_at_x=0.0;

	long double starting_x=1.8;

	while (y_at_x<E_guess && starting_x>lower_bound){
		starting_x-=step_size;
		y_at_x=func.splint(starting_x)+func.addition_adjustment;		

	}

	return starting_x;
}

long double find_custom_classical_upper(long double E_guess, spline func, long double step_size, long double lower_bound, long double upper_bound){

	long double y_at_x=0.0;

        long double starting_x=1.8;

        while (y_at_x<E_guess && starting_x<upper_bound){
                starting_x+=step_size;
                y_at_x=func.splint(starting_x)+func.addition_adjustment;
        }

        return starting_x;

}

long double find_classical_limit(long double E_guess, int pot_choice){

	if (pot_choice==1) return (std::sqrt(2*E_guess));

	else if (pot_choice == 2) return (std::sqrt(std::sqrt(4*E_guess)));

}


int main(int argc, char * argv[]){

	if (argc<5){
		std::cerr<<"Usage: ./Numerov-Analysis num_energy_levels system_mass x_min x_max < potential_input.dat \n";
		return 1;
	}

	std::stringstream ss(argv[1]);
	int num_energy_levels=0;
	ss>>num_energy_levels;

	int potential_choice=3;

	long double mass=0.0;
	std::stringstream ii(argv[2]);
	ii>>mass;

	long double x_min, x_max;
	std::stringstream si(argv[3]);
	si>>x_min;
	std::stringstream sis(argv[4]);	
	sis>>x_max;

	spline custom_function;

	if (potential_choice==3) custom_function.initialize(std::cin);

	long double hbar=1.0;

	long double A_param= hbar/std::sqrt(mass);
	long double B_param= std::sqrt(hbar)/(std::sqrt(std::sqrt(mass)));

	long double s;
	long double red_x_max, red_x_min;

	s=0.0001;

	custom_function.adjust_zero(s, x_min, x_max);

	long double psi_length;

	long double * wavefunction;
	long double * domain;

	int num_intervals=0;

	long double wave_energy=0;

	long double reduced_s=s/B_param;
	long double * reduced_wavefunction;
	long double * reduced_domain;

	long double red_E_guess=0.1L;

	long double increment=0.1L;

	int energy_level_being_computed=0;

	std::ofstream output_stream;

	while (energy_level_being_computed<num_energy_levels){

		bool wavefunction_found=false;
		bool guess_high=false;

		increment=0.1L;

		red_E_guess+=increment;

		while (wavefunction_found == false){
			int node_counter=0;

			num_intervals= (x_max-x_min)/s;

			red_x_max=x_max/B_param;
			red_x_min=x_min/B_param;

			reduced_domain=new long double[num_intervals];
			reduced_wavefunction=new long double[num_intervals];

			reduced_wavefunction[0]=0;
			reduced_wavefunction[1]=0.00001;

			for (int i=0; i<num_intervals; ++i){

				reduced_domain[i]=red_x_min + (i*reduced_s);
			
				if (i>1){

					long double prev_G_red=compute_reduced_G(reduced_domain[i-1], red_E_guess, potential_choice, custom_function, A_param, B_param);
                	                long double current_G_red=compute_reduced_G(reduced_domain[i], red_E_guess, potential_choice, custom_function, A_param, B_param);
                        	        long double next_G_red=compute_reduced_G((red_x_min+((i+1)*reduced_s)), red_E_guess, potential_choice, custom_function, A_param, B_param);

                                	reduced_wavefunction[i]=compute_next_psi(reduced_wavefunction[i-1], reduced_wavefunction[i-2], current_G_red, prev_G_red, next_G_red, reduced_s);
                        
				}

			}

			//wavefunction must be normalized before boundary conditions are checked

			long double reduced_total_sum=0.0;
	  		for (int i=0; i<num_intervals; ++i){
        	       		reduced_total_sum+=(reduced_wavefunction[i]*reduced_wavefunction[i]*reduced_s);
	        	}

        		long double reduced_normalization_factor=1.0/std::sqrt(reduced_total_sum);
		        for (int i=0; i<num_intervals; ++i){
        	        	reduced_wavefunction[i]*=reduced_normalization_factor;
        		}


			for (int i=0; i<(num_intervals-100); ++i){
				if (reduced_wavefunction[i]*reduced_wavefunction[i+1]<0) ++node_counter;
			}

			if ((reduced_wavefunction[num_intervals-1]*reduced_wavefunction[num_intervals-1]) < (0.005) && node_counter==energy_level_being_computed){

        	                wavefunction_found=true;
                	        wave_energy=red_E_guess;

	                }

			if (node_counter>energy_level_being_computed){
			
				red_E_guess-=increment;
				guess_high=true;
				increment/=10.0L;
			}

			else if (node_counter<=energy_level_being_computed){
	
				red_E_guess+=increment;
				guess_high=false;

			}

			if (wavefunction_found==false){

				delete reduced_domain;
				delete reduced_wavefunction;

			}

		}
		std::cerr<<"Energy Level: "<<energy_level_being_computed<<'\n';

		std::cerr<<"Reduced-Energy eigenvalue (in a.u.): "<<std::setprecision(10)<<(wave_energy-custom_function.addition_adjustment)*A_param<<'\n';

		std::cerr<<"Dimensional Energy eigenvalue (in kcal/mol): "<<std::setprecision(10)<<(((wave_energy-custom_function.addition_adjustment)*A_param)*627.503)<<"\n\n";

		int zero=0;

		char str[10];
		strcpy(str, "E_level");
                std::stringstream sss;
                sss<<energy_level_being_computed;
                std::string hold;
                sss>>hold;
		std::stringstream ssh;
		ssh<<mass;
		hold+="_m";
		std::string hold1;
		ssh>>hold1;
		hold+=hold1;

                strcat(str, hold.c_str());

                output_stream.open(str);

		for (int i=0; i<num_intervals; ++i){

			output_stream
				<<std::setw(18)<<(reduced_domain[i]*B_param)
				<<std::setw(18)<<(reduced_wavefunction[i]/std::sqrt(B_param))
				<<std::setw(18)<<(reduced_wavefunction[i]/std::sqrt(B_param))*(reduced_wavefunction[i]/std::sqrt(B_param))
				<<'\n';
		}
		output_stream.close();


		delete reduced_domain;
		delete reduced_wavefunction;

		++energy_level_being_computed;

	}

	return 0;

}
