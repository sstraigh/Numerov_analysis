#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "wavefunction.h"
#include "spline.h"
#include <cstring>

const long double electron_proton_mass_ratio = 1836.15267245;
const long double hartree_to_wavenumbers = 219474.63;
const long double kcal_to_hartree = 0.00159362;
const long double angstroms_to_bohr = 1.889725989;

void propagate_computation(int num_energy_levels, long double mass, long double x_min,
                           long double x_max, long double total_time, int num_points);

//Computes the next value in the wavefunction using the current value [psi(x)]
//and the previous valu [psi(x-s)]
long double compute_next_psi(long double current_psi, long double prev_psi, long double current_G, long double prev_G, long double next_G, long double s){

    long double next_psi= 2.0*current_psi - prev_psi;
    next_psi += (5.0 * current_G * current_psi * s * s / 6.0);
    next_psi += (prev_G * prev_psi * s * s /12.0);
    next_psi /= (1- (next_G * s * s / 12.0));

    return next_psi;

}

//computes the next G-value (for further reference and an overview of the
//Numerov algorithm, see Levine, Quatum Chemistry, 6th ed)
long double compute_reduced_G(long double x_reduced, long double red_E_guess, 
		long double potential, long double A_param){

    long double red_potential = potential/A_param;
    
    return 2*(red_potential - red_E_guess);
    
}

void read_in_frame(long double * scan_position,
  long double * pot, std::ifstream &potential, int num_points){

    for(int i=0; i<num_points; ++i){

	//readin relevant x-y values from data files needed
	//x-values the same for all files;

    std::string potential_line;

	std::getline(potential, potential_line);

	std::istringstream potential_ss(potential_line);

    double junk;
	//for potential, dealing with a two column file, x vs y
	potential_ss >> junk >> scan_position[i] >> pot[i];
	
	scan_position[i] *= angstroms_to_bohr; //convert from A to Bohr
	pot[i] *= kcal_to_hartree;

    }

}

void read_input_file(int &num_levels, long double &mass, long double &min_x, long double &max_x, long double &total_time){

    /* INPUT.dat organized as in the following example:

            num_levels = 3
            mass       = M (for OD, = 1.78882340923978815770)
            min_x      = 0.6 (minimum in angstroms)
            max_x      = 1.5 (maximum in angstroms)
            total_time = 27.0 (last time-step in picoseconds)
    */

    std::ifstream input_file("INPUT.dat");
    std::string level_line, mass_line, min_x_line, max_x_line, time_line, tmp;
    
    getline(input_file, level_line);
    getline(input_file,  mass_line);
    getline(input_file, min_x_line);
    getline(input_file, max_x_line);
    getline(input_file,  time_line);
    
    std::stringstream num_levels_ss(level_line);
    num_levels_ss >> tmp >> tmp >> num_levels;
    
    std::stringstream mass_stream(mass_line);
    mass_stream >> tmp >> tmp >> mass;
    
    std::stringstream min_x_stream(min_x_line);
    min_x_stream >> tmp >> tmp >> min_x;
    
    std::stringstream max_x_stream(max_x_line);
    max_x_stream >> tmp >> tmp >> max_x;
    
    std::stringstream total_time_stream(time_line);
    total_time_stream >> tmp >> tmp >> total_time;

    mass *= electron_proton_mass_ratio;
    min_x *= angstroms_to_bohr;
    max_x *= angstroms_to_bohr;
}



int main(int argc, char * argv[]){

    if (argc != 2){
	std::cerr<<"Usage: ./Numerov-Analysis num_points (read other input from INPUT.dat)\n";
	return 1;
    }

    //initialize variable holders: wavefunctions, values of wavefunctions

    int num_energy_levels=0;
    long double mass=0.0;
    long double x_min, x_max;
    long double total_time=0.0;
    read_input_file(num_energy_levels, mass, x_min, x_max, total_time);

    std::istringstream iss(argv[1]);
    int num_points;
    iss >> num_points;
    if (iss.fail()){
      std::cerr << "cannot parse sstream: " << iss << std::endl;
    	exit(1);
    }

    //main propagation routine through trajectory

    propagate_computation(num_energy_levels, mass, x_min, x_max, total_time, num_points);

    return 0;

}

void propagate_computation(int num_energy_levels, long double mass, long double x_min,
			   long double x_max, long double total_time, int num_points){
  //files required:
  //	STATIS_POT (obtained)
  //	COORD_HOD (i.e. fort.20001)
  //  DIPIND_CMD (i.e. fort.30001)
  //  DIPMOL_CMD (i.e. fort.40001)

  long double step_size = 0.001;
  int num_intervals = (x_max - x_min) / step_size;

  std::string output_freq_filename = "freq_file.dat";
  std::string pot_filename = "fort.10001";

  std::ofstream export_freq;
  std::ifstream import_potential;

  export_freq.open(output_freq_filename.c_str());
  import_potential.open(pot_filename.c_str());

  // Naive implementation of the first line of the potential read-in


  long double scan_potential[num_points];
  long double scan_positions[num_points];

  //arrays needed for computation of eac arrays wavefunction:
  //scan pos: x-values, y-values = magnitude of potential

  read_in_frame(scan_positions, scan_potential, import_potential, num_points);

  //reset the string buffers to capture the next (read blank) line
  std::string blank_line;
  std::getline(import_potential, blank_line);

  spline potential_function(scan_positions, scan_potential, num_points);

  //Bring the "spline" to zero, sets the addition_adjustment parameter
  //within the spline object
  potential_function.adjust_zero(step_size, x_min, x_max);

  //potential_function.print_spline(step_size, x_min, x_max);

  long double A_param = 1.0 / std::sqrt(mass); 
  long double B_param = std::sqrt(A_param);

  long double red_E_guess=0.1L;
  int energy_level_being_computed=0;

  long double reduced_x_min = x_min / B_param;
  long double reduced_x_step = step_size / B_param;

  long double eigenvalues[num_energy_levels];
 
  std::vector<long double> eigenfunctions;
  
   //main body loop-compute x energy eigenfunctions
  while (energy_level_being_computed < num_energy_levels){

    long double eigenfunction[num_intervals];

    long double low_guess = 0.0L;
    long double high_guess = 2.0L;
    bool wavefunction_found=false;

    red_E_guess = (low_guess + high_guess) / 2.0;

	//for the given wavefunction in question, repeat the following
	//loop until a satisfactory solution is found
	while (wavefunction_found == false){
	          
	  //fill the rest of the wavefunction values in according to the
	  //"intelligent" guess for the reduced energy

	  eigenfunction[0] = 0.0L;
	  eigenfunction[1] = 0.00001L;

	  for (int i=0; i<num_intervals; ++i){

	    long double reduced_x = reduced_x_min + i*reduced_x_step;
	    long double prev_reduced_x = reduced_x - reduced_x_step;
	    long double next_reduced_x = reduced_x + reduced_x_step;			
	  
	    long double x = x_min + i*step_size;
	    long double prev_x = x - step_size;
	    long double next_x = x + step_size;			

        if (i > 1){

          long double pot_prev = potential_function.splint(prev_x) 
               			     + potential_function.addition_adjustment;
          long double pot_current = potential_function.splint(x) 
                                  + potential_function.addition_adjustment;
          long double pot_next = potential_function.splint(next_x)
                               + potential_function.addition_adjustment;
                
          long double prev_G_red = compute_reduced_G(prev_reduced_x, red_E_guess, 
                                                     pot_prev, A_param);
          long double current_G_red = compute_reduced_G(reduced_x, red_E_guess, 
                                                        pot_current, A_param);
          long double next_G_red = compute_reduced_G(next_reduced_x, red_E_guess, 
                                                     pot_next, A_param);
                
          eigenfunction[i] = compute_next_psi(eigenfunction[i-1],
                                              eigenfunction[i-2], current_G_red, 
                                              prev_G_red, next_G_red, reduced_x_step);
        }
      }

 	    //normalize wavefunction
 	  long double normalization = 0.0L;
  	  for (size_t l = 0; l < num_intervals; ++l)
        normalization += eigenfunction[l] * eigenfunction[l] * reduced_x_step;
	  	    
  	  normalization = std::sqrt(normalization);

  	  for (size_t l = 0; l < num_intervals; ++l)
        eigenfunction[l] /= normalization;

      bool endpoint_check = false;
      //boolean check to see if endpoint meets requirements after normalization
      if (eigenfunction[num_intervals - 1] * eigenfunction[num_intervals - 1] < 0.01)
        endpoint_check = true;

      int node_counter = 0;
	  //check number of nodes to see if the right eigenfunction has been found
	  for (size_t l = 0; l < num_intervals - 1; ++l)
	    if (eigenfunction[l] * eigenfunction[l+1] < 0) ++node_counter;

      if (node_counter > energy_level_being_computed){
        //guessed too high; high guess becomes reduced guess,
        //reduced guess becomes the average
        high_guess = red_E_guess;
        red_E_guess = (low_guess + high_guess) / 2.0;
      }

      else if ((node_counter < energy_level_being_computed) ||
        (node_counter == energy_level_being_computed && endpoint_check == false)){
        //guessed too low; low guess becomes reduced guess, 
        //...or the guess has the correct number of nodes but isn't converged
              
        // reduced guess reset to the average
        low_guess = red_E_guess;
        red_E_guess = (low_guess + high_guess) / 2.0;
      }  

      else if (node_counter == energy_level_being_computed && endpoint_check == true)
        wavefunction_found = true;

    }

//	      std::cout<<red_E_guess<<"\t" << red_E_guess * A_param << std::endl;
	  
    eigenvalues[energy_level_being_computed] = red_E_guess * A_param;
    ++energy_level_being_computed;
    for (int i = 0; i < num_intervals; ++i) eigenfunctions.push_back(eigenfunction[i]);

    }
    //write difference of eigenenergies to transition frequency file
  export_freq<<(eigenvalues[1] - eigenvalues[0])*hartree_to_wavenumbers<<'\t';
  export_freq<<(eigenvalues[2] - eigenvalues[1])*hartree_to_wavenumbers<<std::endl;

  for (int i=0; i<num_points; ++i){
    scan_potential[i] = 0.0;
    scan_positions[i] = 0.0;
  }
  
  export_freq.close();
}
