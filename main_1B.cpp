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
                           long double x_max, long double total_time);

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

void read_in_frame(long double * scan_position, long double * reduced_dipole, long double * pot, 
	std::ifstream &potential, std::ifstream &coordinates, 
	std::ifstream &inddip, std::ifstream &permdip, int num_points,
	long double &total_dip_x, long double &total_dip_y, long double &total_dip_z){

    for(int i=0; i<num_points; ++i){

	//readin relevant x-y values from data files needed
	//x-values the same for all files;

        std::string potential_line;
	std::string location_line;
	std::string permdip_line;
	std::string induceddip_line;

	std::getline(potential, potential_line);
	std::getline(coordinates, location_line);
	std::getline(permdip, permdip_line);
	std::getline(inddip, induceddip_line);

	std::istringstream potential_ss(potential_line);
	std::istringstream location_ss(location_line);
	std::istringstream permdip_ss(permdip_line);
	std::istringstream induceddip_ss(induceddip_line);

	//for potential, dealing with a two column file, x vs y
	potential_ss >> scan_position[i] >> pot[i];
	
	scan_position[i] *= angstroms_to_bohr; //convert from A to Bohr
	pot[i] *= kcal_to_hartree;

	//for location, dealing with 11 cols:
	// t, scan_pos, Ox, Oy, Oz, Dx, Dy, Dz, Hx, Hy, Hz

	long double x0, x1, x2, y1, y2, y0, z0, z1, z2;
	long double temp;
	location_ss >> temp >> temp >> x0 >> y0 >> z0 >> x1 >> y1 >> z1 >> x2
		 >> y2 >> z2;

	//compute the OD unit vector
	long double red_x, red_y, red_z, normalization;
	red_x = x1 - x0;
	red_y = y1 - y0;
	red_z = z1 - z0;
	normalization = red_x*red_x + red_y*red_y + red_z*red_z;
	normalization = std::sqrt(normalization);
	red_x = red_x / normalization;
	red_y = red_y / normalization;
	red_z = red_z / normalization;

	//read-in fort.30001 and fort.40001, output from DLPOLY scan
	//this is the 1B+NB dipole surface

	//col format: t, scan_pos, x_cm, y_cm, z_cm, x_dip, y_dip, z_dip, dip_mag
	long double x_cm, y_cm, z_cm;
	long double x_inddip, y_inddip, z_inddip;
	long double x_prmdip, y_prmdip, z_prmdip;
	long double prm_mag, ind_mag;

	induceddip_ss >> temp >> temp >> x_cm >> y_cm >> z_cm 
		      >> x_inddip >> y_inddip >> z_inddip >> ind_mag;
	permdip_ss >> temp >> temp >> temp >> temp >> temp
		   >> x_prmdip >> y_prmdip >> z_prmdip >> prm_mag;

	long double dip_mag, t_dip_x, t_dip_y, t_dip_z;

	t_dip_x = /* x_inddip*/ + x_prmdip;
	t_dip_y = /* y_inddip*/ + y_prmdip;
	t_dip_z = /* z_inddip*/ + z_prmdip;
	dip_mag = ind_mag + prm_mag;

 	reduced_dipole[i] = t_dip_x * red_x + t_dip_y * red_y + t_dip_z * red_z;

	total_dip_x += red_x;
	total_dip_y += red_y;
	total_dip_z += red_z;

    }

    total_dip_x /= num_points;
    total_dip_y /= num_points;
    total_dip_z /= num_points;

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

    if (argc != 1){
	std::cerr<<"Usage: ./Numerov-Analysis (read input from INPUT.dat)\n";
	return 1;
    }

    //initialize variable holders: wavefunctions, values of wavefunctions

    int num_energy_levels=0;
    long double mass=0.0;
    long double x_min, x_max;
    long double total_time=0.0;
    read_input_file(num_energy_levels, mass, x_min, x_max, total_time);

    //main propagation routine through trajectory

    propagate_computation(num_energy_levels, mass, x_min, x_max, total_time);

    return 0;

}

void propagate_computation(int num_energy_levels, long double mass, long double x_min,
			   long double x_max, long double total_time){
    //files required:
    //	STATIS_POT (obtained)
    //	COORD_HOD (i.e. fort.20001)
    //  DIPIND_CMD (i.e. fort.30001)
    //  DIPMOL_CMD (i.e. fort.40001)

    long double step_size = 0.001;
    int num_intervals = (x_max - x_min) / step_size;

    std::string output_freq_filename = "freq_file.dat";
    std::string output_tdip_filename = "transition_dipole.dat";

    std::string pot_filename = "STATIS_POT";
    std::string coord_filename = "fort.20001";
    std::string inddip_filename = "fort.30001";
    std::string perdip_filename = "fort.40001";

    std::ofstream export_freq;
    std::ofstream export_tdip;

    std::ifstream import_potential;
    std::ifstream import_location;
    std::ifstream import_induceddip;
    std::ifstream import_permdip;

    export_freq.open(output_freq_filename.c_str());
    export_tdip.open(output_tdip_filename.c_str());

    import_potential.open(pot_filename.c_str());
    import_location.open(coord_filename.c_str());
    import_induceddip.open(inddip_filename.c_str());
    import_permdip.open(perdip_filename.c_str());

    // Naive implementation of the first line of the potential read-in
    std::string line;
    std::getline(import_potential, line);

    std::istringstream iss(line);
    int num_points;
    long double old_time=0.0;
    long double current_time=0.0;
    iss >> num_points >> current_time;
    if (iss.fail()){
	std::cerr << "cannot parse sstream: " << iss << std::endl;
  	exit(1);
    }

    //known; number of points per scan
    //known; time of current scan
    //known; time of total scans

    long double scan_potential[num_points];
    long double scan_positions[num_points];
    long double reduced_dipole[num_points];

    while(current_time <= total_time){

	long double total_dip_x = 0.0;
	long double total_dip_y = 0.0;
	long double total_dip_z = 0.0;

	//arrays needed for computation of eac arrays wavefunction:
	//scan pos: x-values, y-values = magnitude of potential

	read_in_frame(scan_positions, reduced_dipole, scan_potential,
		import_potential, import_location, import_permdip,
		import_induceddip, num_points, total_dip_x, total_dip_y,
		total_dip_z);

	std::string blank_line;
	//reset the string buffers to capture the next (read blank) line
	old_time = current_time;

	std::getline(import_potential, blank_line);
        std::istringstream iss(blank_line);
	long double temp;
	iss >> temp >> current_time;
	
	blank_line.clear();
	std::getline(import_location, blank_line);
	blank_line.clear();
	std::getline(import_permdip, blank_line);
	blank_line.clear();
	std::getline(import_induceddip, blank_line);

	spline potential_function(scan_positions, scan_potential, num_points);
	spline dipole_function(scan_positions, reduced_dipole, num_points);

        //Bring the "spline" to zero, sets the addition_adjustment parameter
        //within the spline object
    	potential_function.adjust_zero(step_size, x_min, x_max);

//    	potential_function.print_spline(step_size, x_min, x_max);

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

/*	        std::cout<< red_E_guess << "\t" << red_E_guess * A_param * hartree_to_wavenumbers
			 << "\t" << node_counter << "\t" << energy_level_being_computed << std::endl;*/

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

//	    std::cout<<red_E_guess<<"\t" << red_E_guess * A_param << std::endl;
	
	    eigenvalues[energy_level_being_computed] = red_E_guess * A_param;

	    ++energy_level_being_computed;

	    for (int i = 0; i < num_intervals; ++i) eigenfunctions.push_back(eigenfunction[i]);

        }

	//write difference of eigenenergies to transition frequency file
	export_freq<<old_time<<'\t';
	export_freq<<(eigenvalues[1] - eigenvalues[0])*hartree_to_wavenumbers<<'\t';
	export_freq<<(eigenvalues[2] - eigenvalues[1])*hartree_to_wavenumbers<<std::endl;

	export_tdip<<old_time<<'\t';

	for (int i=0; i<(num_energy_levels-1); ++i){
		
		long double transition_dipole_moment=0.0;

		for (int j=0; j<num_intervals; ++j){

			long double x_value = reduced_x_min + (j*reduced_x_step);

			long double dipole_at_x = dipole_function.splint((x_value*B_param));
//			std::cout<<dipole_at_x<<std::endl;
			transition_dipole_moment += 1/std::sqrt(B_param) * eigenfunctions[i*num_intervals + j]
				* eigenfunctions[(i+1) * num_intervals + j] * 1/std::sqrt(B_param) * dipole_at_x 
				* step_size;			
/*			std::cerr << (x_min + j*step_size) << '\t' << dipole_at_x << '\t' << 1/std::sqrt(B_param) * eigenfunctions[i*num_intervals + j]
				  << '\t' << 1/std::sqrt(B_param) * eigenfunctions[(i+1)*num_intervals + j] << std::endl;*/

		}
		export_tdip<<transition_dipole_moment<<'\t';
	
	}	
	export_tdip << total_dip_x << '\t' << total_dip_y << '\t' <<total_dip_z <<std::endl;

	for (int i=0; i<num_points; ++i){
	
	    scan_potential[i] = 0.0;
	    scan_positions[i] = 0.0;
	    reduced_dipole[i] = 0.0;

	}

    }

    export_freq.close();
    export_tdip.close();

}
