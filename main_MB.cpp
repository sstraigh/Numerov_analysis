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

const long double hartree_to_wavenumbers = 219474.63;
const long double kcal_to_hartree = 0.00159362;
const long double angstroms_to_bohr = 1.889725989;

void propagate_computation(int num_energy_levels, long double mass, long double x_min,
                           long double x_max, long double total_time);

//Computes the next value in the wavefunction using the current value [psi(x)]
//and the previous valu [psi(x-s)]
long double compute_next_psi(long double current_psi,   long double prev_psi,   long double current_G,   long double prev_G,   long double next_G,   long double s){

    long double next_psi= 2.0*current_psi - prev_psi;
    next_psi += (5.0 * current_G * current_psi * s * s / 6.0);
    next_psi += (prev_G * prev_psi * s * s /12.0);
    next_psi /= (1- (next_G * s * s / 12.0));

    return next_psi;

}

//computes the next G-value (for further reference and an overview of the
//Numerov algorithm, see Levine, Quatum Chemistry, 6th ed)
long double compute_reduced_G(  long double x_reduced, long double red_E_guess, 
		spline &func, long double A_param, long double B_param){

    long double x=x_reduced*B_param;
//    std::cerr<<"splining potential function"<<std::endl;
    long double potential=func.splint(x)+func.addition_adjustment;
//    std::cerr<<potential<<std::endl;
    potential/=A_param;
    return (2*potential-2*red_E_guess);
    
}

long double compute_G(long double x, long double E_guess, spline &func, long double mass){

//custom potential, computed using the spline interpolation function adapted from
//Numerical Recipes in C
    long double potential=(func.splint(x)+func.addition_adjustment);
    return (2*potential-2*E_guess)*mass;
}

long double compute_normalization_factor(long double * wavefunction, long double reduced_s,
	int num_intervals){

    long double reduced_total_sum=0.0;
    for (int i=0; i<num_intervals; ++i){
	reduced_total_sum+=(wavefunction[i]*wavefunction[i]*reduced_s);
    }

    return reduced_total_sum;

}

void read_in_frame(long double * scan_position, long double * reduced_dipole, long double * pot, 
	std::ifstream &potential, std::ifstream &coordinates, 
	std::ifstream &inddip, std::ifstream &permdip, int num_points,
	long double &total_dip_x, long double &total_dip_y, long double &total_dip_z,
	std::ifstream &import_NB_mu){

    for(int i=0; i<num_points; ++i){

	//readin relevant x-y values from data files needed
	//x-values the same for all files;

        std::string potential_line;
	std::string location_line;
	std::string permdip_line;
	std::string induceddip_line;
	std::string NB_mu_line;

	std::getline(potential, potential_line);
	std::getline(coordinates, location_line);
	std::getline(permdip, permdip_line);
	std::getline(inddip, induceddip_line);
	std::getline(import_NB_mu, NB_mu_line);

	std::istringstream potential_ss(potential_line);
	std::istringstream location_ss(location_line);
	std::istringstream permdip_ss(permdip_line);
	std::istringstream induceddip_ss(induceddip_line);
	std::istringstream NBdip_ss(NB_mu_line);

	//for potential, dealing with a two column file, x vs y
	potential_ss >> scan_position[i] >> pot[i];
	
	scan_position[i] *= angstroms_to_bohr; //convert from A to Bohr
	pot[i] *= kcal_to_hartree;

//	std::cout<<pot[i]<<std::endl;

	//for location, dealing with 11 cols:
	// t, scan_pos, Ox, Oy, Oz, Dx, Dy, Dz, Hx, Hy, Hz

	long double x0, x1, x2, y1, y2, y0, z0, z1, z2;
	long double temp;
	location_ss >> temp >> temp >> x0 >> y0 >> z0 >> x1 >> y1 >> z1 >> x2
		 >> y2 >> z2;
	long double red_x, red_y, red_z, normalization;
	red_x = x1 - x0;
	red_y = y1 - y0;
	red_z = z1 - z0;
	normalization = red_x*red_x + red_y*red_y + red_z*red_z;
	normalization = std::sqrt(normalization);
	red_x = red_x / normalization;
	red_y = red_y / normalization;
	red_z = red_z / normalization;

	long double total_dip, t_dip_x, t_dip_y, t_dip_z;
	NBdip_ss >> temp >> t_dip_x >> t_dip_y >> t_dip_z >> total_dip;

 	reduced_dipole[i] = t_dip_x*red_x+t_dip_y*red_y+t_dip_z*red_z;

	total_dip_x = t_dip_x;
	total_dip_y = t_dip_y;
	total_dip_z = t_dip_z;
/*
	//for dipole files, dealing with 9 cols: 
 		// t, scan_pos, x_cm, y_cm, z_cm, dip_x, dip_y, dip_z, mag
	long double ind_dip, perm_dip, x_cm, y_cm, z_cm, indx, indy, indz;
	long double total_dip, t_dip_x, t_dip_y, t_dip_z, permx, permy, permz;
	permdip_ss >> temp >> temp >> x_cm >> y_cm >> z_cm >> permx >> permy
		   >> permz >> perm_dip;
	induceddip_ss >> temp >> temp >> x_cm >> y_cm >> z_cm >> indx 
		      >> indy >> indz >> ind_dip;
	total_dip = ind_dip + total_dip;
	t_dip_x = indx + permx;
	t_dip_y = indy + permy;
	t_dip_z = indz + permz;
 	reduced_dipole[i] = t_dip_x*red_x+t_dip_y*red_y+t_dip_z*red_z;

//	std::cout<<reduced_dipole[i]<<std::endl;

	total_dip_x = t_dip_x;
	total_dip_y = t_dip_y;
	total_dip_z = t_dip_z;

*/	

	}

}

void read_command_line(char * arguments[], int &num_levels, long double &mass, long double &min_x, long double &max_x, long double &total_time){

    std::stringstream num_levels_ss(arguments[1]);
    num_levels_ss>>num_levels;

    std::stringstream mass_stream(arguments[2]);
    mass_stream>>mass;

    std::stringstream min_x_stream(arguments[3]);
    min_x_stream>>min_x;

    std::stringstream max_x_stream(arguments[4]);
    max_x_stream>>max_x;

    std::stringstream total_time_stream(arguments[5]);
    total_time_stream>>total_time;

  //  std::cout<<num_levels<<'\t'<<mass<<'\t'<<min_x<<'\t'<<max_x<<'\t'
//	     <<total_time<<std::endl;

}

int main(int argc, char * argv[]){

    if (argc<6){
	std::cerr<<"Usage: ./Numerov-Analysis num_energy_levels system_mass x_min x_max total_scan_time\n";
	return 1;
    }

//initialize variable holders: wavefunctions, values of wavefunctions

    int num_energy_levels=0;
    long double mass=0.0;
    long double x_min, x_max;
    long double total_time=0.0;
    read_command_line(argv, num_energy_levels, mass, x_min, x_max, total_time);

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

    long double step_size = 0.0001;
    int num_intervals = (x_max - x_min) / step_size;

    std::string output_freq_filename = "freq_file.dat";
    std::string output_tdip_filename = "transition_dipole.dat";

    std::string pot_filename = "STATIS_POT";
    std::string coord_filename = "fort.20001";
    std::string inddip_filename = "fort.30001";
    std::string perdip_filename = "fort.40001";
    std::string NB_mu_filename = "1B+NB_dipole_surface.dat";

    std::ofstream export_freq;
    std::ofstream export_tdip;

    std::ifstream import_potential;
    std::ifstream import_location;
    std::ifstream import_induceddip;
    std::ifstream import_permdip;
    std::ifstream import_NB_mu;

    export_freq.open(output_freq_filename.c_str());
    export_tdip.open(output_tdip_filename.c_str());

    import_potential.open(pot_filename.c_str());
    import_location.open(coord_filename.c_str());
    import_induceddip.open(inddip_filename.c_str());
    import_permdip.open(perdip_filename.c_str());
    import_NB_mu.open(NB_mu_filename.c_str());

// Naive implementation of the first line of the potential read-in
    std::string line;
    std::getline(import_potential, line);

    std::istringstream iss(line);
    int num_points;
    long double old_time=0.0;
    long double current_time=0.0;
    iss >> num_points >> current_time;
    if (iss.fail()){
	std::cerr << "cannot parse sstream: " << iss.str() << std::endl;
  	exit(1);
    }

//known; number of points per scan
//known; time of current scan
//known; time of total scans
   
    wavefunction * eigenfunctions = new wavefunction[num_energy_levels];

    long double * potential_at_scan = new long double[num_points];
    long double * scan_positions = new long double[num_points];
    long double * reduced_dipole = new long double[num_points];

    for (int i=0; i<num_energy_levels; ++i){
	eigenfunctions[i].initialize_wavefunction(x_min, x_max, step_size, mass);
    }

//    std::cout << "eigenfunctions initialized 1st time" << std::endl;

    while(current_time <= total_time){

	long double total_dip_x = 0.0;
	long double total_dip_y = 0.0;
	long double total_dip_z = 0.0;

//	std::cout<<"computing frame : "<<current_time<<std::endl;

	//arrays needed for computation of eac arrays wavefunction:
	//scan pos: x-values, y-values = magnitude of potential

	//<u> needs total dip (x,y,z) and ind+perm dip(x,y,z)
		//will need to interpolate these

	read_in_frame(scan_positions, reduced_dipole, potential_at_scan,
		import_potential, import_location, import_permdip,
		import_induceddip, num_points, total_dip_x, total_dip_y,
		total_dip_z, import_NB_mu);

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

//	std::cout << "splining dip and pot" << std::endl;

	spline potential_function(scan_positions, potential_at_scan, num_points);
	spline dipole_function(scan_positions, reduced_dipole, num_points);

        //Bring the "spline" to zero, sets the addition_adjustment parameter
        //within the spline object
    	potential_function.adjust_zero(step_size, x_min, x_max);

        long double red_E_guess=0.1L;
        long double increment=0.1L;
        int energy_level_being_computed=0;

	long double reduced_x_min = 0.0;
	long double reduced_x_step = 0.0;

        //main body loop-compute x energy eigenfunctions
        while (energy_level_being_computed<num_energy_levels){

//	    eigenfunctions[energy_level_being_computed].initialize_wavefunction(x_min, x_max, step_size,
  //                                                    mass);

	//    std::cout<<"Computing energy level "<< (energy_level_being_computed+1)
	//	     <<" of "<<num_energy_levels<<std::endl;

	    bool wavefunction_found=false;

	    int node_counter = 0;

            increment=0.1L;
	    red_E_guess+=increment;

	    //for the given wavefunction in question, repeat the following
	    //loop until a satisfactory solution is found
	    while (wavefunction_found == false){

//		std::cout<<"attempting init..."<<std::endl;
//		eigenfunctions[energy_level_being_computed].initialize_wavefunction(x_min, x_max, step_size,
  //                                                    mass);

//		std::cout<<" Calling eigenfunctions for xmin and xstep.."<<std::endl;

		reduced_x_min = eigenfunctions[energy_level_being_computed].get_reduced_x_min();
		reduced_x_step = eigenfunctions[energy_level_being_computed].get_reduced_x_step();

	        //fill the rest of the wavefunction values in according to the
	        //"intelligent" guess for the reduced energy
	        for (int i=0; i<num_intervals; ++i){

		    long double reduced_x = reduced_x_min + i*reduced_x_step;
		    long double prev_reduced_x = reduced_x - reduced_x_step;
		    long double next_reduced_x = reduced_x + reduced_x_step;			

		    long double A_param = eigenfunctions[energy_level_being_computed].get_A_param();
		    long double B_param = eigenfunctions[energy_level_being_computed].get_B_param();

//		    std::cout<<reduced_x<<std::endl;

		    if (i > 1){

	                long double prev_G_red=compute_reduced_G(prev_reduced_x,
			    red_E_guess, potential_function, A_param, B_param);

                        long double current_G_red=compute_reduced_G(reduced_x,
      			    red_E_guess, potential_function, A_param, B_param);

                        long double next_G_red=compute_reduced_G(next_reduced_x, 
			    red_E_guess, potential_function, A_param, B_param);
                        eigenfunctions[energy_level_being_computed].compute_next_psi(
			    current_G_red, prev_G_red, next_G_red, i);
		//	std::cout<<"done with pot"<<std::endl;
              	    }
                }

		wavefunction_found = eigenfunctions[energy_level_being_computed].check_wavefunction(energy_level_being_computed);

		node_counter = eigenfunctions[energy_level_being_computed].get_node_num();

	        //If we guessed to high (which can only happen after a previous guess was too
	        //low), then reduce the energy guess to what it was previously and divide the
	        //increment by 10. This "intelligently" narrows down to a more correct energy
	        //guess. For example, if the eigenvalue is 0.315, then the sequence of guesses
	        //proceeds as following:
	        //	0.1 (too low, guess=  guess + increment)
	        //	0.2 (too low, guess=  guess + increment)
	        //	0.3 (too low, guess=  guess + increment)
	        //	0.4 (too hi, guess=  guess + increment, increment/=10)
	        //	0.3 (too low, guess=  guess + increment)	 
	        //	0.31 (too low, guess= guess+increment)
	        //	0.32, increment /=10,
	        //	...0.31, 0.311, 0.312, 0.313, 0.314, 0.315
	        if (node_counter>energy_level_being_computed && wavefunction_found == false){
		    red_E_guess-=increment;
		    increment/=10.0L;
		    eigenfunctions[energy_level_being_computed].deinitialize();
	        }

	        else if (node_counter<=energy_level_being_computed && wavefunction_found == false){
		    red_E_guess+=increment;
//		    std::cout<<eigenfunctions[0].get_psi_at_index(100)<<std::endl;
		    eigenfunctions[energy_level_being_computed].deinitialize();
//		    std::cout<<eigenfunctions[energy_level_being_computed].get_A_param()<<std::endl;
	        }
//		std::cout<<energy_level_being_computed<<'\t'<<red_E_guess<<std::endl;
//		std::cout<<red_E_guess<<std::endl;
	    }

//	    std::cout<<red_E_guess<<std::endl;

	    eigenfunctions[energy_level_being_computed].set_energy(red_E_guess, potential_function.addition_adjustment);
	    //Print Energy levels and eigenvalues to the screen
//	    std::cout<<wave_energy<<'\t'<<potential_function.addition_adjustment<<'\t'<<A_param<<std::endl;

	    ++energy_level_being_computed;
//  	    std::cout<<"Computing energy level : " <<energy_level_being_computed;      

        }

//	std::cout<<"accessing energies..."<<std::endl;

	export_freq<<old_time<<'\t';
	export_freq<<(eigenfunctions[1].energy - eigenfunctions[0].energy)*hartree_to_wavenumbers<<'\t';
	export_freq<<(eigenfunctions[2].energy - eigenfunctions[1].energy)*hartree_to_wavenumbers<<std::endl;

//	std::cout<<"energies written"<<std::endl;

//	std::cout<<dipole_function.splint(1.0);

	export_tdip<<old_time<<'\t';
	long double B_param = eigenfunctions[0].get_B_param();

//	std::cout<<dipole_function.splint(1.0)<<std::endl;

	for (int i=0; i<(num_energy_levels-1); ++i){
		
		long double transition_dipole_moment=0.0;

		for (int j=0; j<num_intervals; ++j){

			long double x_value = reduced_x_min + (j*reduced_x_step);

//			std::cout<<x_value<<'\t'<<reduced_x_min<<'\t'<<reduced_x_step<<'\t'<<B_param<<std::endl;

			long double dipole_at_x = dipole_function.splint((x_value*B_param));
//			std::cout<<dipole_at_x<<std::endl;
			transition_dipole_moment += eigenfunctions[i].get_psi_at_index(j)
				* eigenfunctions[i+1].get_psi_at_index(j) * dipole_at_x 
				* reduced_x_step;			
		}
		export_tdip<<transition_dipole_moment<<'\t';
	
	}	
	for (int i=0; i<num_energy_levels; ++i) eigenfunctions[i].deinitialize();
	export_tdip<<total_dip_x<<'\t'<<total_dip_y<<'\t'<<total_dip_z<<std::endl;

//	std::cout<< "eigenfunctions de-initialized" << std::endl;

//	potential_function.addition_adjustment = 0;
	for (int i=0; i<num_points; ++i){
	
	    potential_at_scan[i]= 0.0;
	    scan_positions[i] = 0.0;
	    reduced_dipole[i] = 0.0;

	}

    }

//    for (int i=0; i<num_energy_levels; ++i) eigenfunctions[i].free_holder();
    delete[] eigenfunctions;

    export_freq.close();
    export_tdip.close();

}
