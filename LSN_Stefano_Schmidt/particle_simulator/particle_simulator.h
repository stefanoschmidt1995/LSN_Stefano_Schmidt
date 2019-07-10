#ifndef __particle_simulator_class_
#define __particle_simulator_class_

//class meant to sample from configuration of a system of particle interactia via LJ potential.
//The thermodinamical system can be defined either by NVE either by NVT thermodinamical variables 

#include"../standard_header.h"

typedef std::vector<double> avg_type;

class
particle_simulator{

	public:
		particle_simulator(std::string filename, std::string initial_conf = "config.0", std::string intial_old ="");
		~particle_simulator();
		void config_from_file(std::string filename);		//configures the simulator from file setting the relevant params
		void set_particle_configuration(std::string filename, std::string intial_old ="");
									//put the current state in a new particle configuration
		void write_to_file(std::string filename, std::string filename_old ="");
									//writes to file the current state of the system
		void write_observable(std::string filename, unsigned int index);
									//append to file the istantaneuos value of the observable 										corresponding to index
		void write_avg_observable(std::string filename);	//write the averages of all the observable to file filename
		void write_hist(std::string filename);			//write averaged histogram to file with uncertainties
		void reset_measures(); 					//cleans the buffer of measures
		void move(); 						//makes a move (Verlet integration in NVE and Metropolis in NVT)
		void measure(); 					//adds the measure of configuration to the set of measures 
		std::vector<avg_type> get_measures(); 			//returns a vector with all the averages
		std::vector<std::vector<double>> get_g_hist(); 		//returns averaged g_hist with format
									//[r_i, h_i, simga_i, std_dev_i] with i being the i-th bin
		std::vector<double> get_config_parameters();		//return config parameters in the same order of config_from_file
		void set_new_temperature(double new_temp);		//set a new temperature for the system
		double get_success_rate();				//for NVT returns the success rate of the Metropolis moves
		bool find_eq_temperature(double temp, int n_iter);	//in NVE looks iteratively for a system with an avg kinetic 										energy matching the required temperature temp within the 										statistical error obtained averaging n_iter measures
	private:
		void move_NVT(); 
		void move_NVE(); 
		double get_temperature();				//for private use only... 
		double Pbc(double r);					//periodic boundary conditions
		double particle_energy(R3_point x, int ip); 		//return the energy of particle index ip as if it was located 										in position x
		R3_point force(int ip);					//computes the force as -grad_ip V(x) required for integration 										of EOM in NVE


		//data members
		std::vector<R3_point> state;			//holds the position of the particles
		std::vector<R3_point> old_state, velocities;	//hold the previous state and the particle's velocities (only in NVE)
		unsigned int n_props = 6;
		int N; 						//number of particles
		int attempted=0, accepted = 0; 
		double beta,temp,vol, rho, box, rcut; 		//model parameters
		std::string type;				//hold the type of system (NVE or NVT are implemented)
		int n_meas; 					//number of measure to put in each block when averaging
		std::vector<std::list<double>*> measure_vector;	//[[E_tot], [E_kin], [E_pot], [virial], [Temp], [P]]
		std::list<std::vector<int>*> hist_list;		//holds different g_hist measured at many times
		Random rng;
		double vtail,ptail; 				//correction for finite size box to pressure and potential energy
		double delta; 					//time step for NVE; spatial step for a metropolis move in NVT
		
		//std::vector<int> g_hist; 			//hold the histogram of counting for g
		int n_bins;					//number of bins to be used for the computation of g
		double bin_size;				//size of a bin (dimension of a distance)
		int moves = 0;					//keep track of how many updates has experienced the system
};


#endif







