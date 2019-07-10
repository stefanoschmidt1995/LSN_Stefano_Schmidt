/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
CLASS ISING SIMULATOR - DECLARATIONS
**************************************************************************
Class to sample from configuration of 1D Ising model at fixed temperature using sampling with Metropolis algorithm or Gibbs sampling
The hamiltonian is H = - J Sum s_i*s_j - h * Sum s_i
*************************************************************************/
#ifndef __ising_sim_
#define __ising_sim_

#include "../standard_header.h"

typedef std::vector<double> avg_type;

class
ising_simulator{
	public:
		ising_simulator(std::string filename);//constructor: takes in input a configuration file
		~ising_simulator(){};
		void config_from_file(std::string filename);
		void set_spin_configuration(std::string filename); //put in the state vector a new spin configuration
		void write_to_file(std::string filename); //write to file the current state of the system
		void reset_measures(); //clears the buffer of measures
		void set_temperature(double new_temperature); //set new temperature to the Ising system
		void move(); //create a new sample of ising model
		void measure(); //add the measure of configuration to the whole set of measures 
		std::vector<avg_type> get_measures(); //return a vector with all the averages
		std::vector<double> get_config_parameters(); //returns the configuration parameters: [temp, N, J, h, metropolis,n_meas]
		avg_type get_autocorrelation(unsigned int time_step, int N_meas = 20); //compute the autocorrelation after time_step
	private:
		int Pbc(int i); //periodic boundary conditions
		double get_energy_difference(int spin_index); //returns the change in energy if spin spin_index is flipped
		int draw_from_conditional(int spin_index);
			//draw a spin orientation for spin of index k=spin_index drawn from p(s_k|s_i,i!=k)

	//data members
		int n_props = 4; //n_properties
		int N; //number of spin
		int moves=0, accepted = 0; 
		std::vector<int> s; //holds current spin configurations
		double beta,temp,J,h; //ising hamiltonian parameters
		bool metropolis; //true metropolis; false gibbs sampling
		int n_meas; //number of measures in each block
		std::vector<std::list<double>> measure_vector; //[[energy], [heat], [mag], [chi]]
		Random rng; //random generator
		

};

#endif
