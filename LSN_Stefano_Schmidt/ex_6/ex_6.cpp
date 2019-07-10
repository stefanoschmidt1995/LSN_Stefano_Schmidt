/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 6
*************************************************************************/
#include "../standard_header.h"
#include "ising_simulator.h"

int 
main(){

		//************* EX 06.1 - ISING MODEL
	//deifining some parameters for simulation
	int n_step= 20000, eq_steps=1000; //n_step = steps of simulation at each temperature; eq_steps = equilibration steps
	int sparse_avg_step = 4;
	int temp_steps = 50;

	double T_min = 0.1, T_max = 5.;

	//defining the simulator and some stuff related to ouput files
	ising_simulator simulator("input.dat");
	std::vector<double> params = simulator.get_config_parameters();
	std::string file_prefix = "./out/gibbs_";
	if(params[4] == 1) file_prefix = "./out/met_";

	std::ofstream write_energy, write_heat, write_chi, write_mag;
	write_energy.open(file_prefix+"energy.dat");
	write_heat.open(file_prefix+"heat.dat");
	if(params[3]!=0) write_mag.open(file_prefix+"mag.dat"); //write only if h is nonzero
	write_chi.open(file_prefix+"chi.dat");

	//starting actual simulation
	std::cout <<"Starting simulation of Ising Model"<<std::endl;
	for (int t = 0; t <= temp_steps; t++){
		simulator.set_spin_configuration("initial.config");
		simulator.reset_measures();
		double T = T_min+t*(T_max-T_min)/(double)temp_steps;
		simulator.set_temperature(T);
		std::cout << "**Simulating system at temperature "<< T << std::endl;
		for(int i=0; i< eq_steps; i++) //equilibration
			simulator.move();

		for(int i=0; i< sparse_avg_step*n_step; i++){ //actual simulation
			simulator.move();
			if (i%sparse_avg_step==0)
				simulator.measure();
		}
		std::vector<avg_type> averages = simulator.get_measures(); //gets the saved measures averaged with L specified in input file
		std::cout << "\tInternal energy: "<< averages[0][0] << " " << averages[0][2] << std::endl;
		std::cout << "\tSpecific heat: "<< averages[1][0] << " " << averages[1][2] << std::endl;
		std::cout << "\tMagnetization: "<< averages[2][0] << " " << averages[2][2] << std::endl;
		std::cout << "\tSusceptivity: "<< averages[3][0] << " " << averages[3][2] << std::endl;

			//handling file printing
		write_energy << T << "\t"<<averages[0][0] << "\t" << averages[0][2] << std::endl;
		write_heat << T << "\t"<<averages[1][0] << "\t" << averages[1][2] << std::endl;
		if(params[3]!=0) //write only if h is nonzero
			write_mag << T << "\t"<<averages[2][0] << "\t" << averages[2][2] << std::endl;
		write_chi << T << "\t"<<averages[3][0] << "\t" << averages[3][2] << std::endl;

	}
	write_energy.close();
	write_heat.close();
	if(params[3]!=0) write_mag.close();
	write_chi.close();
	return 0;
}
