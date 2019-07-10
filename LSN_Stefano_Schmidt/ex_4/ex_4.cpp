/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 4
*************************************************************************/
#include "../particle_simulator/particle_simulator.h"
#include "../standard_header.h"

void EmptyFiles(){//function to empty some output files before using them
	std::ofstream ofs;
	std::string folder = "./out/";
	std::vector<std::string> files_to_read = {"output_epot.dat","output_ekin.dat","output_temp.dat","output_etot.dat","output_temp.dat", "output_press.dat"};
	for(auto it = files_to_read.begin(); it!= files_to_read.end(); it++){
		ofs.open(folder+(*it), std::ios::trunc);
		ofs << std::endl;
	ofs.close();
	}
	return;
}

int
main (int argc, char *argv[]){

		//************* EX 04.4 - LJ SIMULATION
	bool set_temp = true;
	double temp_to_set;
	int n_step = 5000;
	int iprint =100;

	std::string phase_type= "liquid";
	if (argc>1)
		phase_type = argv[1];

	particle_simulator simulator("input."+phase_type);
	std::ofstream initialize_file("./out/avg_"+phase_type+".dat");
	initialize_file << "[rho, [E_tot], [E_kin], [E_pot], [virial], [Temp], [P], time]" << std::endl;
	initialize_file.close();

	std::vector<double> params = simulator.get_config_parameters();
	int measure_step =  params[6]; //a measure is performed every n_meas step
	temp_to_set = params[1];
	if (set_temp) //finding eq temperature for the system
		if(!simulator.find_eq_temperature(temp_to_set, 50)) return -1;
	simulator.reset_measures();

		//starting the actual simulation
	std::cout << std::endl << "----STARTING SIMULATION----"<<std::endl;
	params = simulator.get_config_parameters();
	std::cout << "\tTemp: " << params[1] << std::endl;
	std::cout << "\trho: " << params[3] << std::endl;
	std::cout << "\tr_cut: " << params[4] << std::endl <<std::endl;

	for(int istep=1; istep <= n_step; ++istep){
		simulator.move();
		simulator.measure();
		if(istep%iprint == 0)
			std::cout << "Number of time-steps: " << istep << std::endl;
		if(istep%measure_step == 0 && istep >= 2*measure_step) //if there is a new block, data are averaged and saved
			simulator.write_avg_observable("./out/avg_"+phase_type+".dat");

			//snippet to print to file istanteneous values of different observables (not used here, but might be useful)
		/*if(istep%measurestep == 0){
			simulator.measure();
			/*simulator.write_observable("./out/output_epot.dat",2); 
			simulator.write_observable("./out/output_ekin.dat",1); 
			simulator.write_observable("./out/output_etot.dat",0); 
			simulator.write_observable("./out/output_temp.dat",4);
			simulator.write_observable("./out/output_press.dat",5);
		}*/
	}
	simulator.write_hist("./out/hist_NVE_"+phase_type+".dat");
	return 0;
}








