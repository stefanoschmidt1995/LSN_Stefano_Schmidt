/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 7
*************************************************************************/
#include "../particle_simulator/particle_simulator.h"
#include "../standard_header.h"

void EmptyFiles(std::string phase_type){
	std::ofstream ofs;
	std::vector<std::string> files_to_read = {"./out/p_"+phase_type+".out","./out/e_"+phase_type+".out", "./out/avg_"+phase_type+".dat"};
	for(auto it = files_to_read.begin(); it!= files_to_read.end(); it++){
		ofs.open(*it, std::ios::trunc);
		ofs << std::endl;
		ofs.close();
 	}
  return;
}

int
main (int argc, char *argv[]){

	std::string phase_type= "liquid"; //phase to simulate (liquid default)
	if (argc>1)
		phase_type = argv[1];

		//************* EX 07.1/07.4 - AUTOCORRELATION, DATA BLOCKING, SIMULATION & HISTOGRAM
		//This code computes in one single simulation many things:	istant values of e/p
		//								averaged data as a function of time
		//								histogram of g(r)
		//								plot of errors as a function L (date per block)

	int n_step= 29000, eq_steps=1000; //n_step + eq_step = 5e4
	EmptyFiles(phase_type);
	particle_simulator simulator("input."+phase_type);

	std::cout << std::endl << "----STARTING SIMULATION----"<<std::endl;
	std::vector<double> params = simulator.get_config_parameters();
	std::cout << "\tTemp: " << params[1] << std::endl;
	std::cout << "\trho: " << params[3] << std::endl;
	std::cout << "\tr_cut: " << params[4] << std::endl <<std::endl;
	std::cout << "Writing to files: \n\t./out/e_"+phase_type+".out; ./out/p_"+phase_type+".out\n\t./out/avg_eerr_"+phase_type+".out; ./out/avg_perr_"+phase_type+".out\n\t./out/avg_"+phase_type+".dat \n\t./out/hist_NVT_"+phase_type+".dat" <<std::endl;
	int measure_step = params[6];

	for (int i =0; i< n_step+eq_steps; i++){
		simulator.move();
		simulator.measure();
		simulator.write_observable("./out/p_"+phase_type+".out", 5); //writing istantaneous value of pressure
		simulator.write_observable("./out/e_"+phase_type+".out", 2); //writing istantaneuos value of energy
			
			//equilibration part is over, measure buffer is cleared to compute averages
		if (i == eq_steps) simulator.reset_measures(); 
		if(i%measure_step == 0 && i >= 2*measure_step+eq_steps && i >eq_steps) //computing averages
			simulator.write_avg_observable("./out/avg_"+phase_type+".dat");
		if (i%1000==0) std::cerr <<"Done step: "<< i<<std::endl; //communicate with user
	}

	simulator.write_hist("./out/hist_NVT_"+phase_type+".dat"); //saving histogram

		//preparing files avg_p.out and avg_e.out with uncertainties as a function of L
	std::fstream e_file;
	std::fstream p_file;
	e_file.open("./out/e_"+phase_type+".out", std::ios::in);
	p_file.open("./out/p_"+phase_type+".out",std::ios::in);
	std::vector<double> e_vect(n_step), p_vect(n_step);
	for (int i =0; i< n_step; i++){
		e_file >> e_vect.at(i);
		p_file >> p_vect.at(i);
	}
	e_file.close();
	p_file.close();
	e_file.open("./out/avg_eerr_"+phase_type+".out",std::ios::out);
	p_file.open("./out/avg_perr_"+phase_type+".out",std::ios::out);
	int n_blocks = 10;
	for(int L = 10; L< (int)2e3; L++){
		std::vector<double> e_res = compute_statistical_error(e_vect, n_blocks, L);
		std::vector<double> p_res = compute_statistical_error(p_vect, n_blocks, L);
		e_file << L << "\t" << e_res[2] <<std::endl;
		p_file << L << "\t" << p_res[2] <<std::endl;
	}
	e_file.close();
	p_file.close();

	return 0;
}







