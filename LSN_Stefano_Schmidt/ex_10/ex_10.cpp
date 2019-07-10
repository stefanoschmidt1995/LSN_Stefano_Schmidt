/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 10
*************************************************************************/
#include "../standard_header.h"
#include "../GA_helper/GA_helpers.h"

bool accept_metropolis_move(individual_type*, individual_type*, double, Random*);

int
main(){
		//************* EX 10.1 - SIMULATED ANNEALING
	int N = 30; //path lenght
	double l = 1.; //size of the box
	int count =0; //to count how many iterations without a change in path lenght
	int N_max = 100; //maximum number of useless iteration to have before stopping the algorithm
	int n_iter = 0; //keep track of the iteration that are done...
	int N_step = 1000; //# of step at every temperature
	double c = .05;
	double T=1000;

	Random rng = initialize_random();	

	destinations_list SLT_path = generate_destinations_box(N,l, &rng);
	std::string best_path_filename = "./out/best_path_box.dat";
	//destinations_list SLT_path = generate_destinations_circle(N,l, &rng);
	//std::string best_path_filename = "./out/best_path_circle.dat";	

	individual_type individual(N, &SLT_path);

	std::ofstream distance_file("./out/distances.dat");

	double best_lenght = individual.path_lenght(), old_lenght=0;
	individual_type best_individual(individual);
	
		//starting with the updates in temperatures
	while(count<N_max){
		double beta = 1./T;
		std::cout << "Sistem at temperature "<< T << std::endl;
		for(int i =0; i< N_step; i++){//doing a N_step metropolis algorithm
			individual_type try_individual(individual);
			try_individual.perform_mutation(&rng);
			if(accept_metropolis_move(&individual, &try_individual, beta, &rng)){
				individual = try_individual;
				if (individual.path_lenght() < best_lenght){ //if a new best individual is found
					best_lenght = individual.path_lenght();
					best_individual = individual;
				}
			}
		}
		distance_file << n_iter << "\t"<< T << "\t" << best_lenght <<std::endl; //printing out current best lenght
		if (old_lenght == best_lenght)
			count++;
		else{
			count =0;
			old_lenght = best_lenght;
		}
		T= T/(c+1);
		n_iter++;
	}
	distance_file.close();
	std::cout << "Best path " << best_individual.path_lenght() << std::endl;
	std::cout << "Path is " <<std::endl;	
	individual.print();
	individual.print_path_to_file(best_path_filename);

	return 0;
}


bool
accept_metropolis_move(individual_type* individual, individual_type* new_individual, double beta, Random* rng){
	double old_cost = individual->path_lenght();
	double new_cost = new_individual->path_lenght();
	double acceptance = exp(-beta*(new_cost-old_cost));
	if(acceptance>=1)
		return true;
	else{
		double rand = rng->Rannyu();
		return (rand<=acceptance);
	}
}















