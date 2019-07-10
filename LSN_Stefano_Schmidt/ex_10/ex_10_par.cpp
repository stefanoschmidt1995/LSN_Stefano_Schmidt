/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 10 - PARALLEL CODE
*************************************************************************/
#include "../standard_header.h"
#include "../GA_helper/GA_helpers.h"
#include "mpi.h"

//to run the code: mpiexec -np 10 ./ex_10_par.exe 

bool accept_metropolis_move(individual_type*, individual_type*, double, Random*);

int
main(int argc, char* argv[]){

		//************* EX 10.2 - PARALLEL SIMULATED ANNEALING
	MPI::Init(argc,argv);

	int size = MPI::COMM_WORLD.Get_size(); //getting the size of the parallel computing system 
	int rank = MPI::COMM_WORLD.Get_rank(); //getting the tid of the current knot

	int N = 30; //path lenght
	double l = 1.; //size of the box
	int count = 0; //to count how many iterations without a change in path lenght
	int N_max = 100; //maximum number of iteration to have before stopping the algorithm
	int n_iter = 0; //keep track of the iteration that are done...
	int N_step = 1000; //# of step at every temperature
	double c = .05;
	double T=1000;

	Random rng = initialize_random();	

	destinations_list SLT_path = generate_destinations_box(N,l, &rng);
	std::string best_path_filename = "./out/best_path_parallel_box.dat";
	//destinations_list SLT_path = generate_destinations_circle(N,l, &rng);
	//std::string best_path_filename = "./out/best_path_parallel_circle.dat";

	rng = initialize_random(rank*10);//all threads have the same path but different random seeds for the metropolis
	
	individual_type individual(N, &SLT_path);

	double best_lenght = individual.path_lenght(), old_lenght=0;
	std::vector<double> best_lenght_vector;
	individual_type best_individual(individual);

	std::cerr << "Process "<< rank << " start its job"<<std::endl;	

		//starting with the updates in temperatures
	while(count<N_max){
		double beta = 1./T;
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
		best_lenght_vector.push_back(best_lenght);
		if (old_lenght == best_lenght)
			count++;
		else{
			count =0;
			old_lenght = best_lenght;
		}
		T= T/(c+1);
		n_iter++;
	}

	double threads_result[size];
	MPI_Gather(&best_lenght,1,MPI_DOUBLE,&(threads_result[rank]),1,MPI_DOUBLE,0,MPI::COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	int best_thread=2437;

	if(rank==0){//choosing the best path
		double best_tot_lenght=N*l*l*10;
		for(int i=0; i<size; i++){
			if(best_tot_lenght>threads_result[i]){
				best_tot_lenght = threads_result[i];
				best_thread = i;
			}
			std::cout << threads_result[i] <<std::endl;
		}
		std::cout << "Found best lenght among those of all threads... "<< std::endl << "tid: "<< best_thread << "path lenght: " << best_tot_lenght << std::endl;
	}

	MPI_Bcast(&best_thread, 1, MPI_INTEGER, 0,MPI::COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == best_thread){
		individual.print();
		individual.print_path_to_file(best_path_filename);
	}
	MPI::Finalize();

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
