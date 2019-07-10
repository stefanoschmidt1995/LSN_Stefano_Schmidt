/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 9
*************************************************************************/
#include "../standard_header.h"
#include "../GA_helper/GA_helpers.h"

int
main (int argc, char *argv[]){

		//************* EX 09.1 - GENETIC ALGORITHM for solving Traveling Salesman Problem
	Random rng = initialize_random();

	int N_gen = 5001; //N_gen of generations ot be created
	int M = 200, N = 30; //M population number, N path lenght
	double l =1.; //size of the box
	double fitness_offset = 2*(double)N*l*l; //RW scaling
	double lambda = 1.5*M; //gives the avg number of mutations in a generation

		//N_gen, M, N can be passed to main
	if(argc >=3){
		N_gen = atof(argv[1]);
		M = atof(argv[2]);
		if(argc>4) N = atof(argv[3]);
	}
	if (argc>5 || argc == 2 ){
		std::cerr << "Usage: ./es_9.cpp N_gen M_pop (N_path)" << std::endl;
		return 1;
	}

	std::ofstream distance_file("distances.dat");

	std::default_random_engine poisson_generator;
	std::poisson_distribution<int> poisson_distribution(lambda); 	

	population_type* population = new population_type(M);
	destinations_list SLT_path = generate_destinations_box(N, l, &rng);
//	destinations_list SLT_path = generate_destinations_circle(N, l, &rng);

		//STARTING GENETIC ALGORITHM
			//STEP 1 - Initializing population at random
	for(int i =0; i<M; i++){
		individual_type temp_gene(N, &SLT_path, fitness_offset);
		population->at(i)= temp_gene;
	}

	individual_type best_individual(population->at(0)); //it's not important who is the actual best individual at time 0

	for(int n_gen=0; n_gen<N_gen; n_gen++){
		std::vector<double> fitness_rank(M);
		if (n_gen%25==0){
			std::cout << "Generation # " << n_gen <<std::endl;
			best_individual.print_path_to_file("./paths/"+std::to_string(n_gen)+".dat");
		}
			//STEP 2 - Do some mutations
		int how_many_mutations =  poisson_distribution(poisson_generator);
		for(int j=0; j< how_many_mutations; j++){
			//choose here a individual that shall undergo a mutation
			int to_change = ((int) (10000*rng.Rannyu())) % (M);
			(population->at(to_change)).perform_mutation(&rng);
		}

			//STEP 3 - Creating next generations
		population_type* new_population = new population_type(M);

		fitness_rank = get_fitness_rank(population);
		int fitter_index = std::distance(fitness_rank.begin(), std::max_element(fitness_rank.begin(),fitness_rank.end()));
		if((population->at(fitter_index)).path_lenght() < best_individual.path_lenght()) //saving best individual in the new population
			best_individual = population->at(fitter_index);

		new_population->at(0) = best_individual; //first of the population is always the best individual

		for(int j=1; j<M; j++){
				//choosing two parents for the child
			int a= rng.sample_from_discrete_p(&fitness_rank);
			int b= rng.sample_from_discrete_p(&fitness_rank); //this means that can happen partenogenesis
				//creating a child
			new_population->at(j) = create_son(&population->at(a),&population->at(b), &rng);
		}
			//STEP 4 - New generation takes control
		delete population;
		population = new_population;

			//computing best individual and saving to file
		fitness_rank = get_fitness_rank(population);
		fitter_index = std::distance(fitness_rank.begin(), std::max_element(fitness_rank.begin(),fitness_rank.end()));
		if((population->at(fitter_index)).fitness() > best_individual.fitness()) //saving best individual in the new population
			best_individual = population->at(fitter_index);
		distance_file << n_gen << "\t" << (population->at(fitter_index)).path_lenght() << "\t" << get_best_half_average_lenght(population)<<std::endl;

	}
	best_individual.print_path_to_file("best_path.dat");
	delete population;
	distance_file.close();

	return 0;
}
