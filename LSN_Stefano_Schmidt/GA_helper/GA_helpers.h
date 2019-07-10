/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
Genetic Algorithm HELPERS - DECLARATIONS
**************************************************************************
File holding all the relevant declaration of helper functions and classes for solving STP problem with genetic algorithm
	methods of individual_type class
	functions to generate a list of cities to be visited by the salesman
	function useful for generating children from parents and other helpers
*************************************************************************/

#include "../standard_header.h"

typedef std::vector<R3_point> destinations_list;

destinations_list generate_destinations_box(int N, double l, Random* rng, int dim=2); 
destinations_list generate_destinations_circle(int N, double l, Random* rng, int dim=2);

class
individual_type: public std::vector<int> {
	//gene is a sequence of indices in [0,N-1] starting and ending with 0 and having a permutation of 1 - N-1 in the middle
	//it represents the sequence of places to be visited in an ORDERED list of destination (destinations_list)
	public:
		individual_type(std::vector<int> sequence, destinations_list* destinations=NULL, double fitness_offset=0.):
			vector(sequence),N(sequence.size()), destinations(destinations),offset(fitness_offset)
			{random_shuffle(); if(fitness_offset<0) fitness_offset =0.;};
		individual_type(int N=10, destinations_list* destinations=NULL, double fitness_offset=0.):
			vector(N+1),N(N), destinations(destinations),offset(fitness_offset)
			{random_shuffle(); if(fitness_offset<0) fitness_offset =0.;};
		~individual_type(){};
		double fitness(); //pow(offset - path_lenght(),2);
		double path_lenght();
		void random_shuffle();
		void perform_mutation(Random* rng);
		void print();
		void print_path_to_file(std::string filename);
		destinations_list* get_destinations_list(){return destinations;};
		double get_offset(){return offset;};
		bool is_all_right();
	private:
		int N;
		destinations_list* destinations;
		double offset;
};

individual_type create_son(individual_type*, individual_type*, Random*); //creates a son from the two parents individual
individual_type create_son_weak(individual_type* , individual_type*, Random*); //creates a weak son from the two parents individual
typedef std::vector<individual_type> population_type;
double get_best_half_average_lenght(population_type*); //average over the path lenght of the first best half
std::vector<double> get_fitness_rank(population_type*); //return the fitness of every individual of population


