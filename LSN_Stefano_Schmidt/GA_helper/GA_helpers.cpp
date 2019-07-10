/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
Genetic Algorithm HELPERS - DECLARATIONS
**************************************************************************
Definitions of functions declared in GA_helpers.h
	methods of individual_type class
	functions to generate a list of cities to be visited by the salesman
	function useful for generating children from parents and other helpers
*************************************************************************/
#include "GA_helpers.h"
	//methods of individual_type class
double
individual_type::path_lenght(){
	if(destinations == NULL) return -1.;
	double lenght=0;
	for(auto it = this->begin(); it!=this->end()-1; it++){
		lenght += norm(destinations->at(*(it+1)) - destinations->at(*(it)));
	}
	return lenght;
}

double
individual_type::fitness(){
	return pow(offset - path_lenght(),2);
}

void
individual_type::random_shuffle(){
	this->at(N) = 0;
	for(unsigned int i=0; i<this->size()-1; i++)
		this->at(i) = i;
	std::random_shuffle(this->begin()+1,this->end()-1);
	return;
}


void 
individual_type::perform_mutation(Random* rng){
	//different kinds of mutations are possible...
	//One of them is selected at random
	int N_mutations = 3; //# possible mutations
	int which_mutation = ((int) (1000*rng->Rannyu())) % (N_mutations);
	int a,b,temp;
	

	switch(which_mutation){
		case 0: //random swap
			a = ((int) (10000*rng->Rannyu())) % (N-1)+1;
			b=0;
			do{
				b = ((int) (10000*rng->Rannyu())) % (N-1)+1;
			}while(a==b);

			temp = this->at(a);
			this->at(a) = this->at(b);
			this->at(b) = temp;
			break;
		case 1: //random swap between neighbours
			a = ((int) (10000*rng->Rannyu())) % (N-1)+1;
			b=a+1;
			if( b>=N) b=1;
			temp = this->at(a);
			this->at(a) = this->at(b);
			this->at(b) = temp;
			break;
		case 2: //permutation beetween group of cities
			int m_step = ((int) (10000*rng->Rannyu())) % (N/2);
			int start_a = ((int) (10000*rng->Rannyu())) % (N/2-m_step) + 1; //index to start the permutation
			//std::cout << start_a << "\t"<< m_step <<std::endl;
			//print(); //DEBUG
			int start_b = start_a + N/2;
			for(int i =0; i< m_step; i++){
				temp = this->at(start_a+i);
				this->at(start_a+i) = this->at(start_b+i);
				this->at(start_b+i) = temp;
			}
			//print(); //DEBUG
			break;
	}

	return;
}

bool
individual_type::is_all_right(){
	//check whether everything in the sequence of city is fine
	if (this->at(0) != 0 ) return false;
	if (this->at(this->size()-1) != 0 ) return false;
	int count =0;
	for (unsigned int i =0; i< this->size()-2; i++){
		if (std::find(this->begin(), this->end(), i)!=this->end()) //looking for i-th city
			count++;
	}
	if (count != (int)this->size()-2) {return false;}
	return true;
}


void
individual_type::print(){
	for(auto it = this->begin(); it!=this->end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	return;
}

void
individual_type::print_path_to_file(std::string filename){
	std::ofstream file(filename);
	destinations_list* destinations = this->get_destinations_list();
	int n_dim = (destinations->at(0)).size();
	file << "lenght: " << this->path_lenght() <<std::endl;
	for(auto it= this->begin(); it!= this->end(); it++){
		for(int j=0; j<n_dim; j++)
			file << (destinations->at(*it))[j] << "\t";
		file << std::endl;
	}
	file.close();
	return;
}

	//some useful functions to generate a list of cities to be visited by the salesman
//generating configuration for the cities...
destinations_list
generate_destinations_box(int N, double l, Random* rng, int dim){
	destinations_list res(N);
	for(int i =0; i<N ; i++){
		R3_point temp_point(dim);
		for(int j=0; j<dim; j++)
			temp_point[j] = rng->Rannyu(-l,l);
		res[i]= temp_point;
	}
	return res;
}

destinations_list
generate_destinations_circle(int N, double l, Random* rng, int dim){
	destinations_list res(N);
	for(int i =0; i<N ; i++){
		R3_point temp_point(dim);
		do{
			for(int j=0; j<dim; j++)
				temp_point[j] = rng->Rannyu(-l,l);	
		}while(norm(temp_point)>l*l);
		res[i]= temp_point;
	}
	return res;
}

	//function useful for generating children from parents and other helpers
std::vector<double>
get_fitness_rank(population_type* pop){
	std::vector<double> fitness_rank(pop->size());
	for(unsigned int jj = 0; jj < pop->size(); jj++)
		fitness_rank[jj] = pop->at(jj).fitness();
	return fitness_rank;
}

individual_type
create_son(individual_type* a, individual_type* b, Random* rng){
	if (a->size()!=b->size()) return *a;
	if (a->fitness()<b->fitness()){
		//swapping a and b so that the individual with more fitness is at first
		individual_type* temp = a;
		a = b;
		b = temp;
	}

	individual_type res(*a);

		//replacing the last half of a with the last half of b s.t. there are no superposition
	int j =res.size()-2;
	int random_separation = (int)(10000*rng->Rannyu()) % (res.size()/4) + (int) (res.size()/4) +1;
	for (auto it = res.end()-2 ; it!= res.begin()+ random_separation; it--)
		*it =0;

	for (auto it = res.end()-2 ; it!= res.begin()+ random_separation; it--){
		for(unsigned int jj=j; jj>0; jj--){
			if(std::find(res.begin(), res.end(), b->at(jj)) == res.end() ){
				//found a new value for *it: replace it
				*it = b->at(jj);
				j=jj-1;
				break;
			}
		}
	}
	if (res.is_all_right())
		return res;
	else
		return *a;
}

individual_type
create_son_weak(individual_type* a, individual_type* b, Random* rng){
	if (a->size()!=b->size()) return *a;
	if (a->fitness()<b->fitness()){
		//swapping a and b so that the individual with more fitness is at first
		individual_type* temp = a;
		a = b;
		b = temp;
	}

	individual_type res(*a);

	//select a sequence of a
	unsigned int start = (int)(10000*rng->Rannyu()) % (res.size()/2);
	unsigned int end = start + (int)(10000*rng->Rannyu()) % (res.size()/2);
	
	for (unsigned int i =0; i<res.size(); i++){
		if((i<start) || (i>end))
			res[i] = 0;
	}
	int j =1;
	for (auto it = res.begin()+1 ; it!= res.end()-1; it++){
		unsigned int index = (int)std::distance(res.begin(), it);
		if((index>=start) && (index<=end))
			continue;
		for(unsigned int jj=j; jj<res.size(); jj++){
			if(std::find(res.begin(), res.end(), b->at(jj)) == res.end() ){
				//found a new value for *it: replace it
				*it = b->at(jj);
				j=jj+1;
				break;
			}
		}
	}
	if (res.is_all_right())
		return res;
	else
		return *a;

}

double get_best_half_average_lenght(population_type* pop){
	std::vector<double> fitness_rank = get_fitness_rank(pop);
	double dist = 0;
	for(unsigned int jj = 0; jj < pop->size()/2; jj++){
		int fitter_index = std::distance(fitness_rank.begin(), std::max_element(fitness_rank.begin(),fitness_rank.end()));
		dist+= pop->at(fitter_index).path_lenght();
		fitness_rank.at(fitter_index) = 0;
	}
	return dist/(double)(pop->size()/2);
}





