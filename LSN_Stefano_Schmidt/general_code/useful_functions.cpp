/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
USEFUL FUNCTIONS
**************************************************************************
Definition of some functions useful for many of the exercises in the course...
*************************************************************************/
#include "../standard_header.h"

Random
initialize_random(){
	Random res = initialize_random(0);
	return res;
}

Random
initialize_random(int rand_seed){
	//Function used to initialize the random generator according to numbers hold in files Primes and seed.in
	Random rnd;
	int seed[4];
	int p1, p2;
	std::ifstream Primes("../general_code/Primes");
	if (Primes.good()){
	  Primes >> p1 >> p2 ;
	} else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
	Primes.close();

	std::ifstream input("../general_code/seed.in");
	std::string property;
	if (input.good()){
	  while ( !input.eof() ){
		 input >> property;
		 if( property == "RANDOMSEED" ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1+rand_seed,p2-rand_seed);
		 }
	  }
	  input.close();
	} else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
	return rnd;
}

double
norm(point x){
	int res =0;
	for(auto it = x.begin(); it!=x.end(); it++)
		res+=(*it)*(*it);
	return (double)res;
}

double
norm(R3_point x){
	double res =0;
	for(auto it = x.begin(); it!=x.end(); it++)
		res+=(*it)*(*it);
	return res;
}

R3_point operator+(const R3_point& a, const R3_point& b) {
	if (a.size() != b.size()) return {};
	R3_point res(a.size());
	for(unsigned int i =0; i< a.size(); i++)
		res[i] = a[i] + b[i];
	return res;
}

R3_point operator-(const R3_point& a, const R3_point& b) {
	if (a.size() != b.size()) return {};
	R3_point res(a.size());
	for(unsigned int i =0; i< a.size(); i++)
		res[i] = a[i] - b[i];
	return res;
}

std::vector<double> compute_statistical_error(const std::vector<double>& values, int N_blocks){
	return compute_statistical_error(values, N_blocks, values.size()/N_blocks);
}

std::vector<double>
compute_statistical_error(const std::vector<double>& values, int N_blocks, int L){
	//Function that uses blocking method with N_blocks to compute AVERAGE, VARIANCE and STANDARD DEVIATION OF MEAN of a set of values.
	//L is the number of element to be put in a block
	if (L<1)
		L = values.size()/N_blocks; //if it's a bad (negative) value L is chosen s.t. all values in values vector are used
	if (values.size() <  L*N_blocks){
		std::cerr << "There are too few data to perform the blocking with " << N_blocks << " blocks and " << L << " data per block." << std::endl;
		return std::vector<double>(3);
	}
	
	std::vector<double> avg(N_blocks);
	for (int i = 0; i<N_blocks; i++){
		for (int j =0; j<L; j++)
			avg[i] += values.at(L*i+j);
		avg[i] = avg[i]/(double)L;
		//std::cout<<"avg: " <<avg[i] << std::endl; //DEBUG
	}
	double var = 0, block_average =0;
	for (int i = 0; i<N_blocks; i++){ //computing variance and average of the blocks
		block_average += avg[i];
		var += avg[i]*avg[i];
	}
	std::vector<double> result(3);
	result[0]= block_average/N_blocks;
	result[1]= var/N_blocks - result[0]*result[0];
	result[2]= sqrt(result[1] / (N_blocks-1));
	if (result[1] < 1e-8){
		result[1]=0;
		result[2]=0;
	}
	return result;
}
