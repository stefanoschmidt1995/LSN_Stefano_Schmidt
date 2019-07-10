/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 5
*************************************************************************/
#include "../standard_header.h"

//METROPOLIS ALGORITHM for sampling hydrogen orbitals

//We can sample the following probability distribution:
//	P_100(x) = N * exp(-2*r)
//	P_210(x) = N * r**2 * cos(theta)**2 * exp(-r)
//	P_310(x) = N * r*(6-r)*exp(-r/3.)*cos_theta
//	P_320(x) = pow(r,4)*pow(3.*cos_theta*cos_theta-1.,2)*exp(-2*r/3.);
//	P_210(x) =  r**2*cos_theta*cos_theta*exp(-r);
//... plus some other bonus probabilities you can easily add

typedef double (*pdf_type)(R3_point); //type for a probability distribution function

double P_100(R3_point x);
double P_200(R3_point x);
double P_310(R3_point x);
double P_320(R3_point x);
double P_210(R3_point x);
void get_orbital_parameters(std::string, pdf_type*, double*); //functions to set a suitable pdf and its RW step by specifying its name

//declaration of helper functions for Metropolis algorithm
R3_point propose_a_move(const R3_point&, double, Random*); //propoese a move for the Metropolis algorithm
bool accept_move(R3_point x, R3_point y, pdf_type pdf, Random* rng); //return whether the move should be accepted
std::vector<double> sample_from_pdf(std::string, int, Random*); //fills a file in ./out with M points sampled from pdf specified in std::string argument and returns a vector with r**2 for each sampled point

int main(){
		//************* EX 05.1 - HYDROGEN ORBITALS SAMPLING
	Random rng = initialize_random();
	int M = 200000; //number of points to be sampled
	int L = 50; //number of data within one block

	R3_point trial(3),pos(3);

	std::vector<std::string> dist_to_sample = {"1s","2p","2s","3d", "3p"};
	std::vector<double> r_vector(M);

	std::ofstream out_file;

	for (auto it = dist_to_sample.begin(); it!=dist_to_sample.end(); it++){
		std::vector<double> r_vector (sample_from_pdf(*it, M, &rng)); //doing the whole Metropolis here
		out_file.open("./out/"+(*it)+"_averages.dat");
		for(int N = 2; N<M/L; N*=1.3){ //saving averages with data blocking
			if(N<=5) N++;
			std::vector<double> avg_r = compute_statistical_error(r_vector, N, L);
			out_file << N << "\t" << avg_r[0] << "\t" <<  avg_r[2] << std::endl;
		}
		out_file.close();
	}
	return 0;
}

//FUNCTIONS USEFUL FOR THE METROPOLIS ALGORITHM
std::vector<double>
sample_from_pdf(std::string pdf_name, int M, Random* rng){
		int M_eq = 1000; //number of events for equilibration (set manually)
		int sparse_avg = 10; //sparse averaging step

		std::string filename = "./out/"+pdf_name + "_orbital_sampling.dat";
		std::cout <<"Sampling from " << pdf_name <<std::endl << "\tFilling file: " << filename <<std::endl;	
		std::fstream file(filename, std::ios::out);
		pdf_type pdf;//holds the address of pdf to be sample from
		int count_accepted =0;
		double l = 1.; //size of the step (default value)
		std::vector<double> r_vector(M);

			//initializing pos
		R3_point pos = {1,-1,2};

		get_orbital_parameters(pdf_name, &pdf, &l); //getting pdf name and a suitable step for it
		if(pdf == NULL){
			std::cerr << "Unable to find PDF for "<< pdf_name <<std::endl;
			return r_vector;
		}

			//starting equilibration & sampling
		int count =0;
		for(int i =0; i<M*sparse_avg+M_eq; i++){
			R3_point trial = propose_a_move(pos, l, rng);
			bool accepted = accept_move(pos, trial, pdf, rng);
			if(accepted){
				pos = trial;
				if(i>M_eq) count_accepted++;
			}
				//dealing with actual sampling
			if(i>M_eq && i%sparse_avg ==0){
				r_vector.at(count) = sqrt(norm(pos));
				file << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t" << std::endl;
				count++;
			}
		}
			//computing r_avg for output & printing accepatance rate
		std::vector<double> res = compute_statistical_error(r_vector, M/100);
		std::cout << "\t r_avg " <<pdf_name <<": " << res[0] << " " << res[2] << std::endl;
		std::cout << "\t Moves accepted " <<count_accepted/(double)(M*sparse_avg) << std::endl;
		file.close();
		return r_vector;
}




R3_point
propose_a_move(const R3_point& x, double l, Random* rng){
	//Returns a vector holding the position to step to given the "current  position" x
	if(x.size() != 3) return {};
	R3_point res(x);
	for(int i=0; i<=2; i++)
			res[i] = x[i] + rng->Rannyu(-l,l); //uniform transition probability
			//res[i] = x[i] + rng->Gauss(0,l/1.7); //gaussian transition probability
	return res;
}

bool
accept_move(R3_point x, R3_point y, pdf_type pdf, Random* rng){
	//given the proposed move from x to y it returns true if the move is to be accepted according to Metropolis criterium
	double w = pdf(y)/pdf(x);
	double r = rng->Rannyu();
	if(r<=w)
		return true;
	else
		return false;
}

//ORBITALS PROBABILITY DISTRIBUTION FUNCTIONS

void
get_orbital_parameters(std::string orbital_name, pdf_type* pdf_name, double* l){
		//defining orbitals pdfs and l parameter
	if (orbital_name == "1s"){
		*pdf_name = P_100;
		*l=1.2;
	}
	else if (orbital_name == "2p"){
		*pdf_name = P_210;
		*l = 3.;
	}
	else if (orbital_name == "2s"){
		*pdf_name = P_200;
		*l = 4.6;
	}
	else if (orbital_name == "3d"){
		*pdf_name = P_320;
		*l = 5.;
	}
	else if (orbital_name == "3p"){
		*pdf_name = P_310;
		*l = 6.1;
	}
	else
		*pdf_name = NULL;

	return;
}

double
P_100(R3_point x){
	double r = sqrt(norm(x));
	return exp(-2.*r);
}

double
P_200(R3_point x){
	double r = sqrt(norm(x));
	return (2-r)*(2-r)*exp(-r);
}

double
P_310(R3_point x){
	double r = sqrt(norm(x));
	double cos_theta = x.back()/r;
	double psi = r*(6-r)*exp(-r/3.)*cos_theta;
	return psi*psi;
}


double
P_320(R3_point x){
	double r = sqrt(norm(x));
	double cos_theta = x.back()/r;
	return pow(r,4)*pow(3.*cos_theta*cos_theta-1.,2)*exp(-2*r/3.);
}

double
P_210(R3_point x){
	double r = sqrt(norm(x));
	double cos_theta = x.back()/r;
	return r*r*cos_theta*cos_theta*exp(-r);
}



