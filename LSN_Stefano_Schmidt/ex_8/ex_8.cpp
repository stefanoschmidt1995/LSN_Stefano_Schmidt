/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 8
*************************************************************************/
#include "../standard_header.h"

//*****DECLARATION OF SOME USEFUL FUNCTIONS
typedef double (*pdf_type)(double, double, double); //type for a probability distribution function 1D with 2 parameters
typedef double (*potential_type)(double); //type for a potential function

double propose_a_move(const double& x, double l, Random* rng);
	//propose a metropolis step with uniform probability
double standard_pdf(double x, double mu, double sigma);
	//holds the standard double exponential pdf to optimize
double standard_potential(double x); //hold the potential in which th particle is moving
std::vector<double> perform_variational_MC(double mu, double sigma, int M, Random* rng, std::string to_sample, pdf_type pdf = standard_pdf, potential_type pot = standard_potential, double mass = 1.);
	//samples from the wavefunction specified in pdf if to_sample = "x"
	//samples from energy configuration of the wavefunction specified in pdf moving in the potential pot
bool accept_move(double x, double y, double mu, double sigma, Random* rng, pdf_type pdf = standard_pdf);
	//returns wheter to accept a metropolis move
void choose_parameters(double* mu, double* sigma, Random* rng);
	//optimises (with a sort of Tabu search) the couple of mu and sigma choosing those who minimises average energy as computed by metropolis algorithm
double get_energy(double x, double mu, double sigma, potential_type pot = standard_potential, double mass = 1.);
	//contribution to energy given by the the point x
double get_avg_energy(double mu, double sigma, int N_iter, Random* rng);
	//returns the energy averaged on N_iter of metropolis samplings
std::list<std::vector<double>> get_neighbours_list(std::vector<double> center, int N_neig, std::vector<double> l_step, Random* rng);
	//return a list of randomly chosen N_neig neighbours to center in parameter space of mu, sigma.

//*****MAIN FUNCTION
int main(){
		//************* EX 08.1 - METROPOLIS ALGORITHM & GROUND STATE OF A 1D PARTICLE
	Random rng = initialize_random();

	int M = 500000; //number of events to be sampled
	double mu = 1., sigma = .1; //initial guess of parameters

		//choosing best parameters (i.e. those who minimize the energy)
	choose_parameters(&mu, &sigma, &rng);
	std::cout <<"Chosen parameter: " << mu << " " << sigma << std::endl;

		//sampling psi**2 
	std::cout << "Sampling "<< M <<" positions according to wavefunction"<<std::endl;
	std::vector<double> sampled_positions = perform_variational_MC(mu, sigma, M, &rng, "x");
	std::ofstream out_file("positions.dat");
	out_file <<mu <<"\t"<<sigma<<std::endl;
	for (auto it= sampled_positions.begin() ; it!= sampled_positions.end(); it++) out_file << *it << "\n";
	out_file.close();

		//writing (as usual) energy for fixed L=50
	const int L =50;
	std::cout << "Sampling "<< M <<" values of energy according to wavefunction in the given potential"<<std::endl;
	std::vector<double> sampled_energies = perform_variational_MC(mu, sigma, M, &rng, "E");
	out_file.open("energies.dat");
	for(int N = 2; N<M/L; N*=1.3){
		if(N<=5) N++;
		std::vector<double> res = compute_statistical_error(sampled_energies,N, L);
		out_file << N << "\t" << res[0] << "\t" << res[2] << std::endl;
	}
	out_file.close();

	return 0;
}

//*****DEFINITION OF SOME USEFUL FUNCTIONS (declared above)
double
propose_a_move(const double& x, double l, Random* rng){
	double res=0;
	res = x + rng->Rannyu(-l,l);
	return res;
}

double
standard_pdf(double x, double mu, double sigma){
	double res = exp(- pow(x-mu,2)/(2*sigma)) + exp(- pow(x+mu,2)/(2*sigma));
//	double norm_factor = 2*(1+exp(-mu*mu/sigma))*sqrt(M_PI*sigma); //useless: it is simplified in the Metropolis algorithm		
//	return res*res/norm_factor;
	return res*res;
}

double
standard_potential(double x){
	return pow(x,4) - 5.*x*x/2.;
}


double
get_energy(double x, double mu, double sigma, potential_type pot, double mass){
	double potential_energy = pot(x);
	double kin_energy = -( pow(x-mu,2) * exp(-pow((x-mu),2)/(2*sigma) ) + pow(x+mu,2) * exp(-pow((x+mu),2)/(2*sigma)) )/ (2.*mass*sigma*sigma);
	kin_energy = kin_energy + ( exp(-pow((x-mu),2)/(2*sigma)) + exp(-pow((x+mu),2)/(2*sigma)) ) / (2.*mass*sigma);
	kin_energy /= ( exp(- pow(x-mu,2)/(2*sigma)) + exp(- pow(x+mu,2)/(2*sigma)) );
	return potential_energy + kin_energy;
}


bool
accept_move(double x, double y, double mu, double sigma, Random* rng, pdf_type pdf){
	//given the proposed move from x to y it returns true if the move is to be accepted according to Metropolis criterium
	double w = pdf(y, mu, sigma)/pdf(x, mu, sigma);
	double r = rng->Rannyu();
	if(r<=w)
		return true;
	else
		return false;
}

std::vector<double>
perform_variational_MC(double mu, double sigma, int M, Random* rng, std::string to_sample, pdf_type pdf, potential_type pot, double mass){
	//This function sample M points from the pdf given (which accepts two parameters mu and sigma)
	//The value are sampled from n_loops markov chains with a sparse averaging step of 10

	if(to_sample != "E" && to_sample!= "x")	to_sample= "E"; 

	int n_loops = 200; //# of new initialization done at every time
	int M_eq = 100; //number of events for equilibration (set manually)
	int sparse_avg_const = 10;
	std::vector<double> r_vector(M); //result vector
	int index =0; //keep track of the last position in r_vector to write in

	int count_accepted =0;
	double extra_factor = .5;
	if (sigma<.1) extra_factor = (-log(sigma));
	double l = sigma*10*extra_factor; //size of the step (default value)

		//starting equilibration & sampling
	for(int loop = 0; loop < n_loops ; loop ++){
		int signum = (2*((int)(rng->Rannyu()*100)%2)) -1;
		double trial,pos = rng->Rannyu(-10*sigma,10*sigma)+(double)signum*mu; //initializing
		for(int i =0; i<M_eq + M*sparse_avg_const / n_loops; i++){
			trial = propose_a_move(pos, l, rng);
			bool accepted = accept_move(pos, trial, mu, sigma, rng, pdf);
				
			if(accepted){
				pos = trial;
				if(i>M_eq) count_accepted++;
			}
				//dealing with actual sampling
			if(i>M_eq && i%sparse_avg_const == 0){ //saving the current energy or position
				if(to_sample == "E")
					r_vector.at(index) = get_energy(pos, mu, sigma, pot, mass); //pay attention to indices!!!!!
				if(to_sample == "x")
					r_vector.at(index) = pos;
				index++;
			}
		}
	}
	//std::cout << "\t Moves accepted " <<count_accepted/(double)(M*sparse_avg_const) << std::endl;
	return r_vector;
}

std::list<std::vector<double>>
get_neighbours_list(std::vector<double> center, int N_neig, std::vector<double> l_step, Random* rng){
	std::list<std::vector<double>> res;
	for(int i=0; i <N_neig; i++){
		std::vector<double> new_s(center);
		for (unsigned int j =0; j< center.size(); j++){
			double step = rng->Rannyu(-l_step[j], l_step[j]);
			new_s[j] += step;
			if (new_s[j]<0) new_s[j] = -(new_s[j]);
		}
		res.push_back(new_s);
	}
	return res;
}

double
get_avg_energy(double mu, double sigma, int N_iter, Random* rng){
	if(N_iter<100) N_iter =100;
	std::vector<double> res = perform_variational_MC(mu,sigma, N_iter, rng, "E");
	res =  compute_statistical_error(res, N_iter/50);
	return res[0];
}

void
choose_parameters(double* mu, double* sigma, Random* rng){
	//this function set mu and sigma s.t. the energy of the wavefunction is minimized.
		//To do so it uses something similar to Tabu search https://en.wikipedia.org/wiki/Tabu_search

	int N_neig =10; //neighbours to be considered
	int N_average= 1000; //metropolis step to sample in order to get acg energy
	int n_useless_iter = 0; //criterion to exit the loop
	int max_iter = 20;
	std::vector<double> l_step = {.1,.05};


	std::vector<double> sBest = {*mu, *sigma};
	double best_energy = get_avg_energy(sBest[0], sBest[1], N_average, rng);
	double old_best_energy = best_energy;

	std::cout << "Choosing parameters which minimize energy" << std::endl;

	while(n_useless_iter<max_iter){
		//searching in the neighbour of the best candidate sBest
		double current_energy =0;
		std::list<std::vector<double>> neighbours_list = get_neighbours_list(sBest, N_neig ,l_step, rng);
		for(auto it = neighbours_list.begin(); it != neighbours_list.end(); it++){
			current_energy = get_avg_energy((*it).at(0), (*it).at(1), N_average, rng);
			if (current_energy<best_energy){
				sBest = *it;
				best_energy = current_energy;
			}
		}
		std::cout <<"\tCurrent best energy: " <<best_energy << std::endl;
		if (best_energy<old_best_energy){
			old_best_energy = best_energy;
			n_useless_iter =0;
		}
		else
			n_useless_iter++;
	}
	*mu= sBest[0];
	*sigma= sBest[1];
	return;
}


