/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
CLASS ISING SIMULATOR - DEFINITIONS
**************************************************************************
Definitions of class ising simulator declared in ising_simulator.h
*************************************************************************/
#include "ising_simulator.h"

ising_simulator::ising_simulator(std::string filename){
	std::cout << "Classic 1D Ising model             " << std::endl;
	std::cout << "Monte Carlo simulation             " << std::endl << std::endl;
	std::cout << "Nearest neighbour interaction      " << std::endl << std::endl;
	std::cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << std::endl << std::endl;
	std::cout << "The program uses k_B=1 and mu_B=1 units " << std::endl;

	config_from_file(filename);

	//Prepare arrays for measurements
	measure_vector = std::vector<std::list<double>>(n_props);
	for(auto it = measure_vector.begin(); it!= measure_vector.end(); it++){
		std::list<double> temp_list;
		*it = temp_list;		
	}
		

	s = std::vector<int>(N);
	//initial configuration at random
	for (int i=0; i<N; ++i)
	{
		if(rng.Rannyu() >= 0.5) s[i] = 1;
		else s[i] = -1;
	}
	
//Evaluate energy etc. of the initial configuration
	measure();

	std::cout << "Initial energy per spin = " << *(measure_vector[0].begin())/(double)N << std::endl;
	return;
}

void
ising_simulator::config_from_file(std::string filename){
	std::ifstream ReadInput;

	//setting random seed
	rng = initialize_random();

	//Read input informations
	ReadInput.open(filename);

	ReadInput >> temp;
	if (temp<=0)
		beta = 1e10;
	else
		beta = 1.0/temp;
	std::cout << "Temperature = " << temp << std::endl;

	ReadInput >> N;
	std::cout << "Number of spins = " << N << std::endl;

	ReadInput >> J;
	std::cout << "Exchange interaction = " << J << std::endl;

	ReadInput >> h;
	std::cout << "External field = " << h << std::endl << std::endl;
		
	ReadInput >> metropolis; // if=1 Metropolis else Gibbs

	ReadInput >> n_meas;

	if(metropolis) std::cout << "The program perform Metropolis moves" << std::endl;
	else std::cout << "The program performs Gibbs moves" << std::endl;
	std::cout << "Number of measures for each block = " << n_meas << std::endl << std::endl;
	ReadInput.close();
	
	return;
}

void
ising_simulator::set_spin_configuration(std::string filename){
	std::ifstream ReadInput;
	ReadInput.open(filename);	

	int count =0;
	while(ReadInput.good()){
		ReadInput >> s[count]; 
		count++;
		if(count >= N) break;
	}

	ReadInput.close();
	return;
}


int
ising_simulator::Pbc(int i){  //Algorithm for periodic boundary conditions
	if(i >= N) i = i - N;
	else if(i < 0) i = i + N;
	return i;
}

void
ising_simulator::measure(){
	std::vector<double> u(n_props);
	for (int i=0; i<N; ++i){
		u[0] += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]); //energy
		u[2] += s[i]; //mag
	}
	u[1] = u[0]*u[0]; //heat
	u[3] = u[2]*u[2]; //susceptivity
	for(int i=0; i< n_props; i++)
		measure_vector[i].push_back(u[i]);
	return;
}

std::vector<avg_type>
ising_simulator::get_measures(){
	std::vector<avg_type> result(n_props);
	int n_blocks= (int)measure_vector.at(0).size()/n_meas;
	for (int i = 0; i<n_props; i++){
		std::vector<double> temp_vector{ std::begin(measure_vector[i]), std::end(measure_vector[i]) };
		result[i] = compute_statistical_error(temp_vector, n_blocks ,n_meas);
	}
		//dealing with specific heat (index = 1)
	//std::cerr << result[1][0] << "\t"<<result[1][2] << "\t" << result[0][0] << "\t"<<result[0][2] <<std::endl;//DEBUG
	result[1][0] = pow(beta,2) * (result[1][0] - result[0][0]*result[0][0]);
	result[1][1] = pow(beta,4) * (result[1][1] + 2* result[0][1]); //error propagation (assuming null covariance)
	result[1][2] = pow(beta,2) * sqrt(pow(result[1][2],2) + 2*pow(result[0][2],2));
	return result;
}

void 
ising_simulator::move(){ //tries to flip N spins at random
	int o; //holds the index of the spin to be moved
	double acceptance;

	for(int i=0; i<N; ++i){
		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		o = (int)(rng.Rannyu()*N); //index of spin to flip

		if(metropolis) //Metropolis
		{
			acceptance = exp(-beta*get_energy_difference(o));
			//acceptance =0;
			if (rng.Rannyu()<=acceptance){
				s[o] = - s[o];
				accepted++;
			}
			moves++;
		}
		else //Gibbs sampling
		{
			s[o] = draw_from_conditional(o);
		}
	}
	return;
}

double
ising_simulator::get_energy_difference(int spin_index) { //returns the energy difference new - old for the flip of spin with index spin_index
	return 2.*J * s[spin_index] * ( s[Pbc(spin_index-1)] + s[Pbc(spin_index+1)] ) + 2. * h * s[spin_index];
}

int
ising_simulator::draw_from_conditional(int spin_index){
	double sum_factor = s[Pbc(spin_index-1)] + s[Pbc(spin_index+1)] * J + h;
	double p_up = 1./(1. + exp(-2*beta* sum_factor) );
	double p_down = 1./(1. + exp(2*beta*J * sum_factor) );
	p_up = p_up / (p_up+p_down);
	if(rng.Rannyu()<=p_up)
		return 1;
	else
		return -1;
}


void
ising_simulator::reset_measures(){
	for(auto it = measure_vector.begin(); it!= measure_vector.end(); it++)
		(*it).erase((*it).begin(), (*it).end());
	return;
}


void
ising_simulator::write_to_file(std::string filename){
	std::ofstream write_conf;

	std::cout << "Print configuration to file " << filename << std::endl << std::endl;
	write_conf.open(filename);
	for (int i=0; i<N; ++i)
		write_conf << s[i] << std::endl;
	write_conf.close();
	return;
}

std::vector<double>
ising_simulator::get_config_parameters(){
	std::vector<double> res(6);
	res[0]= temp;
	res[1] = N;
	res[2] = J;
	res[3] = h;
	res[4] = (double) metropolis;
	res[5] = n_meas;
	return res;
}

void
ising_simulator::set_temperature(double new_temperature){
	temp = new_temperature;
	if (temp<=0)
		beta = 1e10;
	else
		beta = 1.0/temp;
	for (int i=0; i<N; ++i)
	{
		if(rng.Rannyu() >= 0.5) s[i] = 1;
		else s[i] = -1;
	}
	return;
}

