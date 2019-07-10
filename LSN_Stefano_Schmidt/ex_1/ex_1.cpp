/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 1
*************************************************************************/
#include"../standard_header.h"

double
get_uniform_sintheta(Random rng){
	double x=1,y=1;
	while(x*x+y*y > 1){	
		x = rng.Rannyu(-1,1);
		y = rng.Rannyu();
	}
	return y/sqrt(x*x+y*y);

}

int main(){
	Random rng = initialize_random();

		//************* EX 01.1
			//*****1) COMPUTING AN INTEGRAL
	int M = 1000000; //# of total throws
	int L = 50; //L holds the number of data within one block
	int N_max = M/L;  //maximum number of blocks

	std::vector<double> random_numbers(M);
	for (auto it = random_numbers.begin(); it!= random_numbers.end(); it++) //vector of M random number filled
		(*it) = rng.Rannyu(0.,1.);

	std::ofstream avg_file("integral.dat");
	for(int N = 2; N<N_max; N*=1.3){
		if(N<=5) N++;
		avg_type res(compute_statistical_error(random_numbers, N, L));
		avg_file << N << "\t" << res[0] << "\t" << res[2] << std::endl;
	}
	avg_file.close();

			//*****2) COMPUTING ANOTHER INTEGRAL
	std::ofstream var_file("variances.dat");
	for (auto it = random_numbers.begin(); it!= random_numbers.end(); it++) //filling the vector to average
		(*it)= pow((*it) - .5, 2);

	for(int N = 2; N<N_max; N*=1.3){
		if(N<=5) N++;
		avg_type res(compute_statistical_error(random_numbers, N,L));
		var_file << N << "\t" << res[0] << "\t" << res[2] << std::endl;
	}
	var_file.close();

			//*****3) CHI-SQUARE TEST
	int M_sub = 1e2; //# of subinterval to divide [0,1]
	int n = 1e7; //N_throws
	double chi_squared =0;
	double np = n/(double)M_sub; //avg points in a bin (inverse of the bin size)
	std::ofstream chi_file("chi.dat");
	std::vector<int> chi_vector(M_sub);

	for(int i =0; i<n; i++){ //filling histogram
		double num = rng.Rannyu(0.,1.);
		int bin =0;
		while((double)(bin+1)/M_sub <  num)
			bin++;
		chi_vector.at(bin)++ ;
	}
		//computing chi**2 for every bin
	for (auto it = chi_vector.begin(); it!= chi_vector.end(); it++){
		double chi_i = pow((double)(*it) - np, 2) / np;
		chi_squared += chi_i;
		chi_file << it - chi_vector.begin() <<"\t" << chi_i<< std::endl;
		//chi_file << it - chi_vector.begin() <<"\t" << *it<< std::endl;
	}

	chi_file.close();
	std::cout << "Chi value: " << chi_squared <<std::endl;


		//************* EX 01.2
	int N_throws = 1e4;
	std::vector<int> N = {1,2,10,100};

	std::ofstream standard_sample_file("standard.dat");
	std::ofstream exp_sample_file("exp.dat");
	std::ofstream C_L_sample_file("C_L.dat");

	//for each file:
		//colums different S_N
		//rows different realizations
	for(int i=0; i<N_throws; i++){
		for (auto it = N.begin(); it !=N.end();it++){
			double stand_sum=0., exp_sum =0., CL_sum =0.;
			for(int n = 0; n<*it; n++){ //summing over N = *it realization of the variables
				stand_sum += rng.Rannyu();
				exp_sum += rng.Exponential(1);
				CL_sum += rng.Cauchy_Lorentz(0,1.);
			}
			standard_sample_file << stand_sum/(*it)  << "\t";
			exp_sample_file << exp_sum/(*it) << "\t";
			C_L_sample_file << CL_sum/(*it)  << "\t";
		}
		standard_sample_file << std::endl;
		exp_sample_file << std::endl;
		C_L_sample_file << std::endl;
	}

	standard_sample_file.close();
	exp_sample_file.close();
	C_L_sample_file.close();

		//************* EX 01.3 - BUFFON EXPERIMENT
	int N_experiments = 5000;
	int N_blocks = 100;
	int M_max=500000; //here M is the number of throws after which a single data is computed
	double d=1,l=.7; //spatial parameters...

	std::ofstream buffon_file("Buffon.dat");

	std::vector<double> pi_vector(N_experiments);

	for(M= 100; M<M_max; M*=1.5){
		std::cout <<"Throws for each experiment: " <<M<< std::endl;
		for (int i=0; i<N_experiments; i++){ //doing N_experiments experiments each with M throws
			int N_accepted =0;
			for (int i =0; i<M; i++){ //in each loop we simulate a needle's throw
		//each needle is defined by the position of its center and by the angle theta it forms with the rules
				double center = rng.Rannyu(0.,.5*d);
				double sin_theta = get_uniform_sintheta(rng);
				if (center < .5*l*sin_theta)
					N_accepted++;
			}
			pi_vector[i] = 2*l*((double)M)/(d*N_accepted);//PI_estimate
		}
		avg_type res(compute_statistical_error(pi_vector, N_blocks));
		buffon_file << M << "\t" << res[0] << "\t" << res[2] << std::endl;
	}
	return 0;
}
