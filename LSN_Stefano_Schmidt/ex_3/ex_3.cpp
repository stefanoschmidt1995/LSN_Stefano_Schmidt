/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 3
*************************************************************************/
#include "../standard_header.h"

int main(){
		//************* EX 03.1 - OPTION PRICING
	int M = 1e6; //number of events generated
	int T_step = 100; // #of time step to sample the GBM

		//model parameters
	double S_0 = 100.;
	double T =1.;
	double K =100; 
	double r =0.1 ;
	double sigma =0.25;

	std::vector<double> C_dir(M), C_sam(M), P_dir(M), P_sam(M); //defining vectors to hold the data
	Random rng = initialize_random();

	//option price is Monte Carlo computed by
	//C(S(0),0) = E[ exp(−rT)* S(T) − K ]
	//with S(T) given by a Geometric brownian motion

		//computing S(T) and computing the prices
	std::cout << "Doing option pricing with "<<M<< " events generated."<<std::endl;
	for(int i =0; i<M; i++){
			//direct sampling 
		double S_dir, S_sam = S_0;
		double RW_pos = 0;
		S_dir = S_0 * exp( (r-0.5*sigma*sigma)*T + sigma * rng.Gauss(0.,sqrt(T) ));
			//RW sampling
		for(int t_step = 0; t_step <T_step; t_step++){
			double delta_t = T /T_step;
			RW_pos = rng.Gauss(RW_pos, sqrt(delta_t)); //advancing RW
			S_sam = S_0 * exp( (r-0.5*sigma*sigma)*(t_step+1)*delta_t + sigma * RW_pos );
			//S_sam = S_sam * exp( (r-0.5*sigma*sigma)*delta_t + sigma * rng.Gauss(0., sqrt(delta_t) );
		}
		C_dir[i] = exp(-r*T)*std::max(S_dir - K, 0.); //price for Call Option (direct sampling)
		C_sam[i] = exp(-r*T)*std::max(S_sam - K, 0.); //price for Call Option (stepwise sampling)
		P_dir[i] = exp(-r*T)*std::max(-S_dir + K, 0.); //price for Put Option (direct sampling)
		P_sam[i] = exp(-r*T)*std::max(-S_sam + K, 0.); //price for Put Option (stepwise sampling)
	}
		//computing averages and errors with data blocking and printing to stdout
	std::vector<double> avg_C_dir = compute_statistical_error(C_dir, 100);
	std::vector<double> avg_C_sam = compute_statistical_error(C_sam, 100);
	std::cout << "\tC dir: " << avg_C_dir[0] << "\t" << avg_C_dir[2] <<std::endl;
	std::cout << "\tC sam: " << avg_C_sam[0] << "\t" <<avg_C_sam[2] <<std::endl;

	std::vector<double> avg_P_dir = compute_statistical_error(P_dir, 100);
	std::vector<double> avg_P_sam = compute_statistical_error(P_sam, 100);
	std::cout << "\tP dir: " << avg_P_dir[0] << "\t" <<avg_P_dir[2] <<std::endl;
	std::cout << "\tP sam: " << avg_P_sam[0] << "\t" <<avg_P_sam[2] <<std::endl;

		//printing everything to file
	int L = 50; //L holds the number of data within one block

	std::ofstream out_file;
	out_file.open("prices.dat");

	out_file << "N \t C direct \t sigma C direct \t C samples \t sigma C sample \t P direct \t sigma P direct \t P sample \t sigma P sample" << std::endl; //header for prices.dat 

	for(int N = 2; N<M/L; N*=1.3){
		if(N<=5) N++;
		avg_C_dir = compute_statistical_error(C_dir, N, L);
		avg_C_sam = compute_statistical_error(C_sam, N, L);
		avg_P_dir = compute_statistical_error(P_dir, N, L);
		avg_P_sam = compute_statistical_error(P_sam, N, L);
		out_file << N << "\t" << avg_C_dir[0] << "\t" <<  avg_C_dir[2] << "\t" << avg_C_sam[0] << "\t" <<  avg_C_sam[2]<< "\t" << avg_P_dir[0] << "\t" <<  avg_P_dir[2]<< "\t" << avg_P_sam[0] << "\t" <<  avg_P_sam[2]<< std::endl;
	}
	out_file.close();

	return 0;
}
//In this computation there are two elements: 	the model for the time evolution S(t) of the price of an asset
//						the particular kind of option that is chosen

//OSS: If we use some more difficult model for the evolution of the price (i.e. the value of the price at time T) the Monte Carlo Method for computing the price of the option remains the same.
//OSS2: If we use some different option we compute in a different way the average gain which is suitable with the option considered.

