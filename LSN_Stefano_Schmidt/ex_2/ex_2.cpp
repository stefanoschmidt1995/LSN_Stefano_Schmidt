/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
EXERCISE 2
*************************************************************************/
#include "../standard_header.h"

R3_point
get_random_direction(Random* rng, int ndim = 3){
	R3_point x(ndim);
	double x_norm=0;
	do{
		for(auto it = x.begin(); it!=x.end(); it++)
			*it = rng->Rannyu(-1,1);
		x_norm = norm(x);
	} while(x_norm >1.);
	for(auto it = x.begin(); it!=x.end(); it++)
			*it = *it/sqrt(x_norm);
	return x;
}

double
sample_from_p(Random* rng){ //easy function to sample from probability distribution p(x) = 2./3.* (1-x**2)
	double x,y;
	do{
		y = rng->Rannyu(0,2./3.);
		x = rng->Rannyu();
	} while(y > 2./3. * (1.-x*x));
	return x;
}



int main(){
	Random rng = initialize_random();
		//************* EX 02.1
	int M = 1e7; //# values to be computed in the Montecarlo sum
	int L = 50; //constant values of blocks
	int N_max = M/L;
	std::vector<double> unif_dist(M), imp_sampl(M);
	std::vector<double> unif_result(2), imp_result(2);
		//importance sampling is done with the following PDF
		// p(x) = 2/3 * (1-x**2)
		// F(x) = 3/2 * (x-1/3 * x**3)   cumulative
		// x = (-1 + i sqrt(3) - (1 + i sqrt(3)) (-y + sqrt(-1 + y^2))^(2/3))/(2 (-y + sqrt(-1 + y^2))^(1/3)) 
		//inverse of cumulative is awful!! Impossibile to compute seriously: we thus do a random sampling with rejection

	for(int i =0; i<M ; i++){
		unif_dist[i]= M_PI/2. * cos(M_PI/2. * rng.Rannyu());
		double imp_x = sample_from_p(&rng);
		imp_sampl[i] = M_PI/3. * cos(M_PI/2. * imp_x) /(1.-imp_x*imp_x);
	}

	std::ofstream unif_sample_file("unif_sampling.dat");
	std::ofstream imp_sample_file("imp_sampling.dat");
	std::cout << "Computing the integral" << std::endl;
	for(int N = 2; N<N_max; N*=1.3){
		if(N<=5) N++;
			//choosign subsample of data to average
			//computing & saving averages
		unif_result = compute_statistical_error(unif_dist, N, L);
		unif_sample_file <<N << "\t "  << unif_result[0] <<"\t " <<unif_result[2] <<std::endl;
		imp_result = compute_statistical_error(imp_sampl, N,L);
		imp_sample_file << N << "\t " << imp_result[0] <<"\t " <<imp_result[2] <<std::endl;
	}
	unif_sample_file.close();
	imp_sample_file.close();

		//************* EX 02.2 - 3D RANDOM WALK
	int N_steps=1e3, N_simulations = 2e4;
	int ndim = 3; //dimension for lattice RW

		//defining RW variables
	point pos(ndim);
	R3_point cont_pos(ndim);
	std::vector<point> RW_pos(N_simulations, pos);
	std::vector<R3_point> cont_RW_pos(N_simulations, cont_pos);
	std::vector<double> r_RW(N_simulations), r_cont_RW(N_simulations); //vectors to hold at each time step the values of r for each simulation (vector to average at each step)
	int N_blocks = 100;

		//dealing with files
	std::ofstream RW_file("RW.dat"); //RW.dat holds a sample of walker's positions in the lattice RW
	std::ofstream cont_RW_file("continuos_RW.dat"); //RW.dat holds a sample of walker's positions in the continuos RW
	std::ofstream r_file("r_average.dat"); //hold the square root of avg value of r**2 with uncertainties
	r_file << "time step\t" << "r_RW\t" << "sigma_r_RW\t" << "r_cont_RW\t" << "sigma_r_cont_RW" <<std::endl; //header of r_file

	std::cout << "Performing Random Walk" << std::endl;
	for(int i =0; i<N_steps; i++){
		if (i%25 == 24 )
			std::cout << "\tRW: doing time step " << i+1 << std::endl;
		for(int n = 0; n <N_simulations; n++){ //loops on simulations
				//extracting steps
			int axis = ((int) (10000*rng.Rannyu())) % (ndim); //choosing an axis at random
			int displacement = (2*(int)(rng.Rannyu()<0.5))-1; //useful also for byased RW
			R3_point direction = get_random_direction(&rng);

				//updating positions
			RW_pos[n][axis] = RW_pos[n][axis] + displacement; //updating position of lattice RW
			for(int j=0; j<ndim; j++) //updating position of continuous RW
				cont_RW_pos[n][j] = cont_RW_pos[n][j] + direction[j];
			
				//computing r for each simulation
			r_RW[n] = norm(RW_pos[n]);
			r_cont_RW[n] = norm(cont_RW_pos[n]);

		}
			//computing average over r
		std::vector <double> RW_result = compute_statistical_error(r_RW,N_blocks);
		std::vector <double> cont_RW_result = compute_statistical_error(r_cont_RW,N_blocks);
		if (i%5==0)
			r_file << i+1 <<" " << sqrt(RW_result[0]) << " " << RW_result[2] << " " << sqrt(cont_RW_result[0]) << " " << cont_RW_result[2] <<std::endl;

			//writing to file the first RW
		for(int j=0; j<ndim; j++){ //updating position of continuous RW
			cont_RW_file << cont_RW_pos[0][j] << " ";
			RW_file << RW_pos[0][j] << " ";
		}

		RW_file <<std::endl;
		cont_RW_file << std::endl;

	}
	r_file.close();
	RW_file.close();
	cont_RW_file.close();
	return 0;
}










