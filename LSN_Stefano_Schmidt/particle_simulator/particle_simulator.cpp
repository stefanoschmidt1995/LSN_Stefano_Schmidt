/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
CLASS PARTICLE SIMULATOR - DEFINITIONS
**************************************************************************
Definitions of class particle simulator declared in particle_simulator.h
*************************************************************************/
#include "../standard_header.h"
#include "particle_simulator.h"

particle_simulator::~particle_simulator(){
	for(auto it = measure_vector.begin(); it!= measure_vector.end(); it++)
		delete *it;
	for(auto it = hist_list.begin(); it!= hist_list.end(); it++)
		delete *it;
}

particle_simulator::particle_simulator(std::string filename, std::string initial_conf, std::string intial_old){
	std::cout << "Classic Lennard-Jones fluid        " << std::endl;
	std::cout << "Monte Carlo simulation             " << std::endl;
	std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << std::endl;
	std::cout << "The program uses Lennard-Jones units " << std::endl <<std::endl;

	//setting random generator
	rng = initialize_random(0);

	config_from_file(filename);

	//Prepare arrays for measurements
	measure_vector = std::vector<std::list<double>*>(n_props);
	for(auto it = measure_vector.begin(); it!= measure_vector.end(); it++){
		*it = new std::list<double>();
	}

	
	hist_list = std::list<std::vector<int>*>();
	
	//initialising variables for g measurement
	bin_size = (box/1.05)/(double)n_bins;

	//R3_point standard_point(3);
	state = std::vector<R3_point>(N, R3_point(3)); //CHECK IF IT'S RIGHT
	old_state = std::vector<R3_point>(N,R3_point(3));
	velocities = std::vector<R3_point>(N,R3_point(3));
	set_particle_configuration(initial_conf, intial_old);


	return;
}

void
particle_simulator::config_from_file(std::string filename){
	std::ifstream ReadInput;
	ReadInput.open(filename);
	
	ReadInput >> type;	
	std::cout << "Simulation type = " << type << std::endl;
	if(type == "NVT")
		std::cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << std::endl;

	ReadInput >> temp;
	if (temp<=0) beta = 1e10;
	else beta = 1.0/temp;
	std::cout << "Temperature = " << temp << std::endl;

	ReadInput >> N;
	std::cout << "Number of particles = " << N << std::endl;

	ReadInput >> rho;
	std::cout << "Density of particles = " << rho << std::endl;
	vol = (double)N/rho;
	box = pow(vol,1.0/3.0);
	std::cout << "Volume of the simulation box = " << vol << std::endl;
	std::cout << "Edge of the simulation box = " << box << std::endl;

	ReadInput >> rcut;
	std::cout << "Cutoff of the interatomic potential = " << rcut << std::endl;
		
	//Tail corrections for potential energy and pressure
	vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
	ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
	ptail = (double)N * ptail / vol;
	std::cout << "Tail correction for the potential energy = " << vtail << std::endl;
	std::cout << "Tail correction for the pressure         = " << ptail << std::endl; 

	ReadInput >> delta;

	ReadInput >> n_meas;

	ReadInput >> n_bins;
	
	std::cout << "Number of measure per block = " << n_meas << std::endl;
	std::cout << "It makes an histogram of g(r) with "<< n_bins<<" bins"<< std::endl<< std::endl;

	if(type=="NVT"){
		std::cout << "The program perform Metropolis moves with uniform translations" << std::endl;
		std::cout << "Moves parameter = " << delta << std::endl << std::endl;
	}
	if(type=="NVE"){
		std::cout << "The program integrates equation of motion with Verlet algorithm" << std::endl;
		std::cout << "Time step = " << delta << std::endl << std::endl;
	}

	
	ReadInput.close();

	
	return;
}

void
particle_simulator::set_particle_configuration(std::string filename, std::string initial_old){
	//Read configuration from file
	std::cout << "Read initial configuration from file config.0 " << std::endl << std::endl;
	std::ifstream ReadConf;
	ReadConf.open("config.0");
	for (int i=0; i<N; ++i){
		for (int j=0; j<3; j++){
			ReadConf >> (state[i]).at(j); //!!
			state[i][j] = Pbc( state[i][j] * box );
		}
	}
	ReadConf.close();

	if(type=="NVE" && initial_old!=""){
		std::cout << "Read initial old configuration from file "<< initial_old << std::endl << std::endl;
		std::ifstream ReadConf;
		ReadConf.open(initial_old);
		for (int i=0; i<N; ++i){
			for (int j=0; j<3; j++){
				ReadConf >> (old_state[i]).at(j); //!!
				old_state[i][j] = Pbc( old_state[i][j] * box );
			}
		}
		ReadConf.close();
	}

	if(type=="NVE" && initial_old==""){//preparing also velocities and old state
			//Prepare initial velocities
		std::cout << "Prepare random velocities with center of mass velocity equal to zero " << std::endl << std::endl;
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i=0; i<N; ++i){
			for(int j =0; j<3; j++){
				velocities[i].at(j) = rng.Rannyu() - 0.5;
				sumv[j] += velocities[i].at(j);
			}
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)N; //normalizing sumv
		double sumv2 = 0.0, fs;
		for (int i=0; i<N; ++i){
			for(int j =0; j<3; j++)
		     		velocities[i].at(j) -= sumv[j];
			sumv2 += norm(velocities[i]);
		   }
		sumv2 /= (double)N;
			//matching the right temperature for velocities and setting old_state
		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
		for (int i=0; i<N; ++i){
			for(int j =0; j<3; j++){
				velocities[i].at(j) *= fs;
				old_state[i].at(j) = Pbc(state[i].at(j) - velocities[i].at(j) * delta);
			}
		}
	}
	return;
}

void
particle_simulator::move(){
	moves++;
	if (type=="NVT")
		move_NVT();
	if (type=="NVE")
		move_NVE();
	return;

}

void
particle_simulator::move_NVT(){
	int o; //index of the particle to try to move
	double p, energy_old, energy_new;
	R3_point old(3), new_state(3);

	for(int i=0; i<N; ++i){
		o = (int)(rng.Rannyu()*N); //randomly selected particle

			//Old
		for (int i =0; i< 3; i++)
			old[i] = state[o][i];
		energy_old = particle_energy(old,o);

			//New
		for (int i =0; i< 3; i++) //try a metropolis random step
			new_state.at(i) = Pbc( old.at(i) + delta*(rng.Rannyu() - 0.5) ); //Eccolo l'inghippo!!!!
		energy_new = particle_energy(new_state,o);

			//Metropolis test
		p = exp(beta*(energy_old-energy_new)); //acceptance
		if(p >= rng.Rannyu()){ //Update the state if the move is accepted
			for(int j=0; j<3; j++)		
				state[o][j] = new_state[j];
			accepted++;
		}
	}
	//std::cout <<" success rate: "<< get_success_rate() << std::endl;
}
double
particle_simulator::get_success_rate(){
	if (type=="NVT")
		return accepted/(double)(moves*N);
	else
		return 1.;
}

double
particle_simulator::particle_energy(R3_point x, int ip){
	double ene=0.;
	double dr;

	for (int i=0; i<N; ++i){
		if(i != ip){
		// distance ip-i in pbc
			dr =0;
			for(int k=0; k<3; k++)
				dr+= pow(Pbc(x[k]- state[i][k]),2);
			dr = sqrt(dr);
			if(dr < rcut)
				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
		}
	}

	return 4.0*ene;
}

double
particle_simulator::Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}

void
particle_simulator::measure(){
	int bin;
	std::vector<int>* g_hist = new std::vector<int>(n_bins);
	double v = 0., w = 0., t=0., p=0., meas_temp=0.,meas_etot;
	double vij, wij;
	double dr=0;

		//cycle over pairs of particles for computing v, w
	for (int i=0; i<N-1; ++i){
		for (int j=i+1; j<N; ++j){

			// distance i-j in pbc
			dr =0;
			for(int k=0; k<3; k++)
				dr+= pow(Pbc(state[j][k]- state[i][k]),2);
			dr = sqrt(dr);

			//update of the histogram of g(r)
			if(dr < n_bins*bin_size){
				bin =0;
				while((double)(bin+1)*bin_size < dr)
					bin++;
				g_hist->at(bin) += 2; //incresing the counting on the bin in which dr fit
			}

			if(dr < rcut){
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
				// contribution to energy and virial
				v += vij;
				w += wij;
			}
		}					
	}
		//kinetic energy & temperature
	if (type=="NVE"){
		for (int i=0; i<N; ++i) 
			for (int j=0; j<3; ++j) 
				t += 0.5 * pow(velocities[i][j],2);
		meas_temp = (2.0 / 3.0) * t/(double)N; //Temperature
	}
	if (type == "NVT"){
		meas_temp = temp;
		t = 3./2. * N*temp; //equipartion theorem for canonical ensemble
	}

	v = v/(double)N + vtail; 		//Potential energy per particle
	w = 48*w/(double)N; 			//Virial per particle (without tail correction yet)
	t = t/(double)N; 			//Kinetic energy per particle
	meas_etot = (t+v); 			//Total energy per particle
	p = rho*temp + w*N/(3.*vol) + ptail; 	//pressure

	measure_vector[0]->push_back(meas_etot); 	//E_tot
	measure_vector[1]->push_back(t);	 	//E_kin
	measure_vector[2]->push_back(v);	 	//E_pot
	measure_vector[3]->push_back(w+ptail/rho);	//virial
	measure_vector[4]->push_back(meas_temp);	//temperature
	measure_vector[5]->push_back(p);		//pressure

	hist_list.push_back(g_hist);
	return;
}

std::vector<std::vector<double>>
particle_simulator::get_g_hist(){
	std::vector<avg_type> res(n_bins, avg_type(3));
	int n_blocks = hist_list.size()/n_meas;
	if (n_blocks <= 1){
		std::cerr << "Too few measures ("<<hist_list.size() <<") to provide the required precision of " << n_meas << " data per block"<<std::endl;
		return res;
	}
	std::vector<double> bin_vector(hist_list.size());
	for (int i =0; i <n_bins; i++){
		std::fill(bin_vector.begin(), bin_vector.end(), 0);
		int j=0;
		for(auto it = hist_list.begin(); it!= hist_list.end(); it++){
			bin_vector.at(j) = (double) (*it)->at(i);
			j++;
		}
		double r =  ((double)i+.5)*bin_size;
		double norm_factor = ((4.*M_PI/3.)*(pow(r+bin_size,3)-pow(r,3)));

		avg_type avg = compute_statistical_error(bin_vector, n_blocks, n_meas);
		
		res[i].at(0) = r;
		res[i].at(1) = avg[0]/ norm_factor;
		res[i].at(2) = avg[2] /norm_factor;
	}

	return res;
	
}

std::vector<avg_type>
particle_simulator::get_measures(){
	std::vector<avg_type> result(n_props, avg_type(3));
	int n_blocks = (measure_vector[0])->size()/n_meas;
	if (n_blocks <= 1){
		std::cerr << "Too few measures ("<<hist_list.size() <<") to provide the required precision of " << n_meas << " data per block"<<std::endl;
		return result;
	}
	for (unsigned int i = 0; i<n_props; i++){
		std::vector<double> temp_vector{ std::begin(*measure_vector[i]), std::end(*measure_vector[i]) };
		result[i] = compute_statistical_error(temp_vector, n_blocks);
	}
	return result;
}

void
particle_simulator::reset_measures(){
	for(auto it = hist_list.begin(); it!= hist_list.end(); it++)
		delete *it;
	hist_list.erase(hist_list.begin(), hist_list.end());
	for(auto it = measure_vector.begin(); it!= measure_vector.end(); it++)
		(*it)->erase((*it)->begin(), (*it)->end());
	moves =0;
	return;
}

void
particle_simulator::move_NVE(){ //Move particles with Verlet algorithm
	double temp_state;
	std::vector<R3_point> forces(N, R3_point(3));

	for(int i=0; i<N; ++i) //initializing forces
		forces.at(i) = force(i);

	for(int i=0; i<N; ++i){ //Verlet integration scheme
		for(int j =0; j<3; j++){ //creating new state & updating variables
			temp_state = Pbc( 2.0 * (state[i]).at(j) - (old_state[i]).at(j) + (forces[i]).at(j) * pow(delta,2) );
			(velocities[i]).at(j) = Pbc(temp_state - (old_state[i]).at(j))/(2.0 * delta);
			(old_state[i]).at(j) = (state[i]).at(j);
			(state[i]).at(j) = temp_state;

		}
	}
	return;
}



R3_point
particle_simulator::force(int ip){ //Compute forces as -Grad_ip V(r)
	double dr;
	R3_point dvec(3), f(3);
	

	for (int i=0; i<N; ++i){
		if(i != ip){
			for(int j =0; j<3; j++) // distance ip-i in pbc
				dvec.at(j) = Pbc( (state[ip]).at(j) - (state[i]).at(j) );
			dr = sqrt(norm(dvec));

			if(dr < rcut)
				for(int j =0; j<3; j++)	f[j] += dvec[j] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
		}
	}
  
	return f;
}

void
particle_simulator::set_new_temperature(double new_temp){
	if (type=="NVT"){
		temp = new_temp;
	}
	if(type == "NVE"){
		measure();
		double old_temp = (measure_vector[4])->back(); //recovering last measure of temperature
		double fs = sqrt(new_temp/old_temp); // v'/v = (T'/T)**(.5)

		for (int i=0; i<N; ++i){
			for(int j =0; j<3; j++){
				velocities[i].at(j) *= fs;
				old_state[i].at(j) = Pbc(state[i].at(j) - velocities[i].at(j) * delta);
			}
		}
		move();
		measure();
		std::cout << "   Temperature changed from " << old_temp << " to "<<  (measure_vector[4])->back() <<std::endl;
	}
	return;
}

bool
particle_simulator::find_eq_temperature(double temp_to_find, int n_iter){
	if (type=="NVT"){
		set_new_temperature(temp_to_find);
		return true;
	}
	if (type!="NVE") return false;
		//starting actual algorith for NVE
	double actual_temp=get_temperature(), old_temp=0, err, jump_factor=1.1;
	int count = 0;
	std::cout << "Doing equilibration: trying ot set T=" << temp_to_find <<" at equilibrium----"<< std::endl;
	do{
		if(abs(actual_temp - temp_to_find)<1) jump_factor =1.;
		set_new_temperature(temp_to_find*jump_factor);
		do{ //wait for the sistem to be in equilibrium, i.e. for the temperature not to change (much)...
			old_temp = get_temperature();
			for(int i =0; i<100; i++) move_NVE();
			actual_temp = get_temperature();
		}while(abs(actual_temp - old_temp)>1e-1);
		count ++;
		if (count>=50){
			std::cerr << "Unable to set the right temperature after 50 iterations... Sorry :(" <<std::endl;
			return false;
		}
			//checking compatibilty among the result and the request
		std::vector<double> stat_check(n_iter);
		for(int i =0; i<n_iter; i++){
			move_NVE();  
			stat_check[i] = get_temperature();
		}
		std::vector<double> avg_vect = compute_statistical_error(stat_check, n_iter/10);
		actual_temp = avg_vect[0];
		err = avg_vect[2];
		std::cout << "\tTemperature of the sistem averaged on "<<n_iter <<" iterations: " <<actual_temp <<" " << err 				<<std::endl;
		std::cout << "\tdistance from desired temperature: " << fabs(-actual_temp + temp_to_find)<<std::endl;
	}while(fabs(actual_temp - temp_to_find) > 2.*err); //fabs!! not abs 
	temp = actual_temp;
	return true;
}

double
particle_simulator::get_temperature(){
	if (type=="NVT")
		return temp;
	if (type =="NVE"){
		double t =0.;
		for (int i=0; i<N; ++i) 
			for (int j=0; j<3; ++j) 
				t += 0.5 * pow(velocities[i][j],2);
		return (2.0 / 3.0) * t/(double)N; //Temperature
	}
	return -1;//awful number to say that there is something wrong
}



void
particle_simulator::write_to_file(std::string filename, std::string filename_old){
	std::ofstream write_conf;

	std::cout << "Print configuration to file " << filename << std::endl << std::endl;
	write_conf.open(filename);
	for (int i=0; i<N; ++i){
		for(int j=0; j<3; j++)
			write_conf << Pbc(state[i][j]/box) << "\t";
		write_conf << std::endl;
	}
	write_conf.close();
	if (filename_old!=""){
		write_conf.open(filename_old);
		for (int i=0; i<N; ++i){
			for(int j=0; j<3; j++)
				write_conf << Pbc(old_state[i][j]/box) << "\t";
			write_conf << std::endl;
		}
		write_conf.close();
	}
	return;
}

void
particle_simulator::write_observable(std::string filename, unsigned int index){
	if (index>=n_props) return;
	std::fstream write_conf;
	write_conf.open(filename, std::ios::app);
	write_conf << measure_vector[index]->back() << std::endl;
	write_conf.close();
	return;

}

void
particle_simulator::write_hist(std::string filename){
	std::ofstream write_conf;
	write_conf.open(filename);
	std::vector<std::vector<double>> hist_to_write(get_g_hist());
	
	for (int i = 0; i<n_bins; i++)
		write_conf << (hist_to_write[i]).at(0) << " " << (hist_to_write[i]).at(1) << " " << (hist_to_write[i]).at(2) <<std::endl;

	write_conf.close();
	return;
}

std::vector<double>
particle_simulator::get_config_parameters(){ 
	std::vector<double> result(8);
	if (type=="NVT")
		result[0] = 1.;
	if (type=="NVE")
		result[0] = 0.;
	result[1] = temp;
	result[2] = N;
	result[3] = rho;
	result[4] = rcut;
	result[5] = delta;
	result[6] = n_meas;
	result[7] = n_bins;	
	return result;
}

void
particle_simulator::write_avg_observable(std::string filename){
	std::vector<avg_type> averages = get_measures();
	std::fstream file_handler(filename, std::ios::app);
	//header: [rho, [E_tot], [E_kin], [E_pot], [virial], [Temp], [P], moves/time]
	file_handler << rho << "\t";
	for (unsigned int i=0; i < averages.size(); i++){
		file_handler << averages.at(i)[0] << "\t"<< averages.at(i)[2] << "\t";
	}
	if (type=="NVE")
		file_handler << moves*delta;
	if (type=="NVT")
		file_handler << moves;
	file_handler << std::endl;
	file_handler.close();
	return;

}









