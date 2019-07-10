/*************************************************************************
STEFANO SCHMIDT - NUMERICAL SIMULATION LABORATORY
STANDARD HEADER
**************************************************************************
Declaration of some useful functions used in many exercises declared in ./general_code/useful_functions.cpp
Standard #include command of some c++ standard library and of class random
*************************************************************************/

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <numeric>
#include <list>
#include <random>
#include <algorithm>
#include"./general_code/random.h"

typedef std::vector<int> point;
typedef std::vector<double> R3_point;
typedef std::vector<double> avg_type;

Random initialize_random(); //helper to initialize the random generator with a standard seed taken from files in general_code
Random initialize_random(int);  //above but with the option of specifying a seed
double norm(point); //L 2 modulus of the point (n dimensions) 
double norm(R3_point); //L2 modulus of the point (n dimensions)
R3_point operator+(const R3_point& , const R3_point&); //sum for R3_points
R3_point operator-(const R3_point& , const R3_point&); //difference for R3_points
std::vector<double> compute_statistical_error(const std::vector<double>& values, int N_blocks);
	//compute statistical errors with blocking methods using N_blocks and all avaible data in values
	//returns a vector with avg|var|avg stddev
std::vector<double> compute_statistical_error(const std::vector<double>& values, int N_blocks,  int L);
	//compute statistical errors with blocking methods using N_blocks and L data per block
