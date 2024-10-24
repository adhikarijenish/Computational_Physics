#pragma once
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>

/*
    ******** Generate Pairs of random number (x,y) ********************


*/ 
void random_number(std::mt19937 & engine,std::vector<double> &random_2D);

/*
 ***************** Take Mean of 1D-Vector *************
*/ 
double mean(std::vector<double> &vector_1D);

double sum(std::vector<double> &vector_1D);

void GetUserParam( int argc, char *argv[],std::string & filename, int & Pairs, int & X, bool & save);