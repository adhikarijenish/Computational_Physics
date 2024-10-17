#pragma once
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>

/*
    ******** Generate Pairs of random number (x,y) ********************


*/ 
std::vector<double> random_number( std::mt19937 & engine);

/*
 ***************** Take Mean of 1D-Vector *************
*/ 
double mean(std::vector<double> &vector_1D);