#pragma once
#include<vector>
#include<random>
/* Generate Random number in intervall [0,1]*/
void random_number_generator(std::vector<double> &random_number,std::mt19937 & engine, const int & len);
/* Generate random Spin-Distribution in 2D */
void random_2D_Spin(std::vector<std::vector<int>>& Spin_2D,std::mt19937 & engine,const int & N);
/* Generate random Spin-Distribution in 3D */
void random_3D_Spin(std::vector<std::vector<std::vector<int>>> & Spin_3D,std::mt19937 & engine,const int & N);

/* Calculate Energy of 2D Ising-Model */
double Energy_sum_2D(const std::vector<std::vector<int>>& Spin_config);
/* Calculate Energy of 3D Ising-Model */
double Energy_sum_3D(const std::vector<std::vector<std::vector<int>>> & Spin_3D);
/* Calculate Energy difference in 2D Ising-Model */
int delta_E_2D(const std::vector<std::vector<int>>& Spin_2D,const int & i,const int & j); 
/* Calculate Energy difference in 3D Ising-Model */
int delta_E_3D(const std::vector<std::vector<std::vector<int>>> & Spin_3D,const int & i,const int & j,const int & k);

/* Calculate Sum of a Integer-vector as Integer*/
int sum_vec(const std::vector<int>& Spin_1D);

/* Calculate Sum of a Integer-2D-vector as Integer*/
int vector_2D_sum(const std::vector<std::vector<int>> & Spin_2D,std::vector<float> & z);

/* Calculate Sum of a Double-vector as Double*/
double sum_vec(const std::vector<double>& Spin_1D);
/* Calculate mean of 1-D Vector */
double mean(const std::vector<double>& data);

/* Calculate mean of 1-D Vector */
double mean(const std::vector<int>& data);

/* Calculate mean of 2-D Vector */
double mean(const std::vector<std::vector<int>> & data);

/* Calculate mean of 2-D Vector */
double mean(const std::vector<std::vector<std::vector<int>>> & data);

/* Calculate mean-square of 1-D Vector */
double meansquare(const std::vector<double>& data);

/* Calculate variance of 1-D Vector */
double variance(const std::vector<double>& data);

/* Compute a Sweep of 2D-Ising-Model as discribe in Algorithm 1 */
void Ising_2D_Sweep(std::vector<std::vector<int>>& Spin_2D,std::mt19937 &engine,const double & B,double& E);

/* Compute a Sweep of 3D-Ising-Model as discribe in Algorithm 1 */
void Ising_3D_Sweep(std::vector<std::vector<std::vector<int>>> & Spin_3D,std::mt19937 &engine,const double & B,double& E);