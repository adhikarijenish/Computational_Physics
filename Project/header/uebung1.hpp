#pragma once
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<complex>
#include<functional>

using complex = std::complex<double>;

void creat_lattice(std::vector<std::vector<std::vector<complex>>> &lattice,int L, std::mt19937 &gen);
complex plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,std::vector<int> n,int mu,int nu);
complex g_plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,const double beta);
void generate_Phi(std::vector<std::vector<std::vector<complex>>> & Phi,std::mt19937 &gen);
int delta_f(int n,int m);     
int delta_f(std::vector<int> n,std::vector<int> m);                 
complex fermion(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0);
complex fermion_dagger(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0);
std::vector<std::vector<std::vector<complex>>> f_M(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0);
std::vector<std::vector<std::vector<complex>>> f_Mdag(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0);
std::vector<std::vector<std::vector<complex>>> f_M_f_Mdag(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0);
double normsquared(const std::vector<std::vector<std::vector<complex>>> & psi);
std::vector<double> scalar_product(std::vector<std::vector<std::vector<complex>>> p, std::vector<std::vector<std::vector<complex>>> t);
void assign_add_mul(std::vector<std::vector<std::vector<complex>>> &x,const std::vector<std::vector<std::vector<complex>>> &p,const std::vector<double> alpha);
void assign_mul_add(std::vector<std::vector<std::vector<complex>>> &p,const std::vector<std::vector<std::vector<complex>>> &r,const std::vector<double> beta);
std::vector<std::vector<std::vector<complex>>> cg(const std::function< std::vector<std::vector<std::vector<complex>>>(const std::vector<std::vector<std::vector<complex>>> &,const std::vector<std::vector<std::vector<complex>>> &,const double) > & M, 
    std::vector<std::vector<std::vector<complex>>> & psi, const std::vector<std::vector<std::vector<complex>>> & U_gauge, double m0, size_t Max_Iterations, double epsilon);





void GetUserParam( int argc, char *argv[],std::string & filename, int & Pairs, int & X, bool & save);