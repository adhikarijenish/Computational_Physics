#include "../header/uebung1.hpp"
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<functional>

using complex = std::complex<double>;

std::string filename;

int main(int argc, char *argv[]){
    size_t seed = 123456789;
    std::mt19937 gen(seed);
    int L = 3;
    double m0 = 1;

    std::vector<std::vector<std::vector<complex>>> lattice(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    creat_lattice(lattice,L,gen);
    double beta = 1;
    std::vector<std::vector<std::vector<complex>>> test_phi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> test_chi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    generate_Phi(test_phi,gen);
    generate_Phi(test_chi,gen);
    //std::vector<std::vector<std::vector<complex>>> inverse = cg(f_M_f_Mdag,test_phi,lattice,m0,10000,1e-15);
    // check diagonal;

    std::vector<std::vector<std::vector<complex>>> ferm = f_M(test_chi,lattice,m0);
    std::vector<std::vector<std::vector<complex>>> fermd = f_Mdag(test_chi,lattice,m0);

    for(int i = 0;i<L;i++){
        for(int j = 0;j<L;j++){
                std::cout<<"Fermion Matrix: " <<ferm[i][j][0] <<" Hermitian conj: "<<std::conj(fermd[i][j][0])<<" \n";
                }
            }
        
    
    
    return 0;
}