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
    size_t seed = 1234567;
    std::mt19937 gen(seed);
    int L = 4;
    double m0 = 1;

    std::vector<std::vector<std::vector<complex>>> lattice(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    creat_lattice(lattice,L,gen);
    
    std::vector<std::vector<std::vector<complex>>> test_chi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> test_phi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));

    generate_Phi(test_chi,gen);
    test_phi = f_M(test_chi,lattice,m0);

    std::vector<std::vector<std::vector<complex>>> inverse = cg(f_M_f_Mdag,test_phi,lattice,m0,10,1e-5);
    
/*/ check if M and M_dagger diagonal element match


    std::vector<std::vector<std::vector<complex>>> v(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2,complex(1.0,0.0))));
    std::vector<std::vector<std::vector<complex>>> ferm ;
    std::vector<std::vector<std::vector<complex>>> fermd;
    for (int t = 0; t < L; t++) {      
        for (int x = 0; x < L; x++) {  
            for (int a = 0; a < 2; a++) {  
                ferm = f_M(v,lattice,m0);
                fermd = f_Mdag(v,lattice,m0);
                if (t == x){std::cout<< "Fermion: "<<ferm[t][x][a]<< " Anti-Fermion: "<<fermd[t][x][a]<<" \n";}
            }
        }
    } /**/       
    
    
    return 0;
}