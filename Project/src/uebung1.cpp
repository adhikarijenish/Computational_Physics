#include "../header/uebung1.hpp"
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include <cstring>
#include <complex>
#include <functional>

using complex = std::complex<double>;

void creat_lattice(std::vector<std::vector<std::vector<complex>>> &lattice,int L, std::mt19937 &gen){

    std::normal_distribution<double> dist(0.0,1.0);
    for (int x=0;x<L;x++){
        for(int y=0;y<L;y++){
            for(int z=0;z<2;z++){
                lattice[x][y][z] = std::exp(complex(dist(gen),dist(gen)));
            }
        }
    }
}

void generate_Phi(std::vector<std::vector<std::vector<complex>>> & Phi,std::mt19937 &gen){
    int L = Phi.size();
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
            for(int z = 0;z<2;z++){
                std::normal_distribution<double> dist(0.0,1.0);
                Phi[x][y][z] = complex(dist(gen),dist(gen));
            }
        }
    }

}

complex plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,std::vector<int> n,int mu,int nu){
    int x = n[0];
    int y = n[1];
    int L = lattice.size();
    complex p;
    p = lattice[x][y][mu]*lattice[(x+1)%L][y][nu]*std::conj(lattice[x][(y+1)%L][mu])*std::conj(lattice[x][y][nu]);
    
    return p;
    
}

complex g_plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,const double beta){

    complex sum = 0;
    int N = lattice.size();
    std::vector<int> n;
    int mu = 0; int nu = 1;

    for (int x=0;x<N;x++){
        for(int y=0;y<N;y++){
            n.push_back(x);n.push_back(y);
            sum += plaquette(lattice,n,mu,nu);
            n.clear();
        }

    }
    sum = beta*(1-sum.real());
    return sum;
}
int delta_f(int n,int m){
    if (n==m){
    return 1;}
    return 0;
}

int delta_f(std::vector<int> n,std::vector<int> m){
    int x = n[0];
    int y = n[1];
    int x_p = m[0];
    int y_p = m[1];
    if (x==x_p && y==y_p){
        return 1;}
    return 0;
}

complex fermion(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0){
    std::vector<std::vector<std::vector<complex>>> sigma(3,std::vector<std::vector<complex>>(2,std::vector<complex>(2)));

    sigma[0] ={{{0.0,0.0},{1.0,0.0}},{{1.0,0.0},{0.0,0.0}}};
    sigma[1] ={{{0.0,0.0},{0.0,-1.0}},{{0.0,1.0},{0.0,0.0}}};
    sigma[2] ={{{1.0,0.0},{0.0,0.0}},{{0.0,0.0},{-1.0,0.0}}};
    int L = U_gauge.size();
    complex out;
    int x = n[0];
    int y = n[1];
    int x_p = n_prime[0];
    int y_p = n_prime[1];
    int mu = 0; int nu =1;
    auto new_n_prime = n_prime;

    out = complex((m0+2)*delta_f(alpha,beta)*delta_f(n_prime,n),0.0);
    new_n_prime[mu] = n_prime[mu]+1;
    out -=0.5*(complex(1.0,0.0)-sigma[mu][alpha][beta])*U_gauge[x][y][mu]*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[mu] = n_prime[mu]-1;
    out -=0.5*(complex(1.0,0)+sigma[mu][alpha][beta])*std::conj(U_gauge[(x_p-1+L)%L][y_p][mu])*complex(delta_f(new_n_prime,n),0.0);
    
    new_n_prime = n_prime;
    new_n_prime[nu] = n_prime[nu]+1;
    out -=0.5*(complex(1.0,0.0)-sigma[nu][alpha][beta])*U_gauge[x][y][nu]*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[nu] = n_prime[nu]-1;
    out -=0.5*(complex(1.0,0.0)+sigma[nu][alpha][beta])*std::conj(U_gauge[(x_p)%L][(y_p-1+L)%L][nu])*complex(delta_f(new_n_prime,n),0.0);

    return out;

}

complex fermion_dagger(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0){
    
    std::vector<std::vector<std::vector<complex>>> sigma(3,std::vector<std::vector<complex>>(2,std::vector<complex>(2)));
    sigma[0] ={{{0.0,0.0},{1.0,0.0}},{{1.0,0.0},{0.0,0.0}}};
    sigma[1] ={{{0.0,0.0},{0.0,-1.0}},{{0.0,1.0},{0.0,0.0}}};
    sigma[2] ={{{1.0,0.0},{0.0,0.0}},{{0.0,0.0},{-1.0,0.0}}};
    int L = U_gauge.size();
    complex out;
    int x = n[0];
    int y = n[1];
    int x_p = n_prime[0];
    int y_p = n_prime[1];
    int mu = 0; int nu =1;

    auto new_n_prime = n_prime;

    out = complex((m0+2)*delta_f(alpha,beta)*delta_f(n_prime,n),0.0);

    new_n_prime[mu] = n_prime[mu]+1;
    out -=0.5*(complex(1.0,0.0)-sigma[mu][alpha][beta])*std::conj(U_gauge[x][y][mu])*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[mu] = n_prime[mu]-1;
    out -=0.5*(complex(1.0,0.0)+sigma[mu][alpha][beta])*U_gauge[(x_p-1+L)%L][y_p][mu]*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[nu] = n_prime[nu]+1;
    out -=0.5*(complex(1.0,0.0)-sigma[nu][alpha][beta])*std::conj(U_gauge[x][y][nu])*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[nu] = n_prime[nu]-1;
    out -=0.5*(complex(1.0,0.0)+sigma[nu][alpha][beta])*U_gauge[(x_p)%L][(y_p-1+L)%L][nu]*complex(delta_f(new_n_prime,n),0.0);
    
    return out;

}

std::vector<std::vector<std::vector<complex>>> f_M(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0){
    int L = U_gauge.size();
    std::vector<std::vector<std::vector<complex>>> fermion_matrix(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    complex temp;
    std::vector<int> n;
    std::vector<int> n_prime;
    for(int x_p = 0; x_p<L;x_p++){
        for(int y_p = 0; y_p<L;y_p++){
            for(int alpha = 0; alpha<2;alpha++){
                temp = complex(0.0,0.0);
                for(int x = 0;x<L;x++){
                    for(int y = 0;y<L;y++){
                        for(int beta = 0;beta<2;beta++){
                            n_prime.push_back(x_p);
                            n_prime.push_back(y_p);
                            n.push_back(x);
                            n.push_back(y);

                            temp += fermion(U_gauge,n,n_prime,alpha,beta,m0)*Phi[x][y][beta];

                            n_prime.clear();n.clear();
                        }
                    }
                }
            fermion_matrix[x_p][y_p][alpha] = temp;
            }
        }
    }

    return fermion_matrix;

}
std::vector<std::vector<std::vector<complex>>> f_Mdag(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0){
    int L = U_gauge.size();
    std::vector<std::vector<std::vector<complex>>> fermion_matrix(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    complex temp;
    std::vector<int> n;
    std::vector<int> n_prime;
    for(int x_p = 0; x_p<L;x_p++){
        for(int y_p = 0; y_p<L;y_p++){
            for(int alpha = 0; alpha<2;alpha++){
                temp = complex(0.0,0.0);
                for(int x = 0;x<L;x++){
                    for(int y = 0;y<L;y++){
                        for(int beta = 0;beta<2;beta++){
                            n_prime.push_back(x_p);
                            n_prime.push_back(y_p);
                            n.push_back(x);
                            n.push_back(y);

                            temp += fermion_dagger(U_gauge,n,n_prime,beta,alpha,m0)*Phi[x][y][beta];

                            n_prime.clear();n.clear();
                        }
                    }
                }
            fermion_matrix[x_p][y_p][alpha] = temp;
            }
        }
    }

    return fermion_matrix;
    
}

std::vector<std::vector<std::vector<complex>>> f_M_f_Mdag(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0){

    return f_M(f_Mdag(Phi,U_gauge,m0),U_gauge,m0);
}

std::vector<double> normsquared(const std::vector<std::vector<std::vector<complex>>> & psi){
    int L = psi.size();
    std::vector<double> normal;
    double norm_x = 0.0;
    double norm_y = 0.0;
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
                norm_x += (psi[x][y][0]*std::conj(psi[x][y][0])).real();
                norm_y += (psi[x][y][1]*std::conj(psi[x][y][1])).real();
            }
        }
    normal.push_back(std::abs(norm_x));
    normal.push_back(std::abs(norm_y));
    
    return normal;
}

std::vector<double> scalar_product(std::vector<std::vector<std::vector<complex>>> p, std::vector<std::vector<std::vector<complex>>> t){
    int L = p.size();
    complex x_dir = complex(0.0,0.0);
    complex y_dir = complex(0.0,0.0);
    std::vector<double> product;
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
            x_dir += (std::conj(p[x][y][0])*t[x][y][0]);
            y_dir += (std::conj(p[x][y][1])*t[x][y][1]);
        }
    }

    product.push_back((x_dir).real());
    product.push_back((y_dir).real());
    return product;
}

void assign_add_mul(std::vector<std::vector<std::vector<complex>>> &x,const std::vector<std::vector<std::vector<complex>>> &p,const std::vector<double> alpha){
    int L = x.size();
    for(int x_p = 0;x_p<L;x_p++){
        for(int y_p = 0;y_p<L;y_p++){
            x[x_p][y_p][0] = x[x_p][y_p][0] + alpha[0]*p[x_p][y_p][0];
            x[x_p][y_p][1] = x[x_p][y_p][1] + alpha[1]*p[x_p][y_p][1];
        }
    }

}

void assign_mul_add(std::vector<std::vector<std::vector<complex>>> &p,const std::vector<std::vector<std::vector<complex>>> &r,const std::vector<double> beta){
    int L = p.size();
    for(int x_p = 0;x_p<L;x_p++){
        for(int y_p = 0;y_p<L;y_p++){ // p{i+1} = r{i+1} + beta{i+1} * p{i}
            p[x_p][y_p][0] = r[x_p][y_p][0] + beta[0]*p[x_p][y_p][0];
            p[x_p][y_p][1] = r[x_p][y_p][1] + beta[1]*p[x_p][y_p][1];
        }
    }

}


std::vector<std::vector<std::vector<complex>>> cg(const std::function< std::vector<std::vector<std::vector<complex>>>(const std::vector<std::vector<std::vector<complex>>> &,const std::vector<std::vector<std::vector<complex>>> &,const double) > & M, 
    std::vector<std::vector<std::vector<complex>>> & psi, const std::vector<std::vector<std::vector<complex>>> & U_gauge, double m0, size_t Max_Iterations, double epsilon){
	
    std::vector<double> rsqr, err;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> max_err;
    std::vector<double> temp_sc_product;
    int L = psi.size();

	// maximum norm squared of r for which the CG stops (remember stopping condition ||r{i}||/||psi|| < epsilon => (r{i},r{i}) < epsilon^2 * (psi,psi))
    max_err =normsquared(psi);
	double max_err_x = std::pow(epsilon,2.0)*max_err[0];
	double max_err_y = std::pow(epsilon,2.0)*max_err[1];

	// vectors r,p,t,x have the same meaning as defined in lecture
    std::vector<std::vector<std::vector<complex>>> r(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> p(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> t(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> x(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2,complex(0.0,0.0))));  // initial guess for x, zero here
    

    
	//initial residual r{0} = psi - A * x{0} (x{0} = 0 here)
    r = psi;

	// p{0} = r{0}
    p = psi;

    // rsqr = (r{0},r{0})
    rsqr = normsquared(r);

    for(size_t i = 0; i < Max_Iterations; i++) {
    	// t = A * p
    	t = M(p,U_gauge,m0);

    	//alpha{i} = (r{i},r{i})/(p{i},t{i}) (remember we stored (r{i},r{i}) in rsqr)
        temp_sc_product = scalar_product(p,t);
    	alpha.push_back(rsqr[0]/temp_sc_product[0]);
    	alpha.push_back(rsqr[1]/temp_sc_product[1]);

    	// x{i+1} = x{i} + alpha{i} * p{i}
    	assign_add_mul(x,p,alpha);

    	// r{i+1} = r{i} - alpha{i} * t{i}
        alpha[0] = -1.0*alpha[0];
        alpha[1] = -1.0*alpha[1];

    	assign_add_mul(r,t,alpha);

    	// err = (r{i+1},r{i+1})
    	err = normsquared(r);

		 std::cout << "Error_x[" << i << "] = " << err[0]<< " Error_y: "<<err[1] << std::endl;

    	// check for convergence
    	if(err[0] < max_err_x && err[1] < max_err_y) {
    		std::cout << "Required Precision Reached: " << std::endl;
    		return x ;
    	}

    	// beta{i+1} = (r{i+1},r{i+1})/(r{i},r{i}) = err/rsqr
    	beta.push_back(err[0]/rsqr[0]);
    	beta.push_back(err[1]/rsqr[1]);

    	// p{i+1} = r{i+1} + beta{i+1} * p{i}
		assign_mul_add(p,r,beta);

    	// ensure that rsqr is (r{i},r{i}) for the next iteration (i.e. i -> i+1)
    	rsqr = err;

        alpha.clear();
        beta.clear();
    }

    throw std::runtime_error("Error: CG did not converged: Error_x = " + std::to_string(rsqr[0])+"  Error_y = " + std::to_string(rsqr[1]));
}




void GetUserParam( int argc, char *argv[],std::string & filename, int & Pairs, int & X, bool & save){

/* Variablen: */
    int i;
    char* endptr;
    const char usage[] = 
        "boundary [-f <File Name> -p <Pairs> -x <Experiment> -s <Save Radius {0/1}>]";
    const char error_message[] =
        "# FEHLER(GetuserParam): falsche Option: ";

    if (argc>1) { /* falls es ueberhaupt Parameter gibt ... */
        for (i=1; i<argc; i++){
            /* parameter 2 Charakter lang und sollte mit '-' anfaengen ... */
            if ( (std::strlen(argv[i])==2) && (argv[i][0] == '-') ) { 
                switch (argv[i][1]) { 
                    case 'f':
                        filename = argv[++i];
                        
                        break;
                    case 'p':
                        Pairs = strtod( argv[++i], &endptr);
                        
                        break;
                    case 'x':
                        X = strtod( argv[++i], &endptr);
                        
                        break;
                    case 's':
                        save = strtod( argv[++i], &endptr);
                        
                        break;
                    default:
                        break;
                }
            } else {
                std::cout << error_message << std::endl << usage << std::endl;
                exit(1);
            } /* end of: if-then-else */
        } /* end-of: for */ 
    } /* end-of: if */
}
