#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include <cstring>


#include<random>

void random_number_generator(std::vector<double> &random_number,std::mt19937 & engine, const int & len){
    std::uniform_real_distribution<> dis(0, 1); // distribute uniformly in interval 0 and 1
    for(int i = 0; i<=len;i++){
        random_number.push_back(dis(engine)); // creat array of uniformly dis. random number
    }

}


void random_2D_Spin(std::vector<std::vector<int>>& Spin_2D,std::mt19937 & engine,const int & N){
    std::uniform_real_distribution<> dis(0, 1);
     for (int i = 0; i <N; i++) 
     {
          for (int j = 0; j<N; j++)
          {
            if(dis(engine)>= 0.5){
                Spin_2D[i][j] = 1;
            }
            else{
                Spin_2D[i][j] = -1;
            }
               
          }
     }
}

void random_3D_Spin(std::vector<std::vector<std::vector<int>>> & Spin_3D,std::mt19937 & engine,const int & N){
    std::uniform_real_distribution<> dis(0, 1);
     for (int i = 0; i <N; i++) 
     {
          for (int j = 0; j<N; j++)
          {
            for(int k = 0; k< N;k++){
                if(dis(engine)>= 0.5){
                Spin_3D[i][j][k] = 1;
            }
               else{
                Spin_3D[i][j][k] = -1;
            }
            }
               
          }
     }
}

double Energy_sum_2D(const std::vector<std::vector<int>>& Spin_config){
        double Energy = 0;
        int N = Spin_config.size();
        for(int i = 0 ; i<N;i++){
            for(int j = 0 ; j <N;j++){
                Energy += Spin_config[i][j]*(Spin_config[(i-1+N)%N][j] + Spin_config[(i+1+N)%N][j] + Spin_config[i][(j-1+N)%N]+ Spin_config[i][(j+1+N)%N]);
            }
        }
        return -0.5*Energy;
    }

double Energy_sum_3D(const std::vector<std::vector<std::vector<int>>> & Spin_3D){
        double Energy = 0;
        int N = Spin_3D.size();
        for(int i = 0 ; i<N;i++){
            for(int j = 0 ; j <N;j++){
                for( int k = 0; k<N ; k++){
                     Energy += 
                     Spin_3D[i][j][k]*(Spin_3D[(i-1+N)%N][j][k] 
                     + Spin_3D[(i+1+N)%N][j][k] + 
                     Spin_3D[i][(j-1+N)%N][k]+ 
                     Spin_3D[i][(j+1+N)%N][k]+
                     Spin_3D[i][j][(k+1+N)%N]+
                     Spin_3D[i][j][(k-1+N)%N]);

                }
            }
        }
        return -0.5*Energy;
    }

int delta_E_2D(const std::vector<std::vector<int>>& Spin_2D,const int & i,const int & j){
    int N = Spin_2D.size();
    int dE;
    dE = Spin_2D[i][j]*(Spin_2D[i][(j-1+N)%N]+Spin_2D[(i+1+N)%N][j]+Spin_2D[i][(j+1)%N]+Spin_2D[(i-1+N)%N][j]);
    return 2*dE;
}

int delta_E_3D(const std::vector<std::vector<std::vector<int>>> & Spin_3D,const int & i,const int & j,const int & k){
    int N = Spin_3D.size();
    int dE;
    dE = Spin_3D[i][j][k]*(Spin_3D[(i-1+N)%N][j][k] 
            + Spin_3D[(i+1+N)%N][j][k] + 
            Spin_3D[i][(j-1+N)%N][k]+ 
            Spin_3D[i][(j+1+N)%N][k]+
            Spin_3D[i][j][(k+1+N)%N]+
            Spin_3D[i][j][(k-1+N)%N]);
    return 2*dE;
}

int sum_vec(const std::vector<int>& Spin_1D){
    int sum = 0;
    for(const int &s :Spin_1D){
        sum +=s;
    }
    return sum;
}

double sum_vec(const std::vector<double>& Spin_1D){
    double sum = 0;
    for(const double &s :Spin_1D){
        sum +=s;
    }
    return sum;
}

double mean(const std::vector<double>& data){
    double S = 0;
    int N = data.size();
    for(const double & in :data ){
        S +=in;
    }
    return (double) S/N;
}

double mean(const std::vector<int>& data){
    int S = 0;
    int N = data.size();
    for(const int & in :data ){
        S +=in;
    }
    return (double) S/N;
}

double mean(const std::vector<std::vector<int>> & data){
    int N = data.size();
    
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            sum += data[i][j];
        }
    }
    N = N*N;
    return (double) sum/N;
}

double mean(const std::vector<std::vector<std::vector<int>>> & data){
    int sum = 0;
    int N = data.size();
    
    for(const auto & spin: data){
        for(const auto & array :spin){
            for(const int & in : array){
                sum += in;
            }
        }
    }
    N = N*N*N;
    return (double) sum/N;
}

double meansquare(const std::vector<double>& data){
    double S = 0.0;
    int N = data.size();
    for(const double & in:data){
        S +=in*in;
    }
    return  (double) S/N;
}

double variance(const std::vector<double>& data){
    double var;
    double temp = mean(data);
    var = meansquare(data) - temp*temp;
    return var;

}

int vector_2D_sum(const std::vector<std::vector<int>> & Spin_2D,std::vector<float> & z){
    int sum = 0;
    int N = Spin_2D.size();
    for(int i = 0;i<N;i++){
        for(int j = 0; j<N;j++){
            sum += Spin_2D[i][j];
            z.at(N*i+j) = Spin_2D[i][j];
        }
    }
    return sum;
}

void Ising_2D_Sweep(std::vector<std::vector<int>>& Spin_2D,std::mt19937 &engine,const double & B,double& E){
    int N = Spin_2D.size();
    std::uniform_real_distribution<double> real_dis(0.,1.);
    int dE;
    
    for(int i =0;i<N;i++){
        for(int j = 0;j<N;j++){            
        dE = delta_E_2D(Spin_2D,i,j);
        if(dE<=0 || real_dis(engine)<std::exp(-B*dE))
        {
            Spin_2D[i][j] *= -1;
            E +=dE;
            } 
        }
        

    }
    
}

void Ising_3D_Sweep(std::vector<std::vector<std::vector<int>>> & Spin_3D,std::mt19937 &engine,const double & B,double& E){
    int N = Spin_3D.size();
    std::uniform_real_distribution<double> real_dis(0.,1.);
    int dE;
    
    for(int i =0;i<N;i++){
        for(int j = 0;j<N;j++){ 
            for(int k = 0;k<N;k++){           
                dE = delta_E_3D(Spin_3D,i,j,k);
                if(dE<=0 || real_dis(engine)<std::exp(-B*dE))
                {
                    Spin_3D[i][j][k] *= -1;
                    E +=dE;
                    } 
                }
        

        }   
    
    }
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
