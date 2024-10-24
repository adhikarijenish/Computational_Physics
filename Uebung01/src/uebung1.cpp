#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include <cstring>


void random_number(std::mt19937 & engine,std::vector<double> &random_2D){
    
    std::uniform_real_distribution<double> dis(-1.0,1.0);

    random_2D.push_back(dis(engine));
    random_2D.push_back(dis(engine));
}

double mean(std::vector<double> &vector_1D){
    int length = vector_1D.size();
    double sum = 0;
    for(const double &s : vector_1D){
        sum += s;
    }
    return (double) sum/length;
} 

double sum(std::vector<double> &vector_1D){
    double sum = 0;
    for(const double &s : vector_1D){
        sum += s;
    }
    return sum;
}

void GetUserParam( int argc, char *argv[],std::string & filename, int & Pairs, int & X, bool & save){

/* Variablen: */
    int i;
    char* endptr;
    const char usage[] = 
        "boundary [-f <File Name> -p <Pairs> -x <Experiment> -s <Save Radius>]";
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
                        save = argv[++i] == "True";
                        
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
