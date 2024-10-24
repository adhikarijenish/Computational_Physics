#include "../header/uebung1.hpp"
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>


std::string filename;

int main(int argc, char *argv[]){
    size_t seed = 123456789;
    std::random_device device; // Random Seed
    std::mt19937 engine(device()); // generate random number

     // Only Vectors //
    std::vector<double> random_pair;
    std::vector<double> radius;

    std::ofstream file;

    //Parameters//
    bool save;
    int Pairs;
    int X ;
    
    // Variables //
    double x_square;
    double y_square;
    double radii;
    double radii_square;
    double estamator;

    GetUserParam(argc,argv,filename,Pairs,X,save);

    if(save){
        file.open(filename);
    }

    for(int i =0;i<X;i++){
        for(int j=0;j<Pairs;j++){
            random_pair.clear();
            random_number(engine,random_pair);
            x_square = random_pair[0];
            x_square = x_square*x_square;
            y_square = random_pair[1];
            y_square = y_square *y_square;

            radii_square =x_square+y_square;
            radii = std::sqrt(radii_square);
            
            if(save){
                file<<radii<<" "<<radii_square<<"\n";
                }
            
            if (radii_square <=1){
                radius.push_back(1);
            }
            
        }
        estamator = 4.0*sum(radius)/ (double) Pairs;
        std::cout << estamator << "\n";
        radius.clear();

    
        if(file.is_open()){file.close();}

    }
    return 0;
}