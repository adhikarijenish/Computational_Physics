#include "../header/uebung1.hpp"
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>

int main(){
    size_t seed = 123456789;
    std::random_device device; // Random Seed
    std::mt19937 engine(seed); // generate random number
    std::vector<double> random_pair;
    std::ofstream file;
    file.open("../data/oneBigexperiment.txt");



    double x_square;
    double y_square;
    double radii;
    double radii_square;

    int Pairs = 10000;
    int X = 1;

    for(int i =0;i<X;i++){
        for(int j=0;j<Pairs;j++){

            random_pair = random_number(engine);
            x_square = random_pair[0];
            x_square = x_square*x_square;
            y_square = random_pair[1];
            y_square = y_square *y_square;

            radii = std::sqrt(x_square+y_square);
            radii_square = radii *radii;
            file<<radii<<" "<< radii_square<<"\n";
            
        }
    }


    file.close();
    random_pair.clear();
    return 0;
}