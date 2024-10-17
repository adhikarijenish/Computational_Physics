#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>


std::vector<double> random_number( std::mt19937 & engine){
    std::uniform_real_distribution<double> dis(-1,1);

    std::vector<double> random_2D;
    random_2D.push_back(dis(engine));
    random_2D.push_back(dis(engine));

    return random_2D;
}

double mean(std::vector<double> &vector_1D){
    int length = vector_1D.size();
    double sum = 0;
    for(const double &s : vector_1D){
        sum += s;
    }
    return (double) sum/length;
} 