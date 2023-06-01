//including headers
#include "RandomGen.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <random> 

double RandomGen::Uniform( double min, double max )
{
    if ( min > max ) 
    {
        std::cout << "Minimum " << min << " is larger than maximum " << max << ".\n" ;
        abort();
    }

    std::random_device rd; //a generator with poor performance, used as a random seed
    std::mt19937 mt_generator{rd()}; //a generator with good performance
    std::uniform_real_distribution<double> distribution(min, max); //define a type of distribution

    double random_number = distribution(mt_generator);
    return random_number;
}


double RandomGen::Gaussian( double mean, double SD)
{
    std::random_device rd;
    std::mt19937 generator{rd()};
    std::normal_distribution<double> distribution(mean, SD);

    double random_number = distribution(generator);
    return random_number;
}




