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


double RandomGen::BreitWigner( double mass, double width, double min, double max )
{
   //simplify the life a bit
   double M2 = pow( mass, 2.0);
   double MTau = mass * width;
   
   double K;
   K = MTau / ( atan((max-M2)/MTau) - atan((min-M2)/MTau) );

   std::random_device rd; //a generator with poor performance, used as a random seed
   std::mt19937 mt_generator{rd()}; //a generator with good performance
   std::uniform_real_distribution<double> distribution(0.0, 1.0); //define a type of distribution

   double u = distribution(mt_generator);
   double random_number_s;
   
   random_number_s = M2 + MTau * tan( u*MTau/K + atan((min-M2)/MTau) );

   return random_number_s;
}


   


