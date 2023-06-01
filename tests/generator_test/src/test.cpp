#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <random> //c++11 is needed
#include <map>

double Uniform( double min, double max);
double Gaussian( double mean, double SD);

 
int main ( int argc, char* argv[] )
{  
   std::cout << "Random number using <random> library" << ".\n";

   for(int i=0; i<10; i++)
   {
      double Random_number_11 = Uniform(0.,1.);
      std::cout << Random_number_11 << ".\n";
   }

   std::cout << "Random Gaussian number using <random> library" << ".\n";

   for(int i=0; i<10; i++)
   { 
      double trial_number = Gaussian(0.0, 1.0);
      std::cout << trial_number << ".\n";
   }
   

   return 0;

}


  
//Functions//

double Uniform( double min, double max )
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


double Gaussian( double mean, double SD)
{
    std::random_device rd;
    std::mt19937 mt_generator{rd()};
    std::normal_distribution<double> distribution(mean, SD);

    double random_number = distribution(mt_generator);
    return random_number;
}

