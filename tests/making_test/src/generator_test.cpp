//include the header files
#include "RandomGen.hpp"

//include the c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>


int main ( int argc, char* argv[] )
{
   RandomGen* generate = new RandomGen; //define object by using RandomGen class, and name is as ``generate''  
  
   std::cout << "Random number using <random> library" << ".\n";

   for(int i=0; i<10; i++)
   {
      double random_number = generate->Uniform(0.,1.); //function Uniform(min, max)
      std::cout << random_number << ".\n";
   }

   std::cout << "Random Gaussian number using <random> library" << ".\n";

   for(int i=0; i<10; i++)
   { 
      double random_gaussian = generate->Gaussian(0.0, 1.0); //function Gaussian(mean, SD)
      std::cout << random_gaussian << ".\n";
   }
   
   delete generate; //every ``new'' has a ``delete''

   return 0;

}
