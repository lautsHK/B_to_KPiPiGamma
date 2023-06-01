//This is a header file for generating random numbers

#ifndef RANDOMGEN_HPP
#define RANDOMGEN_HPP

class RandomGen 
{
public:

   //generate Uniform 
   static double Uniform( double min, double max );

   //generate Gaussian
   static double Gaussian( double mean, double SD );

   //generate Breit-Wigner
   static double BreitWigner( double mass, double width, double min, double max );

private:
};

#endif
