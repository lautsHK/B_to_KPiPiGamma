//This is a header file for storing some common functions

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <complex>

class Functions
{
public:

   double Gaussian( double x, double mean, double SD );

   double NormalizedBreitWigner( double s, double mass, double width, double min, double max );
   
   std::complex<double> BWPropagator( double s, double mass, double width );

private:
};

#endif
