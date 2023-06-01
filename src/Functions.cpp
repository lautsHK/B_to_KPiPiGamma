//including headers
#include "Functions.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <complex>

double Functions::Gaussian( double x, double mean, double SD )
{
   double normalize = SD*sqrt(2.0*M_PI);
   double exponent = -pow( (x-mean) , 2.0 ) / ( 2.0*pow( SD , 2.0 ) );

   double gauss;
   gauss = exp(exponent)/normalize;
   return gauss;
}
   
double Functions::NormalizedBreitWigner( double s, double mass, double width, double min, double max )
{
   double M2 = pow( mass, 2.0);
   double MTau = mass * width;
   
   double K;
   K = MTau / ( atan((max-M2)/MTau) - atan((min-M2)/MTau) );

   double BW;
   BW = K/( pow((s-M2),2.0) + pow(MTau,2.0) );
   return BW;
}



std::complex<double> Functions::BWPropagator( double s, double mass, double width )
{
   double M2 = pow( mass, 2.0 );
   double MTau = mass * width;
   std::complex<double> i(0.0,1.0);
   
   std::complex<double> BWP;
   BWP = 1.0 / (s - M2 + i*(MTau));
   return BWP;
}

