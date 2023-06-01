//This is a header file for calculating the elements of the master formula
#ifndef HELICITY_AMPLITUDE_HPP
#define HELICITY_AMPLITUDE_HPP

#include "Kinematics.hpp"
#include "ThreeBodies_Kinematics.hpp"
#include "Form_Factors.hpp"
#include "Particle.hpp"
#include "Coupling_Constants.hpp"
#include "Functions.hpp"
#include <complex>

class Helicity_Amplitude
{
public:
   // k_charge is the electric charge of kaon

   std::complex<double> C1_1270( double s, double s13, double s23, int k_charge);
   std::complex<double> C2_1270( double s, double s13, double s23, int k_charge);

   std::complex<double> C1_1400( double s, double s13, double s23, int k_charge);
   std::complex<double> C2_1400( double s, double s13, double s23, int k_charge);

   std::complex<double> C3_1410( double s, double s13, double s23, int k_charge);
   
   std::complex<double> C3_1680( double s, double s13, double s23, int k_charge);
   
   std::complex<double> C1_1430( double s, double s13, double s23, double cos_theta, double phi1, double phi2, int k_charge);
   std::complex<double> C2_1430( double s, double s13, double s23, double cos_theta, double phi1, double phi2, int k_charge);
   std::complex<double> C3_1430( double s, double s13, double s23, double cos_theta, double phi1, double phi2, int k_charge);
   
   double J2_1270( double s, double s13, double s23, int k_charge);
   double J2_1400( double s, double s13, double s23, int k_charge);
   double J2_Total( double s, double s13, double s23, int k_charge, double T1K1270, double T1K1400);

   double JcrossJC_1270( double s, double s13, double s23, int k_charge);
   double JcrossJC_1400( double s, double s13, double s23, int k_charge);
   double JcrossJC_Total( double s, double s13, double s23, int k_charge, double T1K1270, double T1K1400);

private:
};

#endif
