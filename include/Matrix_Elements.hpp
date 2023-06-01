//This is a header file for directly calculating the matrix elements
#ifndef MATRIX_ELEMENTS_HPP
#define MATRIX_ELEMENTS_HPP

#include "Kinematics.hpp"
#include "ThreeBodies_Kinematics.hpp"
#include "Form_Factors.hpp"
#include "Particle.hpp"
#include "Coupling_Constants.hpp"
#include "Functions.hpp"
#include <complex>

class Matrix_Elements
{
public:
   // k_charge is the electric charge of kaon

   std::complex<double> Jx_Two_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400);
   std::complex<double> Jy_Two_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400);
   std::complex<double> Jz_Two_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400);

   std::complex<double> Jx_All_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680);
   std::complex<double> Jy_All_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680);
   std::complex<double> Jz_All_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680);

   double J2( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680);

   double Im_N_dot_J_cross_Jc( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680);

   double ME2_LH( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680);

   double ME2_RH( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680);


private:
};

#endif
