//include the header files
#include "RandomGen.hpp"
#include "Dalitz_Region.hpp"
#include "Functions.hpp"
#include "Kinematics.hpp"
#include "ThreeBodies_Kinematics.hpp"
#include "Coupling_Constants.hpp"
#include "Particle.hpp"
#include "Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp"
#include "Quark_Pair_Creation_Model/Kres_to_VP.hpp"
#include "Form_Factors.hpp"
#include "Helicity_Amplitude.hpp"
#include "Matrix_Elements.hpp"

//include the c++ libraries
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <complex>


int main ( int argc, char* argv[] )
{ 
   //Test 1, testing the std::abs and abs
   /*
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);
 
   double s = 2.;
   double s13 = 1.;
   double s23 = 1.;
   double theta = 0.5;
   double phi1 = 0.5;
   int k_charge = 1;

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi2 = acos( p1p2 / (std::abs(p1)*std::abs(p2)) ) + phi1;
 
   Helicity_Amplitude* HA = new Helicity_Amplitude;

   std::complex<double> C1;
   C1 = HA->C1_1270(s, s13, s23, k_charge);
   double C1_sq;
   C1_sq = pow( abs(C1), 2.0 );
   double C1_sq_std = pow( std::abs(C1), 2.0 );
   std::cout << "C1 is " << C1 << ".\n";
   std::cout << "C1^2 is " << C1_sq << ".\n";
   std::cout << "C1^2(using std::) is " << C1_sq_std << ".\n";

   std::complex<double> C2;
   C2 = HA->C2_1270(s, s13, s23, k_charge);
   double C2_sq;
   C2_sq = pow( abs(C2), 2.0 );
   std::cout << "C2 is " << C1 << ".\n";
   //std::cout << "C2^2 is " << C1_sq << ".\n";
   */

   //Test 2, testing the computation of J2
   /*
   std::complex<double> C1C2_z = C1 * std::conj(C2);
   std::complex<double> C2C1_z = C2 * std::conj(C1);
   std::cout << "C1*C2 is " << C1C2_z << ".\n";
   std::cout << "C2*C1 is " << C2C1_z << ".\n";

   double C1C2 = real(C1C2_z);
   double C2C1 = real(C2C1_z);
   
   double J2;
   J2 = C1_sq * p1_sq + C2_sq * p2_sq - (C1C2 + C2C1) * p1p2;  

   double J2_HA;
   J2_HA = HA->J2_1270(s, s13, s23, k_charge);
   

   double p1_x = std::abs(p1) * cos(theta) * cos(phi1);
   double p2_x = std::abs(p2) * cos(theta) * cos(phi2);
   std::cout << "p1_x is " << p1_x << ".\n";
   std::cout << "abs(p1) is " << std::abs(p1) << ".\n";
   std::cout << "cos(theta) is " << cos(theta) << ".\n";
   std::cout << "cos(phi1) is " << cos(phi1) << ".\n";

   std::cout << "p2_x is " << p2_x << ".\n";

   std::complex<double> Ji_2R;
   Ji_2R = C1 * p1_x - C2 * p2_x;
   std::cout << "Jx is " << Ji_2R << ".\n";

   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete HA;
 
   //std::cout << "J2 is " << J2 << ".\n";
   //std::cout << "J2_HA is " << J2_HA << ".\n";
   */

   //Test 3, testing the Breit wigner function
   /*
   Functions* FN = new Functions; 
   double BWR1 = FN->NormalizedBreitWigner(0.5,1.403,0.174,0.1,5.0);
   double BWR2 = FN->NormalizedBreitWigner(1.0,1.403,0.174,0.1,5.0);
   double BWR3 = FN->NormalizedBreitWigner(1.5,1.403,0.174,0.1,5.0);
   double BWR4 = FN->NormalizedBreitWigner(2.5,1.403,0.174,0.1,5.0);
   double BWR5 = FN->NormalizedBreitWigner(4.0,1.403,0.174,0.1,5.0);
   //The result is normalized
   //std::cout << BWR1 << "," << BWR2 << "," << BWR3 << "," << BWR4 << "," << BWR5 << ".\n";
   
   std::complex<double> BWRz1 = FN->BWPropagator(0.5,1.403,0.174);
   std::complex<double> BWRz2 = FN->BWPropagator(1.0,1.403,0.174);
   std::complex<double> BWRz3 = FN->BWPropagator(1.5,1.403,0.174);
   std::complex<double> BWRz4 = FN->BWPropagator(2.5,1.403,0.174);
   std::complex<double> BWRz5 = FN->BWPropagator(4.0,1.403,0.174);
   //The result is a complex number
   //std::cout << BWRz1 << "," << BWRz2 << "," << BWRz3 << "," << BWRz4 << "," << BWRz5 << ".\n";

   double BWR1p = std::abs(BWRz1)*std::abs(BWRz1);
   double BWR2p = std::abs(BWRz2)*std::abs(BWRz2);
   double BWR3p = std::abs(BWRz3)*std::abs(BWRz3);
   double BWR4p = std::abs(BWRz4)*std::abs(BWRz4);
   double BWR5p = std::abs(BWRz5)*std::abs(BWRz5);
   //The result is a real number.
   //std::cout << BWR1p << "," << BWR2p << "," << BWR3p << "," << BWR4p << "," << BWR5p << ".\n";
  
   double NR2 = BWR1/BWR2; 
   double NR3 = BWR1/BWR3; 
   double NR4 = BWR1/BWR4; 
   double NR5 = BWR1/BWR5; 
   std::cout << NR2 << "," << NR3 << "," << NR4 << "," << NR5 << ".\n";

   double R2 = BWR1p/BWR2p; 
   double R3 = BWR1p/BWR3p; 
   double R4 = BWR1p/BWR4p; 
   double R5 = BWR1p/BWR5p; 
   std::cout << R2 << "," << R3 << "," << R4 << "," << R5 << ".\n";

   delete FN; 
   */

   //Test 4, using I(0., 1.)
   //This is correct
   std::complex<double> Im_i(0.,1.);
   std::cout << "Imaginary I1 is " << Im_i << "\n";

   //This is wrong
   std::complex<double> sqrt_m1_1;
   sqrt_m1_1 = sqrt(-1.);
   std::cout << "Imaginary I2 is " << sqrt_m1_1 << "\n";

   //This is correct
   std::complex<double> sqrt_m1_2;
   std::complex<double> m1(-1.,0.);
   sqrt_m1_2 = sqrt(m1);
   std::cout << "Imaginary I3 is " << sqrt_m1_2 << "\n";

   std::complex<double> abs_2;
   abs_2 = Im_i * conj(Im_i);

   double trial;
   trial = pow(abs(Im_i), 2.);

   std::cout << "The result of i * i^*(1) is " << abs_2 << ".\n";
   std::cout << "The result of i * i^*(2) is " << trial << ".\n";
   
   std::complex<double> c(2.,3.);
   std::complex<double> d;
   d = c * Im_i;

   std::complex<double> e;
   e = c * (0., 1.);

   std::cout << "The result of (2.,3.) * I is " << d << ".\n";
   std::cout << "The result of (2.,3.) * (0., 1.) is " << e << ".\n";


   return 0;

}
