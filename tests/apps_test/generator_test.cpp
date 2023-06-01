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
#include <vector>
#include <algorithm>


int main ( int argc, char* argv[] )
{
   //Define final particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_Pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   delete Pion;
   delete Kaon;

   //Define the charge of K_{res}
   int K_charge;
   std::cout << "What is the charge of the final Kaon? 0 for neutral, 1 for charged. (Integer only) \n>>>";
   std::cin >> K_charge;

   if( K_charge != 0 && K_charge != 1 )
   {
      std::cout << "Sorry. We are not ready to deal with charge " << K_charge << " kaons.\n";
   }

   assert( K_charge == 0 || K_charge == 1 ); 

   //Define the maximum S
   double sqrt_s_min = 2.0*mass_Pi + mass_K;

   double sqrt_s_max;
   std::cout << "What is the maximum centre-of-mass energy? (i.e. sqrt(s), not s.)\n";
   std::cout << "Caution : The minimum value is " << sqrt_s_min << " .\n>>>";
   std::cin >> sqrt_s_max; 
 
   double s_max = pow( sqrt_s_max, 2.0 );
   double s_min = pow( sqrt_s_min, 2.0 );

   if( s_max < s_min )
   {
      std::cout << "The maximum centre-of-mass energy is non-physical.\n";
   }
   assert( s_max >= s_min );
   
   int generation_mode;
   std::cout << "Constant centre-of-mass energy or random? 0 for constant, 1 for uniform random.\n>>>";
   std::cin >> generation_mode; 
   if( generation_mode != 0 && generation_mode != 1 )
   {
      std::cout << "Seems that you have entered something else? Bye.\n";
   }
   assert( generation_mode == 0 || generation_mode == 1 ); 
 
   //Define a Dalitz Region, and name it as PiPiK 
   Dalitz_Region PiPiK;
   PiPiK.SetSMaxM1M2M3(s_max, mass_Pi, mass_Pi, mass_K); //Setting the maximum S, and the mass of Pion, Pion, Kaon  

   double s12_min = PiPiK.Get_Sij_Min(1,2); //Get_Sij_Min(i,j)
   double s12_max = PiPiK.Get_Sij_Max(1,2);

   double s23_min = PiPiK.Get_Sij_Min(2,3);
   double s23_max = PiPiK.Get_Sij_Max(2,3);

   double s13_min = PiPiK.Get_Sij_Min(1,3);
   double s13_max = PiPiK.Get_Sij_Max(1,3);
   /*
   std::cout << "s13 min is " << s13_min << " .\n";
   std::cout << "s13 max is " << s13_max << " .\n";
   std::cout << "s23 min is " << s23_min << " .\n";
   std::cout << "s23 max is " << s23_max << " .\n";
   */

   //Setting the particle strength
   double T1K1270;
   double T1K1400;
   double T1K1410;
   double T1K1430;
   double T1K1680;
   std::cout << "This is a program including different Kaonic resonance mode, including K1(1270), K1(1400), K'*(1410), K2(1430), and K''*(1680).\n";
   std::cout << "The default values of the resonance strength for K1(1270) and K1(1400) are `1.36603' and `-0.366025' respectively.\n";
   std::cout << "Please enter the values which you want, ONE BY ONE. The value of K1(1270) is ...\n>>>";
   std::cin >> T1K1270; 
   std::cout << "The value of K1(1400) is ...\n>>>";
   std::cin >> T1K1400; 
   std::cout << "The value of K'*(1410) is ...\n>>>";
   std::cin >> T1K1410; 
   std::cout << "The value of K2(1430) is ...\n>>>";
   std::cin >> T1K1430; 
   std::cout << "The value of K''*(1680) is ...\n>>>";
   std::cin >> T1K1680; 

   //Define photon parameter value Lambda
   double lambda;  
   std::cout << "Please enter the photon polarization parameter (lambda).\n";
   std::cout << "It should be a value between -1 and 1. The more close to 1, the more right-handed polarized events.\n";
   std::cout << "What is the parameter value?\n>>>";
   std::cin >> lambda;

   if( (lambda < -1.) || (lambda > 1) )
   {
      std::cout << "The photon polarization parameter value is non-physical.\n";
   }
   assert( -1. <= lambda <= 1. );

   double c_RH = (1.+lambda)/2.; //And c_LH = (1-lambda)/2. , but we don't need it. 
   /*
   Both c_RH and c_LH here are numbers between 0 and 1. 
   That is we only care the relative strength between the right and left polarization.
   */

   //Define number of Events
   int N_events; 
   std::cout << "How many generation loops do you want to try? (NOT number of events.)\n>>>";
   std::cin >> N_events; 

   //Define vectors to safe the generated values
   std::vector<double> vec_s13; //a vector to save s13
   vec_s13.reserve(N_events); //the size of vector should be at most the sam as N_events
   std::vector<double> vec_s23;
   vec_s23.reserve(N_events);
   std::vector<double> vec_s_inv;
   vec_s_inv.reserve(N_events);

   std::vector<double> vec_CosTheta;
   vec_CosTheta.reserve(N_events);
   std::vector<double> vec_phi1;
   vec_phi1.reserve(N_events);

   std::vector<double> vec_ME2; //For calculating matrix Elements from vector J
   vec_ME2.reserve(N_events);

   //Define the function wanted 
   Matrix_Elements MEs; //Define an object "Matrix_Elements", and name it as ``MEs''.
   double max_ME2 = 0.0;
   int loop = 0;

   //Set up the random generator   
   RandomGen generate; //define object (a Random Number generator) by using RandomGen class, and name it as ``generate''.  
 

   //Find out the maximum
   for(int i=0; i < N_events; i++)
   {
      double s_inv;
      if(generation_mode == 1)
      {
         s_inv = generate.Uniform( s_min , s_max );
      }
      if(generation_mode == 0)
      {
         s_inv = s_max;
      }
      double s13 = generate.Uniform( s13_min , s13_max );
      double s23 = generate.Uniform( s23_min , s23_max );
 
      double cos_theta = generate.Uniform(-1., 1.);
      double phi1 = generate.Uniform(0.0, 2.*M_PI);
      
      double s23_Dmin = PiPiK.Get_S23_Dalitz_Min_wrt_S13(s_inv, s13);
      double s23_Dmax = PiPiK.Get_S23_Dalitz_Max_wrt_S13(s_inv, s13);

      double polarization; //either 1 or -1, should be an integer in theory.
      double polarize = generate.Uniform(0., 1.); 
      if( polarize < c_RH )
      {
         polarization = 1.;
      }
      else
      {
         polarization = -1.;

      }


      if ( remainder(i, 10000) == 0 )
      {
         std::cout << "Loop " << 10000*loop << " .\n"; //Print message for every 10000 generations
         loop++;  
      }


      if( s23_Dmin == 0.0 )//if non-physics value got,  s23_Dmin would be return as 0.0 from the function "Get_S23_Dalitz_Min_wrt_S13(s_inv,s13)"
      {
         continue;
      }


      if( (s23 >= s23_Dmin) && (s23 <= s23_Dmax) )
      {
         //The following values should be calculated within the if-statement(i.e. here), or it will trigger error because of the unphysical s23 value.

         std::complex<double> Jx = MEs.Jx_All_Resonances(s_inv, s13, s23, cos_theta, phi1, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
         std::complex<double> Jy = MEs.Jy_All_Resonances(s_inv, s13, s23, cos_theta, phi1, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
         //std::complex<double> Jz = MEs.Jz_Two_Resonances(s_inv, s13, s23, cos_theta, phi1, K_charge, 1.36603, -0.366025);     
         //Jz is multiplied by zero. That's why we don't calculate it here.
         std::complex<double> ME = (Jx  + Jy * (0., 1.) * polarization )/sqrt(2.); //Using `lambda' instead of `polarization' may lead to wrong result.
         //We use `lambda' to help calculating the percentage of events of right and left polarization.
         std::complex<double> ME2 = ME * std::conj(ME);
      
         double value_ME2 = std::real(ME2);
         if(std::imag(ME2) != 0.)
         {
            std::cout << "Non-zero imaginary part is found in J2.\n"; 
         }
 
         //Saving the values
         vec_s13.push_back(s13);
         vec_s23.push_back(s23);
         vec_s_inv.push_back(s_inv);
         vec_ME2.push_back(value_ME2);
         vec_CosTheta.push_back(cos_theta);
         vec_phi1.push_back(phi1);

         if (value_ME2 > max_ME2)
         {
            max_ME2 = value_ME2;
         }       

      }
   }

   std::cout << "The maximum value of the square of Matrix Element is " << max_ME2 << " .\n";
   std::cout << "There are " << vec_ME2.size() << " events within physical region.\n";
   std::cout << "Now start the hit-or-miss procedure.\n";

   //Prepare the output file
   std::ofstream s13_output("./../output/s13.dat");
   std::ofstream s23_output("./../output/s23.dat");
   std::ofstream s_output("./../output/s.dat");
   std::ofstream s13s23_output("./../output/s13s23.dat"); //For Mathematica use
   std::ofstream CosTheta_output("./../output/CosTheta.dat");
   std::ofstream phi1_output("./../output/phi1.dat");
   std::ofstream CosTheta_P_output("./../output/CosTheta_P.dat"); //only save non-negative values
   std::ofstream CosTheta_N_output("./../output/CosTheta_N.dat"); //only save negative values into positive values

   assert(s13_output.is_open());
   assert(s23_output.is_open());
   assert(s_output.is_open());
   assert(s13s23_output.is_open());
   assert(CosTheta_output.is_open());
   assert(phi1_output.is_open());
   assert(CosTheta_P_output.is_open());
   assert(CosTheta_N_output.is_open());

    
   //The same loop
   int Dalitz_count = 0;
   int size = vec_ME2.size();

   for(int i=0; i < size; i++)
   {
      double weight = vec_ME2[i];
      double u_w = generate.Uniform( 0.0, max_ME2 );     
      if( weight > u_w )//The generated point "u_w" is smaller than the function value "weight"
      { 
         Dalitz_count++;

         //saving the value
         s13_output << vec_s13[i] << "\n";
         s23_output << vec_s23[i] << "\n";
         s_output << vec_s_inv[i] << "\n";
         CosTheta_output << vec_CosTheta[i] << "\n";
         phi1_output << vec_phi1[i] << "\n";
         s13s23_output << vec_s13[i] << "," << vec_s23[i] << "\n"; //This format is mainly for Mathematica fast checking
         if( vec_CosTheta[i] >= 0.)
         {
            CosTheta_P_output << vec_s13[i] << " " << vec_s23[i] << " " << vec_CosTheta[i] << "\n"; //This format is for root
         }
         else if( vec_CosTheta[i] < 0.)
         {
            CosTheta_N_output << vec_s13[i] << " " << vec_s23[i] << " " << std::abs(vec_CosTheta[i]) << "\n";
         }
      }
   }

   std::cout << "There are " << Dalitz_count << " events left after hit-or-miss, out of " << N_events << " generations.\n";
   std::cout << "The generated numbers are saved in the \"./../output\" directory.\n";

   
   s13_output.close();
   s23_output.close();
   s_output.close();
   s13s23_output.close();
   CosTheta_output.close();
   phi1_output.close();
   CosTheta_P_output.close();
   CosTheta_N_output.close();

   return 0;

}
