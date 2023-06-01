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
   int K_charge = atoi(argv[1]); //read from command line in terminal

   assert( K_charge == 0 || K_charge == 1 ); 

   //Define the maximum S
   double sqrt_s_min = 2.0*mass_Pi + mass_K;

   double sqrt_s_max = atof(argv[2]); //read from command line in terminal

   double s_max = pow( sqrt_s_max, 2.0 );
   double s_min = pow( sqrt_s_min, 2.0 );

   if( s_max < s_min )
   {
      std::cout << "The maximum centre-of-mass energy is non-physical.\n";
   }
   assert( s_max >= s_min );

   //Define generation mode   
   int generation_mode = atoi(argv[3]); //read from command line in terminal

   assert( generation_mode == 0 || generation_mode == 1 ); 
 
   //Define a Dalitz Region, and name it as PiPiK 
   ThreeBodies_Kinematics PiPiK_Kine(mass_Pi, mass_Pi, mass_K);  

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
   double T1K1270 = atof(argv[4]); //read from command line in terminal
   double T1K1400 = atof(argv[5]); //read from command line in terminal
   double T1K1410 = atof(argv[6]); //read from command line in terminal
   double T1K1430 = atof(argv[7]); //read from command line in terminal
   double T1K1680 = atof(argv[8]); //read from command line in terminal

   //Define photon parameter value Lambda
   double lambda = atof(argv[9]); //read from command line in terminal

   assert( -1. <= lambda <= 1. );

   
   //c_RH is |c_{R}|^{2}, similar for L.
   double c_RH = (1.+lambda)/2.; //And c_LH = (1-lambda)/2. , but we don't need it. 
   /*
   Both c_RH and c_LH here are numbers between 0 and 1. 
   That is we only care the relative strength between the right and left polarisation.
   */

   //Define number of Events
   int N_events = atoi(argv[10]); //read from command line in terminal

   //Define the directory name
   std::string prifix = argv[11]; //read from command line in terminal

   //Define the function wanted 
   Matrix_Elements MEs; //Define an object "Matrix_Elements", and name it as ``MEs''.
   double max_ME2_RH = 0.0;
   double max_ME2_LH = 0.0;
   int phy_event_numbers = 0;

   //Set up the random generator   
   RandomGen generate; //define object (a Random Number generator) by using RandomGen class, and name it as ``generate''.  
 

   //Find out the maximum
   for(int i=0; i < 10000; i++)
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
      double phi = generate.Uniform(0.0, 2.*M_PI);
      
      double s23_Dmin = PiPiK.Get_S23_Dalitz_Min_wrt_S13(s_inv, s13);
      double s23_Dmax = PiPiK.Get_S23_Dalitz_Max_wrt_S13(s_inv, s13);

      double polarisation; //either 1 or -1, should be an integer in theory.
      double polarise = generate.Uniform(0., 1.); 
      if( polarise < c_RH )
      {
         polarisation = 1.;
      }
      else
      {
         polarisation = -1.;

      }


      if ( remainder(i, 10000) == 0 )
      {
         std::cout << "Checking the maximum value of the Matrix Element : Loop " << i << " .\n"; //Print message for every 10000 generations
      }


      if( s23_Dmin == 0.0 )//if non-physics value got,  s23_Dmin would be return as 0.0 from the function "Get_S23_Dalitz_Min_wrt_S13(s_inv,s13)"
      {
         continue;
      }


      if( (s23 >= s23_Dmin) && (s23 <= s23_Dmax) )
      {
         phy_event_numbers++;
         //The following values should be calculated within the if-statement(i.e. here), or it will trigger error because of the unphysical s23 value.

         std::complex<double> Jx = MEs.Jx_All_Resonances(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
         std::complex<double> Jy = MEs.Jy_All_Resonances(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
         //std::complex<double> Jz = MEs.Jz_All_Resonances(s_inv, s13, s23, cos_theta, phi, K_charge, 1.36603, -0.366025);     
         //Jz is multiplied by zero. That's why we don't calculate it here.
         std::complex<double> ME = (Jx  + Jy * (0., 1.) * polarisation )/sqrt(2.); //Using `lambda' instead of `polarisation' may lead to wrong result.
         //We use `lambda' to help calculating the percentage of events of right and left polarisation.
         std::complex<double> ME2 = ME * std::conj(ME);
      
         double value_ME2 = std::real(ME2);
         if(std::imag(ME2) != 0.)
         {
            std::cout << "Non-zero imaginary part is found in J2.\n"; //For checking purpose only. This senctence should never be printed. 
         }
 
         if ( (value_ME2 > max_ME2_RH) && (polarisation == 1.) )
         {
            max_ME2_RH = value_ME2;
         }       
         if ( (value_ME2 > max_ME2_LH) && (polarisation == -1.) )
         {
            max_ME2_LH = value_ME2;
         }       

      }
   }

   std::cout << "The maximum value of the square of R.H. Matrix Element is " << max_ME2_RH << " .\n";
   std::cout << "The maximum value of the square of L.H. Matrix Element is " << max_ME2_LH << " .\n";
   std::cout << "There are " << phy_event_numbers << " events within physical region.\n";
   std::cout << "Now start the hit-or-miss procedure.\n";

   //Prepare the output file
   std::string s13_file = "./../output/" + prifix + "_s13.dat";
   std::string s23_file = "./../output/" + prifix + "_s23.dat";
   std::string s_file = "./../output/" + prifix + "_s.dat";
   std::string s13s23_file = "./../output/" + prifix + "_s13s23.dat";
   std::string CosTheta_file = "./../output/" + prifix + "_CosTheta.dat";
   std::string Polar_file = "./../output/" + prifix + "_polar.dat";
   std::string phi_file = "./../output/" + prifix + "_phi.dat";
   std::string ME2_file = "./../output/" + prifix + "_ME2.dat";
   std::string Jx2_file = "./../output/" + prifix + "_Jx2.dat";
   std::string Jy2_file = "./../output/" + prifix + "_Jy2.dat";
   std::string Jz2_file = "./../output/" + prifix + "_Jz2.dat";


   std::ofstream s13_output(s13_file.c_str());
   std::ofstream s23_output(s23_file.c_str());
   std::ofstream s_output(s_file.c_str());
   std::ofstream s13s23_output(s13s23_file.c_str()); 
   std::ofstream CosTheta_output(CosTheta_file.c_str());
   std::ofstream Polar_output(Polar_file.c_str());
   std::ofstream phi_output(phi_file.c_str());
   std::ofstream ME2_output(ME2_file.c_str());
   std::ofstream Jx2_output(Jx2_file.c_str());
   std::ofstream Jy2_output(Jy2_file.c_str());
   std::ofstream Jz2_output(Jz2_file.c_str());
   /*
   std::ofstream CosTheta_P_output("./../output/CosTheta_P.dat"); //only save non-negative values
   std::ofstream CosTheta_N_output("./../output/CosTheta_N.dat"); //only save negative values into positive values
   */
   assert(s13_output.is_open());
   assert(s23_output.is_open());
   assert(s_output.is_open());
   assert(s13s23_output.is_open());
   assert(CosTheta_output.is_open());
   assert(Polar_output.is_open());
   assert(phi_output.is_open());
   assert(ME2_output.is_open());
   assert(Jx2_output.is_open());
   assert(Jy2_output.is_open());
   assert(Jz2_output.is_open());
   /*
   assert(CosTheta_P_output.is_open());
   assert(CosTheta_N_output.is_open());
   */
    
   //Event Generation
   int number_of_events = 0;
   while( number_of_events < N_events )
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
      //double cos_theta = -1.;
      double phi = generate.Uniform(0.0, 2.*M_PI);
      
      double s23_Dmin = PiPiK.Get_S23_Dalitz_Min_wrt_S13(s_inv, s13);
      double s23_Dmax = PiPiK.Get_S23_Dalitz_Max_wrt_S13(s_inv, s13);

      double polarisation; //either 1 or -1, should be an integer in theory.
      double polarise = generate.Uniform(0., 1.);
      double u_w; 
      if( polarise < c_RH ) //c_RH is calculated from lambda
      {
         polarisation = 1.;
         u_w = generate.Uniform( 0.0, max_ME2_RH );     
      }
      else
      {
         polarisation = -1.;
         u_w = generate.Uniform( 0.0, max_ME2_LH );     

      }

      if( s23_Dmin == 0.0 )//if non-physics value got,  s23_Dmin would be return as 0.0 from the function "Get_S23_Dalitz_Min_wrt_S13(s_inv,s13)"
      {
         continue;
      }


      if( (s23 >= s23_Dmin) && (s23 <= s23_Dmax) )
      {
         //The following values should be calculated within the if-statement(i.e. here), or it will trigger error because of the unphysical s23 value.

         std::complex<double> Jx = MEs.Jx_All_Resonances(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
         std::complex<double> Jy = MEs.Jy_All_Resonances(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
         std::complex<double> Jz = MEs.Jz_All_Resonances(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);    
         //Jz is multiplied by zero. That's why we don't calculate it here.
         std::complex<double> ME = (Jx  + Jy * (0., 1.) * polarisation )/sqrt(2.); //Using `lambda' instead of `polarisation' may lead to wrong result.
         //We use `lambda' to help calculating the percentage of events of right and left polarisation.
         std::complex<double> ME2 = ME * std::conj(ME);
         std::complex<double> Jx2 = Jx * std::conj(Jx);
         std::complex<double> Jy2 = Jy * std::conj(Jy);
         std::complex<double> Jz2 = Jz * std::conj(Jz);
     
         double value_Jx2 = std::real(Jx2);
         double value_Jy2 = std::real(Jy2);
         double value_Jz2 = std::real(Jz2);
 
         double weight_ME2 = std::real(ME2)*1.;

         if( weight_ME2 > u_w )//The generated point "u_w" is smaller than the function value "weight"
         { 
            number_of_events++;

            //saving the value
            s_output << s_inv << "\n";
            s13_output << s13 << "\n";
            s23_output << s23 << "\n";
            CosTheta_output << cos_theta << "\n";
            Polar_output << polarisation << "\n";
            phi_output << phi << "\n";
            s13s23_output << s13 << " " << s23 << "\n"; //This format is mainly for root
            ME2_output << weight_ME2 << "\n";
            Jx2_output << value_Jx2 << "\n";
            Jy2_output << value_Jy2 << "\n";
            Jz2_output << value_Jz2 << "\n";

            if ( remainder(number_of_events, 500) == 0 )
            {
               std::cout <<  number_of_events << " events are generated .\n"; //Print message for every 10000 generations
            }
         }
      }
   }


   std::cout << "The generated numbers are saved in the \"./../output\" directory.\n";

   
   s13_output.close();
   s23_output.close();
   s_output.close();
   s13s23_output.close();
   CosTheta_output.close();
   Polar_output.close();
   phi_output.close();
   ME2_output.close();
   Jx2_output.close();
   Jy2_output.close();
   Jz2_output.close();

   return 0;

}
