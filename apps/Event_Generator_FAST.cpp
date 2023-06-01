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
#include "Transform.hpp"

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

   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();

   Particle* B_meson = new Particle(511); //511 is the PDGID of B0
   double mass_B = B_meson->Get_Particle_Mass();

   delete Pion;
   delete Kaon;
   delete KStar;
   delete B_meson;

   //Define the charge of K_{res}
   int K_charge = atoi(argv[1]); //read from command line in terminal

   assert( K_charge == 0 || K_charge == 1 ); 

   //Define the maximum S
   double sqrt_s_min_physical = 2.0*mass_Pi + mass_K;

   double sqrt_s_min = atof(argv[2]); //read from command line in terminal

   double sqrt_s_max = atof(argv[3]); //read from command line in terminal

   double s_max = pow( sqrt_s_max, 2.0 );
   double s_min = pow( sqrt_s_min, 2.0 );
   double s_min_physical = pow( sqrt_s_min_physical, 2.0 );

   if( s_max < s_min )
   {
      std::cout << "The maximum centre-of-mass energy is non-physical.\n";
   }
   assert( s_max >= s_min );
   assert( s_min >= s_min_physical );

   //Define generation mode   
   int generation_mode = atoi(argv[4]); //read from command line in terminal

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

   //Define a function can help us to calculate Kinematics variables, such as momentum.
   ThreeBodies_Kinematics PiPiK_Kine(mass_Pi, mass_Pi, mass_K);

   double p1[4]; // 0, 1, 2, 3 are E, px, py, pz
   double p2[4];
   double p3[4];
   double photon[4]; //Always z-direction
   double normal_z[3] = {0., 0., 1.}; 

   typedef double (ThreeBodies_Kinematics::*Momentum_Function)(double, double, double, double, double);

   Momentum_Function P1_entries[3] = {&ThreeBodies_Kinematics::Get_CoP_P1x, &ThreeBodies_Kinematics::Get_CoP_P1y, &ThreeBodies_Kinematics::Get_CoP_P1z};
   Momentum_Function P2_entries[3] = {&ThreeBodies_Kinematics::Get_CoP_P2x, &ThreeBodies_Kinematics::Get_CoP_P2y, &ThreeBodies_Kinematics::Get_CoP_P2z};
   Momentum_Function P3_entries[3] = {&ThreeBodies_Kinematics::Get_CoP_P3x, &ThreeBodies_Kinematics::Get_CoP_P3y, &ThreeBodies_Kinematics::Get_CoP_P3z};

   //Define the transformation for Lorentz boost and Rotation. Call the "Transform" class as "Apply".
   Transform Apply;

   double p1_Lab[4]; //0, 1, 2, 3, are E, px, py, pz
   double p2_Lab[4];
   double p3_Lab[4];
   double photon_Lab[4];

   typedef double (Transform::*Rotating)(double, double, double, double, double, double);
  
   Rotating Euler_rotate[3] = {&Transform::Euler_Rotation_x, &Transform::Euler_Rotation_y, &Transform::Euler_Rotation_z};

   //Setting the particle strength
   double T1K1270 = atof(argv[5]); //read from command line in terminal
   double T1K1400 = atof(argv[6]); //read from command line in terminal
   double T1K1410 = atof(argv[7]); //read from command line in terminal
   double T1K1430 = atof(argv[8]); //read from command line in terminal
   double T1K1680 = atof(argv[9]); //read from command line in terminal

   //Define photon parameter value Lambda
   double lambda = atof(argv[10]); //read from command line in terminal

   assert( -1. <= lambda <= 1. );

   
   //c_RH is |c_{R}|^{2}, similar for L.
   double c_RH = (1.+lambda)/2.; 
   double c_LH = (1.-lambda)/2.; 
   /*
   Both c_RH and c_LH here are numbers between 0 and 1. 
   This is because we only care the relative strength between the right and left polarisation.
   */

   //Define number of Events
   int Total_events = atoi(argv[11]); //read from command line in terminal

   //Define the directory name
   std::string prifix = argv[12]; //read from command line in terminal

   //Set up the random generator   
   RandomGen generate; //define object (a Random Number generator) by using RandomGen class, and name it as ``generate''.  
   
   //Decide how many left-handed and right-handed events
   int polarisation[2] = {1, -1};
   int N_events[2] = {0, 0};
   double polar;
   for(int i=0; i < Total_events; i++ )
   {
      polar = generate.Uniform(0. ,1.);
    
      if( polar < c_RH )
      {
         N_events[0]++;
      }
      else
      {
         N_events[1]++;
      }
   }
 
   std::cout << "There are " << N_events[0] << " RH events going to be generated.\n";
   std::cout << "There are " << N_events[1] << " LH events going to be generated.\n";
   
   //Define the function wanted 
   Matrix_Elements MEs; //Define an object "Matrix_Elements", and name it as ``MEs''.

   std::cout << "Now start the hit-or-miss procedure.\n";

   //Prepare the output file (variables randomlty generated in K_res frame)
   std::vector<std::string> filenames = 
   {
      "./../output/" + prifix + "_s.dat",
      "./../output/" + prifix + "_SqrtS.dat",
      "./../output/" + prifix + "_s13.dat",
      "./../output/" + prifix + "_s23.dat",
      "./../output/" + prifix + "_s13s23.dat",
      "./../output/" + prifix + "_CosTheta.dat",
      "./../output/" + prifix + "_phi.dat",
      "./../output/" + prifix + "_ME2.dat",
      "./../output/" + prifix + "_polar.dat",
      "./../output/" + prifix + "_nJXJc.dat",
      "./../output/" + prifix + "_J2.dat",
      "./../output/" + prifix + "_Jratio.dat"
   };

   std::vector<std::ofstream> outputs(filenames.size());

   for(int i=0; i<filenames.size(); i++) 
   {
      outputs[i].open(filenames[i].c_str());
      assert(outputs[i].is_open());
   }

   //Prepare the output file for the K_res frame particles
   std::vector<std::string> P_CoP_filenames = 
   {
      "./../output/" + prifix + "_photon_CoP.dat",
      "./../output/" + prifix + "_P1_CoP.dat",
      "./../output/" + prifix + "_P2_CoP.dat",
      "./../output/" + prifix + "_P3_CoP.dat",
   };

   std::vector<std::ofstream> Print_4_Momentum_CoP(P_CoP_filenames.size());

   for(int i=0; i<P_CoP_filenames.size(); i++) 
   {
      Print_4_Momentum_CoP[i].open(P_CoP_filenames[i].c_str());
      assert(Print_4_Momentum_CoP[i].is_open());
   }

   //Prepare the output file for the final particles
   std::vector<std::string> P_filenames = 
   {
      "./../output/" + prifix + "_photon_Lab.dat",
      "./../output/" + prifix + "_P1_Lab.dat",
      "./../output/" + prifix + "_P2_Lab.dat",
      "./../output/" + prifix + "_P3_Lab.dat",
   };

   std::vector<std::ofstream> Print_4_Momentum(P_filenames.size());

   for(int i=0; i<P_filenames.size(); i++) 
   {
      Print_4_Momentum[i].open(P_filenames[i].c_str());
      assert(Print_4_Momentum[i].is_open());
   }

   //Event Generation
   for(int side = 0; side < 2; side++)
   {
      int number_of_events = 0;
      while( number_of_events < N_events[side] )
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

         //double u_w = generate.Uniform( 0.0, max_ME2[side] );     

         //randomly generate Euler Rotation angles
         double angle_1 = generate.Uniform(0.0, 2.*M_PI);
         double Cos_angle_2 = generate.Uniform(-1., 1.);
         double angle_3 = generate.Uniform(0.0, 2.*M_PI);

         if( s23_Dmin == 0.0 )//if non-physics value got,  s23_Dmin would be return as 0.0 from the function "Get_S23_Dalitz_Min_wrt_S13(s_inv,s13)"
         {
            continue;
         }


         if( (s23 >= s23_Dmin) && (s23 <= s23_Dmax) )
         {
            //The following values should be calculated within the if-statement(i.e. here), or it will trigger error because of the unphysical s23 value.

            double ME2;
            if(side == 0)
            { 
               ME2 = MEs.ME2_RH(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680); 
            }
            else if(side == 1)
            {
               ME2 = MEs.ME2_LH(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680); 
            }

            double weight_ME2 = ME2*1.; //Actually we don't need to multiply by 1., it is reserved here for statistical reason

            number_of_events++;

            double sqrt_s = sqrt(s_inv);
            
            //saving the value
            outputs[0] << s_inv << "\n";
            outputs[1] << sqrt_s << "\n";
            outputs[2] << s13 << "\n";
            outputs[3] << s23 << "\n";
            outputs[4] << s13 << " " << s23 << "\n";
            outputs[5] << cos_theta << "\n";
            outputs[6] << phi << "\n";
            outputs[7] << weight_ME2 << "\n";
            outputs[8] << polarisation[side] << "\n";
 
            //calculating J2 and Im(n.JxJ*)
            double Im_nJXJc = MEs.Im_N_dot_J_cross_Jc(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);
            outputs[9] << Im_nJXJc << "\n"; 

            double value_J2 = MEs.J2(s_inv, s13, s23, cos_theta, phi, K_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);
            outputs[10] << value_J2 << "\n"; 

            double J_ratio = Im_nJXJc / value_J2;
            outputs[11] << J_ratio << "\n.";

            //calculating the final particles
            p1[0] = PiPiK_Kine.Get_CoP_Energy1(s_inv, s13, s23);
            p2[0] = PiPiK_Kine.Get_CoP_Energy2(s_inv, s13, s23);
            p3[0] = PiPiK_Kine.Get_CoP_Energy3(s_inv, s13, s23);
            double photon_E_CoP = (pow(mass_B, 2.0) - s_inv) / (2. * sqrt_s ); //Lab-frame
            photon[0] = photon_E_CoP;

            for(int i=0; i<3; i++) 
            {
               p1[i+1] = (PiPiK_Kine.*P1_entries[i])(s_inv, s13, s23, cos_theta, phi);
               p2[i+1] = (PiPiK_Kine.*P2_entries[i])(s_inv, s13, s23, cos_theta, phi);
               p3[i+1] = (PiPiK_Kine.*P3_entries[i])(s_inv, s13, s23, cos_theta, phi);
               photon[i+1] = normal_z[i] * photon_E_CoP ;
            }
            /*
            A better way is to use ROOT Lorentz vector. So by defining "TLorentzVector particle1_p(p1_x, p1_y, p1_z, E1);", people can use Lorentz boost and rotation already.
            But since here is the only part to use 4-vector in the generator, direct calculation is used here. 
            */ 
            Print_4_Momentum_CoP[0] << photon[0] << " " << photon[1] << " " << photon[2] << " " << photon[3] << " " << "\n";
            Print_4_Momentum_CoP[1] << p1[0] << " " << p1[1] << " " << p1[2] << " " << p1[3] << "\n";
            Print_4_Momentum_CoP[2] << p2[0] << " " << p2[1] << " " << p2[2] << " " << p2[3] << "\n";
            Print_4_Momentum_CoP[3] << p3[0] << " " << p3[1] << " " << p3[2] << " " << p3[3] << "\n";
     
            //Boosting and Rotating
            double beta = (pow(mass_B, 2.0) - s_inv) / (pow(mass_B, 2.0) + s_inv);

            p1_Lab[0] = Apply.Lorentz_Boost_E(beta, p1[0], p1[3]);
            p2_Lab[0] = Apply.Lorentz_Boost_E(beta, p2[0], p2[3]);
            p3_Lab[0] = Apply.Lorentz_Boost_E(beta, p3[0], p3[3]);
            photon_Lab[0] = Apply.Lorentz_Boost_E(beta, photon[0], photon[3]);

            //Only the z-direction is needed to be boosted becuase the photon is defined as the z-direction
            p1[3] = Apply.Lorentz_Boost_p(beta, p1[0], p1[3]); //p[3] is no longer the same as before. 
            p2[3] = Apply.Lorentz_Boost_p(beta, p2[0], p2[3]);
            p3[3] = Apply.Lorentz_Boost_p(beta, p3[0], p3[3]);
            photon[3] = Apply.Lorentz_Boost_p(beta, photon[0], photon[3]);

            for(int i=0; i<3; i++) //start from 1
            {
               p1_Lab[i+1] = (Apply.*Euler_rotate[i])(angle_1, Cos_angle_2, angle_3, p1[1], p1[2], p1[3]);
               p2_Lab[i+1] = (Apply.*Euler_rotate[i])(angle_1, Cos_angle_2, angle_3, p2[1], p2[2], p2[3]);
               p3_Lab[i+1] = (Apply.*Euler_rotate[i])(angle_1, Cos_angle_2, angle_3, p3[1], p3[2], p3[3]);
               photon_Lab[i+1] = (Apply.*Euler_rotate[i])(angle_1, Cos_angle_2, angle_3, photon[1], photon[2], photon[3]);
            }

            
            Print_4_Momentum[0] << photon_Lab[0] << " " << photon_Lab[1] << " " << photon_Lab[2] << " " << photon_Lab[3] << "\n";
            Print_4_Momentum[1] << p1_Lab[0] << " " << p1_Lab[1] << " " << p1_Lab[2] << " " << p1_Lab[3] << "\n";
            Print_4_Momentum[2] << p2_Lab[0] << " " << p2_Lab[1] << " " << p2_Lab[2] << " " << p2_Lab[3] << "\n";
            Print_4_Momentum[3] << p3_Lab[0] << " " << p3_Lab[1] << " " << p3_Lab[2] << " " << p3_Lab[3] << "\n";
            
            //printing loops
            if ( remainder(number_of_events, 500) == 0 )
            {
               std::cout <<  number_of_events << " events are generated .\n"; //Print message for every 10000 generations
            }


         }//End the s23(s13) max/min check
      }//End a loop of Left or Right
   }//End Left-Right

   std::cout << "The generated numbers are saved in the \"./../output\" directory.\n";

   for(int i=0; i<outputs.size(); i++) 
   {
      outputs[i].close();
   } 

   for(int i=0; i<Print_4_Momentum_CoP.size(); i++) 
   {
      Print_4_Momentum_CoP[i].close();
   } 

   for(int i=0; i<Print_4_Momentum.size(); i++) 
   {
      Print_4_Momentum[i].close();
   } 

   return 0;

}
