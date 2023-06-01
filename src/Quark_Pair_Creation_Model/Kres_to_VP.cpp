//including headers
#include "Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp"
#include "Quark_Pair_Creation_Model/Kres_to_VP.hpp"
#include "Kinematics.hpp"
#include "Particle.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>

double Kres_to_VP::Get_theta_K() const
{
   double th;
   th = theta_K;

   return th;
}

double Kres_to_VP::MatrixElement_S_K1270(double mass, double mass_V, int decay)
{
   //initilize the value
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //double R = K1270->Get_R_SHO();  

   double mass_P[2]; //2 decay modes
   mass_P[0] = mass_pi; //K1 to K0st + Pi
   mass_P[1] = mass_K;  //K1 to rho + K


   QPCM_A_to_VplusP* K1270_S = new QPCM_A_to_VplusP;
   double MS = K1270_S->MSpqc(mass, mass_V, mass_P[decay]);

   //For M to AB
   Kinematics* kine = new Kinematics;
   double E_A = kine->M2ab_CoP_Energy_a(mass, mass_V, mass_P[decay]);
   double E_B = kine->M2ab_CoP_Energy_a(mass, mass_P[decay], mass_V);
   double C;
   C = 8. * pow(M_PI, 3./2.) * sqrt(E_A * E_B * mass);

   // If decay = 0, then it is K1 -> K* + \pi decay
   // For K1 = 1270, and decay = 0,
   // the rotation factor is :
   //    Sqrt(2) * ( Sin[\theta] - Cos[\theta] ).
   // While for K1 -> rho + K decay, i.e. decay = 1,  
   // the rotation factor is :
   //    Sqrt(2) * ( Sin[\theta] + Cos[\theta] ). 
   double R;
   R = sqrt(2.)*sin(theta_K) + pow(-1., 1. + decay) * cos(theta_K); //theta_K is assigned in the header file
   

   double MS_R;
   MS_R = C * MS * R;

   delete Pion;
   delete Kaon;

   delete K1270_S;

   return MS_R;
} 

double Kres_to_VP::MatrixElement_S_K1400(double mass, double mass_V, int decay)
{
   //initilize the value
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //double R = K1400->Get_R_SHO();  

   double mass_P[2]; //2 decay modes
   mass_P[0] = mass_pi; //K1 to K0st + Pi
   mass_P[1] = mass_K;  //K1 to rho + K


   QPCM_A_to_VplusP* K1400_S = new QPCM_A_to_VplusP;
   double MS = K1400_S->MSpqc(mass, mass_V, mass_P[decay]);

   //For M to AB
   Kinematics* kine = new Kinematics;
   double E_A = kine->M2ab_CoP_Energy_a(mass, mass_V, mass_P[decay]);
   double E_B = kine->M2ab_CoP_Energy_a(mass, mass_P[decay], mass_V);
   double C;
   C = 8. * pow(M_PI, 3./2.) * sqrt(E_A * E_B * mass);

   // If decay = 0, then it is K1 -> K* + \pi decay
   // For K1 = 1400, and decay = 0,
   // the rotation factor is :
   //    Sqrt(2) * ( Cos[\theta] + Sin[\theta] ).
   // While for K1 -> rho + K decay, i.e. decay = 1,  
   // the rotation factor is :
   //    Sqrt(2) * ( Cos[\theta] - Sin[\theta] ). 
   double R;
   R = sqrt(2.)*cos(theta_K) + pow(-1., 0. + decay) * sin(theta_K);

   double MS_R;
   MS_R = C * MS * R;

   delete Pion;
   delete Kaon;

   delete K1400_S;

   return MS_R;
} 


double Kres_to_VP::MatrixElement_D_K1270(double mass, double mass_V, int decay)
{
   //initilize the value
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //double R = K1270->Get_R_SHO();  

   double mass_P[2]; //2 decay modes
   mass_P[0] = mass_pi; //K1 to K0st + Pi
   mass_P[1] = mass_K;  //K1 to rho + K


   QPCM_A_to_VplusP* K1270_D = new QPCM_A_to_VplusP;
   double MD = K1270_D->MDpqc(mass, mass_V, mass_P[decay]);

   //For M to AB
   Kinematics* kine = new Kinematics;
   double E_A = kine->M2ab_CoP_Energy_a(mass, mass_V, mass_P[decay]);
   double E_B = kine->M2ab_CoP_Energy_a(mass, mass_P[decay], mass_V);
   double C;
   C = 8. * pow(M_PI, 3./2.) * sqrt(E_A * E_B * mass);

   double R;
   R = -1.*sin(theta_K) - pow(-1., 0. + decay) * sqrt(2.) * cos(theta_K);

   double MD_R;
   MD_R = C * MD * R;

   delete Pion;
   delete Kaon;

   delete K1270_D;

   return MD_R;
}


double Kres_to_VP::MatrixElement_D_K1400(double mass, double mass_V, int decay)
{
   //initilize the value
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //double R = K1270->Get_R_SHO();  

   double mass_P[2]; //2 decay modes
   mass_P[0] = mass_pi; //K1 to K0st + Pi
   mass_P[1] = mass_K;  //K1 to rho + K


   QPCM_A_to_VplusP* K1400_D = new QPCM_A_to_VplusP;
   double MD = K1400_D->MDpqc(mass, mass_V, mass_P[decay]);

   //For M to AB
   Kinematics* kine = new Kinematics;
   double E_A = kine->M2ab_CoP_Energy_a(mass, mass_V, mass_P[decay]);
   double E_B = kine->M2ab_CoP_Energy_a(mass, mass_P[decay], mass_V);
   double C;
   C = 8. * pow(M_PI, 3./2.) * sqrt(E_A * E_B * mass);

   double R;
   R = -1.*cos(theta_K) - pow(-1., 1. + decay) * sqrt(2.) * sin(theta_K);

   double MD_R;
   MD_R = C * MD * R;

   delete Pion;
   delete Kaon;

   delete K1400_D;

   return MD_R;
} 
