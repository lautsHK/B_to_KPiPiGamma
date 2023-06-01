//including headers
#include "Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp"
#include "Quark_Pair_Creation_Model/Kres_to_VP.hpp"
#include "Kinematics.hpp"
#include "Particle.hpp"
#include "Form_Factors.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>

double Form_Factors::A_K1270(double s, double s_ij, int decay)
{ 
   Kres_to_VP* K1270_to_VP = new Kres_to_VP;

   double ME_S;
   double ME_D;

   ME_S = K1270_to_VP->MatrixElement_S_K1270(sqrt(s), sqrt(s_ij), decay);
   ME_D = K1270_to_VP->MatrixElement_D_K1270(sqrt(s), sqrt(s_ij), decay);

   double A;

   A = -( ME_S + ME_D/sqrt(2.) );

   delete K1270_to_VP;

   return A;
}

double Form_Factors::B_K1270(double s, double s_ij, int decay)
{
   //define particles and kinematics
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double mass_P[2]; //2 decay modes
   mass_P[0] = mass_pi; //K1 to K0st + Pi
   mass_P[1] = mass_K;  //K1 to rho + K

   Kinematics* kine = new Kinematics;
   double E_P = kine->M2ab_CoP_Energy_a( sqrt(s), sqrt(s_ij) , mass_P[decay]);
   double p3_P = kine->M2ab_CoP_Momentum( sqrt(s), sqrt(s_ij) , mass_P[decay]);

   //Coefficient calculations
   double C;
   C = E_P / ( sqrt(s) * pow(p3_P, 2.) );
   
   double C_S;
   C_S = 1. - sqrt(s_ij)/E_P;

   double C_D;
   C_D = 1. + 2.*sqrt(s_ij)/E_P;
 
   //Matrix Elements
   Kres_to_VP* K1270_to_VP = new Kres_to_VP;

   double ME_S;
   double ME_D;

   ME_S = K1270_to_VP->MatrixElement_S_K1270(sqrt(s), sqrt(s_ij), decay);
   ME_D = K1270_to_VP->MatrixElement_D_K1270(sqrt(s), sqrt(s_ij), decay);

   //Form Factor Calculation
   double B;
   B = -C * (C_S * ME_S + C_D * ME_D / sqrt(2.));

   delete Pion;
   delete Kaon;

   delete kine;

   delete K1270_to_VP;

   return B;
}

double Form_Factors::A_K1400(double s, double s_ij, int decay)
{ 
   Kres_to_VP* K1400_to_VP = new Kres_to_VP;

   double ME_S;
   double ME_D;

   ME_S = K1400_to_VP->MatrixElement_S_K1400(sqrt(s), sqrt(s_ij), decay);
   ME_D = K1400_to_VP->MatrixElement_D_K1400(sqrt(s), sqrt(s_ij), decay);

   double A;

   A = -( ME_S + ME_D/sqrt(2.) );

   delete K1400_to_VP;

   return A;
}

double Form_Factors::B_K1400(double s, double s_ij, int decay)
{
   //define particles and kinematics
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double mass_P[2]; //2 decay modes
   mass_P[0] = mass_pi; //K1 to K0st + Pi
   mass_P[1] = mass_K;  //K1 to rho + K

   Kinematics* kine = new Kinematics;
   double E_P = kine->M2ab_CoP_Energy_a( sqrt(s), sqrt(s_ij) , mass_P[decay]);
   double p3_P = kine->M2ab_CoP_Momentum( sqrt(s), sqrt(s_ij) , mass_P[decay]);

   //Coefficient calculations
   double C;
   C = E_P / ( sqrt(s) * pow(p3_P, 2.) );
   
   double C_S;
   C_S = 1. - sqrt(s_ij)/E_P;

   double C_D;
   C_D = 1. + 2.*sqrt(s_ij)/E_P;
 
   //Matrix Elements
   Kres_to_VP* K1400_to_VP = new Kres_to_VP;

   double ME_S;
   double ME_D;

   ME_S = K1400_to_VP->MatrixElement_S_K1400(sqrt(s), sqrt(s_ij), decay);
   ME_D = K1400_to_VP->MatrixElement_D_K1400(sqrt(s), sqrt(s_ij), decay);

   //Form Factor Calculation
   double B;
   B = -C * (C_S * ME_S + C_D * ME_D / sqrt(2.));

   delete Pion;
   delete Kaon;

   delete kine;

   delete K1400_to_VP;

   return B;
}
