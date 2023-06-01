//including headers
#include "Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp"
#include "Kinematics.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>

double QPCM_A_to_VplusP::Get_gammaPQC() const
{
   double g;
   g = gammaPQC;

   return g;
}

double QPCM_A_to_VplusP::Get_f2() const
{
   double f;
   f = f2;

   return f;
}

double QPCM_A_to_VplusP::Get_R_SHO() const
{
   double R;
   R = R_SHO;

   return R;
}

double QPCM_A_to_VplusP::I0(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P)
{
   double sum_all_R2 = pow(radius_A, 2.0) + pow(radius_V, 2.0) + pow(radius_P, 2.0);
   double sum_final_R2 = pow(radius_V, 2.0) + pow(radius_P, 2.0);

   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double C = 4.*sqrt(3.0)/pow(M_PI, 5./4.);
   double R1 = pow(radius_A, 5./2.) * pow(radius_V*radius_P, 3./2.) / pow(sum_all_R2, 5./2.);
   double R2 = ( 2.*pow(radius_A, 2.0) + sum_final_R2 ) * (sum_final_R2) * p2 / (4.*sum_all_R2);
   double P = ( pow(radius_A, 2.0) * sum_final_R2 ) * p2 / (8.*sum_all_R2);

   double value_I0;
   value_I0 = -C * R1 * (1. - R2) * exp(-P);

   delete kine;

   return value_I0;
}

//Overloading I0, with R_SHO approximation
double QPCM_A_to_VplusP::I0(double mass_A, double mass_V, double mass_P)
{
   double radius_A = R_SHO; //the value is defined in header file.
   double radius_V = R_SHO;
   double radius_P = R_SHO;

   double sum_all_R2 = pow(radius_A, 2.0) + pow(radius_V, 2.0) + pow(radius_P, 2.0);
   double sum_final_R2 = pow(radius_V, 2.0) + pow(radius_P, 2.0);

   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double C = 4.*sqrt(3.0)/pow(M_PI, 5./4.);
   double R1 = pow(radius_A, 5./2.) * pow(radius_V*radius_P, 3./2.) / pow(sum_all_R2, 5./2.);
   double R2 = ( 2.*pow(radius_A, 2.0) + sum_final_R2 ) * (sum_final_R2) * p2 / (4.*sum_all_R2);
   double P = ( pow(radius_A, 2.0) * sum_final_R2 ) * p2 / (8.*sum_all_R2);

   double value_I0;
   value_I0 = -C * R1 * (1. - R2) * exp(-P);

   delete kine;

   return value_I0;
}

double QPCM_A_to_VplusP::I1(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P)
{
   double sum_all_R2 = pow(radius_A, 2.0) + pow(radius_V, 2.0) + pow(radius_P, 2.0);
   double sum_final_R2 = pow(radius_V, 2.0) + pow(radius_P, 2.0);

   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double C = 4.*sqrt(3.0)/pow(M_PI, 5./4.);
   double R1 = pow(radius_A, 5./2.) * pow(radius_V*radius_P, 3./2.) / pow(sum_all_R2, 5./2.);
   double P = ( pow(radius_A, 2.0) * sum_final_R2 ) * p2 / (8.*sum_all_R2);

   double value_I1;
   value_I1 = C * R1 * exp(-P);

   delete kine;

   return value_I1;
}

//Overloading I0, with R_SHO approximation
double QPCM_A_to_VplusP::I1(double mass_A, double mass_V, double mass_P)
{
   double radius_A = R_SHO; //the value is defined in header file.
   double radius_V = R_SHO;
   double radius_P = R_SHO;

   double sum_all_R2 = pow(radius_A, 2.0) + pow(radius_V, 2.0) + pow(radius_P, 2.0);
   double sum_final_R2 = pow(radius_V, 2.0) + pow(radius_P, 2.0);

   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double C = 4.*sqrt(3.0)/pow(M_PI, 5./4.);
   double R1 = pow(radius_A, 5./2.) * pow(radius_V*radius_P, 3./2.) / pow(sum_all_R2, 5./2.);
   double P = ( pow(radius_A, 2.0) * sum_final_R2 ) * p2 / (8.*sum_all_R2);

   double value_I1;
   value_I1 = C * R1 * exp(-P);

   delete kine;

   return value_I1;
}

double QPCM_A_to_VplusP::MSpqc(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P, double damp_fac)
{
   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double _I0 = I0(mass_A, mass_V, mass_P, radius_A, radius_V, radius_P);
   double _I1 = I1(mass_A, mass_V, mass_P, radius_A, radius_V, radius_P);

   double value_MSpqc;
   value_MSpqc = sqrt(3./2.) * gammaPQC * ( (2.*_I1 - _I0) / 18.) * exp(-damp_fac*p2);

   delete kine;

   return value_MSpqc;
}

//Overloading MSpqc
double QPCM_A_to_VplusP::MSpqc(double mass_A, double mass_V, double mass_P)
{
   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double _I0 = I0(mass_A, mass_V, mass_P);
   double _I1 = I1(mass_A, mass_V, mass_P);

   double damp_fac = f2; //f2 = 3.0 , stored in the header file

   double value_MSpqc;
   value_MSpqc = sqrt(3./2.) * gammaPQC * ( (2.*_I1 - _I0) / 18.) * exp(-damp_fac*p2);

   delete kine;

   return value_MSpqc;
}


double QPCM_A_to_VplusP::MDpqc(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P, double damp_fac)
{
   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double _I0 = I0(mass_A, mass_V, mass_P, radius_A, radius_V, radius_P);
   double _I1 = I1(mass_A, mass_V, mass_P, radius_A, radius_V, radius_P);

   double value_MDpqc;
   value_MDpqc = sqrt(3./2.) * gammaPQC * ( (_I1 + _I0) / 18.) * exp(-damp_fac*p2);

   delete kine;

   return value_MDpqc;
}

//Overloading MQpqc
double QPCM_A_to_VplusP::MDpqc(double mass_A, double mass_V, double mass_P)
{
   Kinematics* kine = new Kinematics;
   double p_final = kine->M2ab_CoP_Momentum(mass_A, mass_V, mass_P);
   double p2 = pow(p_final, 2.0);

   double _I0 = I0(mass_A, mass_V, mass_P);
   double _I1 = I1(mass_A, mass_V, mass_P);

   double damp_fac = f2; //f2 = 3.0 , stored in the header file

   double value_MDpqc;
   value_MDpqc = sqrt(3./2.) * gammaPQC * ( (_I1 + _I0) / 18.) * exp(-damp_fac*p2);

   delete kine;

   return value_MDpqc;
}



