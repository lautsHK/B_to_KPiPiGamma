//include headers
#include "Helicity_Amplitude.hpp"
#include "Kinematics.hpp"
#include "ThreeBodies_Kinematics.hpp"
#include "Form_Factors.hpp"
#include "Particle.hpp"
#include "Coupling_Constants.hpp"
#include "Functions.hpp"
#include "Matrix_Elements.hpp"

//including c++ libraries
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


std::complex<double> Matrix_Elements::Jx_Two_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400)
{
   //define the particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //define the kinematics
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;

   double p1_x = std::abs(p1) * cos_theta * cos(phi1);
   double p2_x = std::abs(p2) * cos_theta * cos(phi2);

   //***Breit-Wigner for K1270 and K1400***
   Functions* FN = new Functions;
 
   Particle* K1270 = new Particle(10313); //10313 is the PDGID of K1(1270)
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //20313 is the PDGID of K1(1400)
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)

   //***C1_res and C2_res***
   Helicity_Amplitude* HAs = new Helicity_Amplitude;
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = HAs->C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = HAs->C1_1400(s, s13, s23, k_charge);
 
   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400;
   
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = HAs->C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = HAs->C2_1400(s, s13, s23, k_charge);
 
   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400;


   //Calculate Ji
   std::complex<double> Ji_2R;
   Ji_2R = C1_res * p1_x - C2_res * p2_x;

   delete PiPiK;
   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete FN;
   delete HAs;

   return Ji_2R;
}


std::complex<double> Matrix_Elements::Jy_Two_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400)
{
   //define the particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //define the kinematics
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;

   double p1_y = std::abs(p1) * sin(phi1);
   double p2_y = std::abs(p2) * sin(phi2);

   //***Breit-Wigner for K1270 and K1400***
   Functions* FN = new Functions;
 
   Particle* K1270 = new Particle(10313); //10313 is the PDGID of K1(1270)
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //20313 is the PDGID of K1(1400)
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)

   //***C1_res and C2_res***
   Helicity_Amplitude* HAs = new Helicity_Amplitude;
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = HAs->C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = HAs->C1_1400(s, s13, s23, k_charge);
 
   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400;
   
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = HAs->C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = HAs->C2_1400(s, s13, s23, k_charge);
 
   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400;


   //Calculate Ji
   std::complex<double> Ji_2R;
   Ji_2R = C1_res * p1_y - C2_res * p2_y;

   delete PiPiK;
   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete FN;
   delete HAs;

   return Ji_2R;
}


std::complex<double> Matrix_Elements::Jz_Two_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400)
{
   //define the particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //define the kinematics

   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;

   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta
   double p1_z = -std::abs(p1) * sin_theta * cos(phi1);
   double p2_z = -std::abs(p2) * sin_theta * cos(phi2);
   

   //***Breit-Wigner for K1270 and K1400***
   Functions* FN = new Functions;
 
   Particle* K1270 = new Particle(10313); //10313 is the PDGID of K1(1270)
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //20313 is the PDGID of K1(1400)
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)

   //***C1_res and C2_res***
   Helicity_Amplitude* HAs = new Helicity_Amplitude;
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = HAs->C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = HAs->C1_1400(s, s13, s23, k_charge);
 
   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400;
   
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = HAs->C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = HAs->C2_1400(s, s13, s23, k_charge);
 
   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400;


   //Calculate Ji
   std::complex<double> Ji_2R;
   Ji_2R = C1_res * p1_z - C2_res * p2_z;

   delete PiPiK;
   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete FN;
   delete HAs;

   return Ji_2R;
}


std::complex<double> Matrix_Elements::Jx_All_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680)
{
   //define the particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //define the kinematics
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
 
   double p1_x = std::abs(p1) * cos_theta * cos(phi1);
   double p2_x = std::abs(p2) * cos_theta * cos(phi2);

   double p1_y = std::abs(p1) * sin(phi1);
   double p2_y = std::abs(p2) * sin(phi2);

   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta
   double p1_z = -std::abs(p1) * sin_theta * cos(phi1);
   double p2_z = -std::abs(p2) * sin_theta * cos(phi2);

   double p1Xp2_x = (p1_y)*(p2_z) - (p1_z)*(p2_y);

   //***Breit-Wigner for K_res***
   Functions* FN = new Functions;
 
   Particle* K1270 = new Particle(10313); //10313 is the PDGID of K1(1270)
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //10323 is the PDGID of K1(1400)
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   Particle* K1410 = new Particle(100313); //100313 is the PDGID of K*(1410)
   double mass_K1410 = K1410->Get_Particle_Mass();
   double width_K1410 = K1410->Get_Particle_Width();

   Particle* K1430 = new Particle(315); //315 is the PDGID of K*2(1430)
   double mass_K1430 = K1430->Get_Particle_Mass();
   double width_K1430 = K1430->Get_Particle_Width();

   Particle* K1680 = new Particle(30313); //30313 is the PDGID of K*(1680)
   double mass_K1680 = K1680->Get_Particle_Mass();
   double width_K1680 = K1680->Get_Particle_Width();

   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)
   std::complex<double> BW_K1410 = FN->BWPropagator(s, mass_K1410 , width_K1410); //(s, mass, width)
   std::complex<double> BW_K1430 = FN->BWPropagator(s, mass_K1430 , width_K1430); //(s, mass, width)
   std::complex<double> BW_K1680 = FN->BWPropagator(s, mass_K1680 , width_K1680); //(s, mass, width)

   //***C1_res, C2_res, and C3_res***
   Helicity_Amplitude* HA_new = new Helicity_Amplitude;

   //C1
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = HA_new->C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = HA_new->C1_1400(s, s13, s23, k_charge);

   std::complex<double> C1_K1430;
   C1_K1430 = HA_new->C1_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400 + T1K1430 * BW_K1430 * C1_K1430;
   
   //C2
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = HA_new->C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = HA_new->C2_1400(s, s13, s23, k_charge);
 
   std::complex<double> C2_K1430;
   C2_K1430 = HA_new->C2_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400 + T1K1430 * BW_K1430 * C2_K1430;

   //C3
   std::complex<double> C3_res;

   std::complex<double> C3_K1410;
   C3_K1410 = HA_new->C3_1410(s, s13, s23, k_charge);
 
   std::complex<double> C3_K1680;
   C3_K1680 = HA_new->C3_1680(s, s13, s23, k_charge);
 
   std::complex<double> C3_K1430;
   C3_K1430 = HA_new->C3_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C3_res = T1K1410 * BW_K1410 * C3_K1410 + T1K1680 * BW_K1680 * C3_K1680 + T1K1430 * BW_K1430 * C3_K1430;

   //std::cout << "The value of C3 in Jx is " << C3_res << " .\n";

   //Calculate Ji
   std::complex<double> Ji_All;
   Ji_All = C1_res * p1_x - C2_res * p2_x + C3_res * p1Xp2_x;

   delete PiPiK;
   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete K1410;
   delete K1680;
   delete K1430;
   delete FN;
   delete HA_new;

   return Ji_All;
}


std::complex<double> Matrix_Elements::Jy_All_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680)
{
   //define the particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //define the kinematics
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;

   double p1_x = std::abs(p1) * cos_theta * cos(phi1);
   double p2_x = std::abs(p2) * cos_theta * cos(phi2);

   double p1_y = std::abs(p1) * sin(phi1);
   double p2_y = std::abs(p2) * sin(phi2);

   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta
   double p1_z = -std::abs(p1) * sin_theta * cos(phi1);
   double p2_z = -std::abs(p2) * sin_theta * cos(phi2);

   double p1Xp2_y = (p1_z)*(p2_x) - (p1_x)*(p2_z);

   //***Breit-Wigner for K_res***
   Functions* FN = new Functions;
 
   Particle* K1270 = new Particle(10313); //10313 is the PDGID of K1(1270)
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //20313 is the PDGID of K1(1400)
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   Particle* K1410 = new Particle(100313); //100313 is the PDGID of K*(1410)
   double mass_K1410 = K1410->Get_Particle_Mass();
   double width_K1410 = K1410->Get_Particle_Width();

   Particle* K1430 = new Particle(315); //315 is the PDGID of K*2(1430)
   double mass_K1430 = K1430->Get_Particle_Mass();
   double width_K1430 = K1430->Get_Particle_Width();

   Particle* K1680 = new Particle(30313); //30313 is the PDGID of K*(1680)
   double mass_K1680 = K1680->Get_Particle_Mass();
   double width_K1680 = K1680->Get_Particle_Width();

   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)
   std::complex<double> BW_K1410 = FN->BWPropagator(s, mass_K1410 , width_K1410); //(s, mass, width)
   std::complex<double> BW_K1430 = FN->BWPropagator(s, mass_K1430 , width_K1430); //(s, mass, width)
   std::complex<double> BW_K1680 = FN->BWPropagator(s, mass_K1680 , width_K1680); //(s, mass, width)

   //***C1_res, C2_res, and C3_res***
   Helicity_Amplitude* HA_new = new Helicity_Amplitude;

   //C1
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = HA_new->C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = HA_new->C1_1400(s, s13, s23, k_charge);

   std::complex<double> C1_K1430;
   C1_K1430 = HA_new->C1_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400 + T1K1430 * BW_K1430 * C1_K1430;
   
   //C2
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = HA_new->C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = HA_new->C2_1400(s, s13, s23, k_charge);
 
   std::complex<double> C2_K1430;
   C2_K1430 = HA_new->C2_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400 + T1K1430 * BW_K1430 * C2_K1430;

   //C3
   std::complex<double> C3_res;

   std::complex<double> C3_K1410;
   C3_K1410 = HA_new->C3_1410(s, s13, s23, k_charge);
 
   std::complex<double> C3_K1680;
   C3_K1680 = HA_new->C3_1680(s, s13, s23, k_charge);
 
   std::complex<double> C3_K1430;
   C3_K1430 = HA_new->C3_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C3_res = T1K1410 * BW_K1410 * C3_K1410 + T1K1680 * BW_K1680 * C3_K1680 + T1K1430 * BW_K1430 * C3_K1430;

   //std::cout << "The value of C3 in Jy is " << C3_res << " .\n";

   //Calculate Ji
   std::complex<double> Ji_All;
   Ji_All = C1_res * p1_y - C2_res * p2_y + C3_res * p1Xp2_y;

   delete PiPiK;
   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete K1410;
   delete K1680;
   delete K1430;
   delete FN;
   delete HA_new;

   return Ji_All;
}


std::complex<double> Matrix_Elements::Jz_All_Resonances( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680)
{
   //define the particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   //define the kinematics

   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;

   double p1_x = std::abs(p1) * cos_theta * cos(phi1);
   double p2_x = std::abs(p2) * cos_theta * cos(phi2);

   double p1_y = std::abs(p1) * sin(phi1);
   double p2_y = std::abs(p2) * sin(phi2);

   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta
   double p1_z = -std::abs(p1) * sin_theta * cos(phi1);
   double p2_z = -std::abs(p2) * sin_theta * cos(phi2);
   
   double p1Xp2_z = (p1_x)*(p2_y) - (p1_y)*(p2_x);

   //***Breit-Wigner for K_res***
   Functions* FN = new Functions;
 
   Particle* K1270 = new Particle(10313); //10313 is the PDGID of K1(1270)
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //20313 is the PDGID of K1(1400)
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   Particle* K1410 = new Particle(100313); //100313 is the PDGID of K*(1410)
   double mass_K1410 = K1410->Get_Particle_Mass();
   double width_K1410 = K1410->Get_Particle_Width();

   Particle* K1430 = new Particle(315); //315 is the PDGID of K*2(1430)
   double mass_K1430 = K1430->Get_Particle_Mass();
   double width_K1430 = K1430->Get_Particle_Width();

   Particle* K1680 = new Particle(30313); //30313 is the PDGID of K*(1680)
   double mass_K1680 = K1680->Get_Particle_Mass();
   double width_K1680 = K1680->Get_Particle_Width();

   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)
   std::complex<double> BW_K1410 = FN->BWPropagator(s, mass_K1410 , width_K1410); //(s, mass, width)
   std::complex<double> BW_K1430 = FN->BWPropagator(s, mass_K1430 , width_K1430); //(s, mass, width)
   std::complex<double> BW_K1680 = FN->BWPropagator(s, mass_K1680 , width_K1680); //(s, mass, width)

   //***C1_res, C2_res, and C3_res***
   Helicity_Amplitude* HA_new = new Helicity_Amplitude;

   //C1
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = HA_new->C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = HA_new->C1_1400(s, s13, s23, k_charge);

   std::complex<double> C1_K1430;
   C1_K1430 = HA_new->C1_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400 + T1K1430 * BW_K1430 * C1_K1430;
   
   //C2
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = HA_new->C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = HA_new->C2_1400(s, s13, s23, k_charge);
 
   std::complex<double> C2_K1430;
   C2_K1430 = HA_new->C2_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400 + T1K1430 * BW_K1430 * C2_K1430;

   //C3
   std::complex<double> C3_res;

   std::complex<double> C3_K1410;
   C3_K1410 = HA_new->C3_1410(s, s13, s23, k_charge);
 
   std::complex<double> C3_K1680;
   C3_K1680 = HA_new->C3_1680(s, s13, s23, k_charge);
 
   std::complex<double> C3_K1430;
   C3_K1430 = HA_new->C3_1430(s, s13, s23, cos_theta, phi1, phi2, k_charge);

   C3_res = T1K1410 * BW_K1410 * C3_K1410 + T1K1680 * BW_K1680 * C3_K1680 + T1K1430 * BW_K1430 * C3_K1430;

   //std::cout << "The value of C3 in Jz is " << C3_res << " .\n";

   //Calculate Ji
   std::complex<double> Ji_All;
   Ji_All = C1_res * p1_z - C2_res * p2_z + C3_res * p1Xp2_z;

   delete PiPiK;
   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete K1410;
   delete K1680;
   delete K1430;
   delete FN;
   delete HA_new;

   return Ji_All;
}

double Matrix_Elements::J2( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680)
{
   std::complex<double> Jx = Jx_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jy = Jy_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jz = Jz_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);    

   std::complex<double> Jc_x = std::conj(Jx); 
   std::complex<double> Jc_y = std::conj(Jy); 
   std::complex<double> Jc_z = std::conj(Jz); 

   std::complex<double> J_square = Jx*Jc_x + Jy*Jc_y + Jz*Jc_z;
   double value_J2 = std::real(J_square);

   return value_J2;
}

double Matrix_Elements::Im_N_dot_J_cross_Jc( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680)
{
   //define the particles
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   delete Pion;
   delete Kaon;

   //The Vector N
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;

   double p1_x = std::abs(p1) * cos_theta * cos(phi1);
   double p2_x = std::abs(p2) * cos_theta * cos(phi2);

   double p1_y = std::abs(p1) * sin(phi1);
   double p2_y = std::abs(p2) * sin(phi2);

   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta
   double p1_z = -std::abs(p1) * sin_theta * cos(phi1);
   double p2_z = -std::abs(p2) * sin_theta * cos(phi2);
   
   double p1Xp2_x = (p1_y)*(p2_z) - (p1_z)*(p2_y);
   double p1Xp2_y = (p1_z)*(p2_x) - (p1_x)*(p2_z);
   double p1Xp2_z = (p1_x)*(p2_y) - (p1_y)*(p2_x);

   double p1Xp2 = PiPiK->p1_Cross_p2(s, s13, s23);

   delete PiPiK;

   //The Vector J
   std::complex<double> Jx = Jx_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jy = Jy_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jz = Jz_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);    

   std::complex<double> Jc_x = std::conj(Jx); 
   std::complex<double> Jc_y = std::conj(Jy); 
   std::complex<double> Jc_z = std::conj(Jz); 

   std::complex<double> JXJc_x = (Jy)*(Jc_z) - (Jz)*(Jc_y);
   std::complex<double> JXJc_y = (Jz)*(Jc_x) - (Jx)*(Jc_z);
   std::complex<double> JXJc_z = (Jx)*(Jc_y) - (Jy)*(Jc_x);

   std::complex<double> nJXJc = (p1Xp2_x * JXJc_x + p1Xp2_y * JXJc_y + p1Xp2_z * JXJc_z) / p1Xp2; 
   double Im_nJXJc = std::imag(nJXJc);            

   return Im_nJXJc;
}

double Matrix_Elements::ME2_LH( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680)
{
   std::complex<double> Jx = Jx_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jy = Jy_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jz = Jz_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);    

   double phase_space = 1. / (64. * pow((2.*M_PI), 4.) * pow( s, 1.5 ));
   std::complex<double> Im_i(0.0 , 1.0);
   std::complex<double> ME_L = (Jx - Im_i * Jy )/sqrt(2.);
   std::complex<double> ME_L2 = ME_L * std::conj(ME_L); 
   double value_ME2 = std::real(ME_L2) * phase_space;

   return value_ME2;
}

double Matrix_Elements::ME2_RH( double s, double s13, double s23, double cos_theta, double phi, int k_charge, double T1K1270, double T1K1400, double T1K1410, double T1K1430, double T1K1680)
{
   std::complex<double> Jx = Jx_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jy = Jy_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);      
   std::complex<double> Jz = Jz_All_Resonances(s, s13, s23, cos_theta, phi, k_charge, T1K1270, T1K1400, T1K1410, T1K1430, T1K1680);    

   double phase_space = 1. / (64. * pow((2.*M_PI), 4.) * pow( s, 1.5 ));
   std::complex<double> Im_i(0.0 , 1.0);
   std::complex<double> ME_R = (Jx + Im_i * Jy )/sqrt(2.); //This line should be the only line leading to a different result between right and left. 
   std::complex<double> ME_R2 = ME_R * std::conj(ME_R); 
   double value_ME2 = std::real(ME_R2) * phase_space;

   return value_ME2;
}
