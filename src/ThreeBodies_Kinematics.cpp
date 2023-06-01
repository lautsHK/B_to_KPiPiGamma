//including headers
#include "ThreeBodies_Kinematics.hpp"
#include "Kinematics.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>

ThreeBodies_Kinematics::ThreeBodies_Kinematics()
{
   m_mass[0] = 0.0;
   m_mass[1] = 0.0;
   m_mass[2] = 0.0;
}

ThreeBodies_Kinematics::ThreeBodies_Kinematics(double mass1, double mass2, double mass3)
{
   assert(mass1 >= 0.0);
   assert(mass2 >= 0.0);
   assert(mass3 >= 0.0);

   m_mass[0] = mass1;
   m_mass[1] = mass2;
   m_mass[2] = mass3;
}

void ThreeBodies_Kinematics::SetM1M2M3(double mass1, double mass2, double mass3)
{
   assert(mass1 >= 0.0);
   assert(mass2 >= 0.0);
   assert(mass3 >= 0.0);

   m_mass[0] = mass1;
   m_mass[1] = mass2;
   m_mass[2] = mass3;
}

double ThreeBodies_Kinematics::Get_s12(double s, double s13, double s23)
{
   double m1 = m_mass[0];
   double m2 = m_mass[1];
   double m3 = m_mass[2];
   
   double s12;
   s12 = s - s13 - s23 + pow(m1, 2.0) + pow(m2, 2.0) + pow(m3, 2.0);

   return s12;
}

double ThreeBodies_Kinematics::GetMass1() const
{
   double mass1 = m_mass[0];
   return mass1;
}

double ThreeBodies_Kinematics::GetMass2() const
{
   double mass2 = m_mass[1];
   return mass2;
}

double ThreeBodies_Kinematics::GetMass3() const
{
   double mass3 = m_mass[2];
   return mass3;
}

double ThreeBodies_Kinematics::Get_CoP_Energy1(double s, double s13, double s23)
{
   double m1 = m_mass[0];

   Kinematics* kine = new Kinematics;
   
   double E1;
   E1 = kine->M2ab_CoP_Energy_a(sqrt(s), m1, sqrt(s23));

   delete kine;

   return E1;
}

double ThreeBodies_Kinematics::Get_CoP_Energy2(double s, double s13, double s23)
{
   double m2 = m_mass[1];

   Kinematics* kine = new Kinematics;
   
   double E2;
   E2 = kine->M2ab_CoP_Energy_a(sqrt(s), m2, sqrt(s13));

   delete kine;

   return E2;
}

double ThreeBodies_Kinematics::Get_CoP_Energy3(double s, double s13, double s23)
{
   double m3 = m_mass[2];
   
   double s12 = Get_s12(s, s13, s23);

   Kinematics* kine = new Kinematics;
   
   double E3;

   E3 = kine->M2ab_CoP_Energy_a(sqrt(s), m3, sqrt(s12));

   delete kine;

   return E3;
}

double ThreeBodies_Kinematics::Get_CoP_Momentum1(double s, double s13, double s23)
{
   double m1 = m_mass[0];

   Kinematics* kine = new Kinematics;
   
   double p1;
   p1 = kine->M2ab_CoP_Momentum(sqrt(s), m1, sqrt(s23));

   delete kine;

   return p1;
}

double ThreeBodies_Kinematics::Get_CoP_Momentum2(double s, double s13, double s23)
{
   double m2 = m_mass[1];

   Kinematics* kine = new Kinematics;
   
   double p2;
   p2 = kine->M2ab_CoP_Momentum(sqrt(s), m2, sqrt(s13));

   delete kine;

   return p2;
}

double ThreeBodies_Kinematics::Get_CoP_Momentum3(double s, double s13, double s23)
{
   double m3 = m_mass[2];
   
   double s12 = Get_s12(s, s13, s23);

   Kinematics* kine = new Kinematics;
   
   double p3;

   p3 = kine->M2ab_CoP_Momentum(sqrt(s), m3, sqrt(s12));

   delete kine;

   return p3;
}

double ThreeBodies_Kinematics::p1_3Dot_p2(double s, double s13, double s23)
{
   double m1 = m_mass[0];
   double m1_sq = pow(m1, 2.0);
   double m2 = m_mass[1];
   double m2_sq = pow(m2, 2.0);
 
   double E1;
   E1 = Get_CoP_Energy1(s, s13, s23);
   
   double E2;
   E2 = Get_CoP_Energy2(s, s13, s23);

   double s12 = Get_s12(s, s13, s23);

   //calculate the dot product 
   double p1p2;
   p1p2 = E1*E2 + (m1_sq + m2_sq - s12)/2.0;

   return p1p2;
}

double ThreeBodies_Kinematics::p1_Cross_p2(double s, double s13, double s23)
{
   double p1 = Get_CoP_Momentum1(s, s13, s23);
   double p2 = Get_CoP_Momentum2(s, s13, s23);
   double p1p2 = p1_3Dot_p2(s, s13, s23);

   //calculate the cross product 
   double p1xp2;
   p1xp2 = p1 * p2 * sin( acos( p1p2/(p1*p2) ) );

   return p1xp2;
}

double ThreeBodies_Kinematics::p1_4Dot_p2(double s, double s13, double s23)
{
   double m1 = m_mass[0];
   double m1_sq = pow(m1, 2.0);
   double m2 = m_mass[1];
   double m2_sq = pow(m2, 2.0);

   double s12 = Get_s12(s, s13, s23);

   //calculate the dot product 
   double p1p2;
   p1p2 = (s12 - m1_sq - m2_sq)/2.0;

   return p1p2;
}

double ThreeBodies_Kinematics::Get_CoP_P1x(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1 = Get_CoP_Momentum1(s, s13, s23);
   double p2 = Get_CoP_Momentum2(s, s13, s23);

   double p1p2 = p1_3Dot_p2(s, s13, s23);
   
   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   
   double p1_x = std::abs(p1) * cos_theta * cos(phi1);

   return p1_x;
}

double ThreeBodies_Kinematics::Get_CoP_P1y(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1 = Get_CoP_Momentum1(s, s13, s23);
   double p2 = Get_CoP_Momentum2(s, s13, s23);

   double p1p2 = p1_3Dot_p2(s, s13, s23);
   
   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   
   double p1_y = std::abs(p1) * sin(phi1);

   return p1_y;
}

double ThreeBodies_Kinematics::Get_CoP_P1z(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1 = Get_CoP_Momentum1(s, s13, s23);
   double p2 = Get_CoP_Momentum2(s, s13, s23);

   double p1p2 = p1_3Dot_p2(s, s13, s23);
   
   double phi1 = phi - acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   
   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta
   
   double p1_z = -std::abs(p1) * sin_theta * cos(phi1);

   return p1_z;
}

double ThreeBodies_Kinematics::Get_CoP_P2x(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1 = Get_CoP_Momentum1(s, s13, s23);
   double p2 = Get_CoP_Momentum2(s, s13, s23);

   double p1p2 = p1_3Dot_p2(s, s13, s23);
   
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   
   double p2_x = std::abs(p2) * cos_theta * cos(phi2);

   return p2_x;
}

double ThreeBodies_Kinematics::Get_CoP_P2y(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1 = Get_CoP_Momentum1(s, s13, s23);
   double p2 = Get_CoP_Momentum2(s, s13, s23);

   double p1p2 = p1_3Dot_p2(s, s13, s23);
   
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   
   double p2_y = std::abs(p2) * sin(phi2);

   return p2_y;
}

double ThreeBodies_Kinematics::Get_CoP_P2z(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1 = Get_CoP_Momentum1(s, s13, s23);
   double p2 = Get_CoP_Momentum2(s, s13, s23);

   double p1p2 = p1_3Dot_p2(s, s13, s23);
   
   double phi2 = phi + acos( p1p2 / (std::abs(p1)*std::abs(p2)) )/2.;
   
   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta
   
   double p2_z = -std::abs(p2) * sin_theta * cos(phi2);

   return p2_z;
}

double ThreeBodies_Kinematics::Get_CoP_P3x(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1_x = Get_CoP_P1x(s, s13, s23, cos_theta, phi);
   double p2_x = Get_CoP_P2x(s, s13, s23, cos_theta, phi);

   double p3_x = -p1_x - p2_x ;

   return p3_x;
}

double ThreeBodies_Kinematics::Get_CoP_P3y(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1_y = Get_CoP_P1y(s, s13, s23, cos_theta, phi);
   double p2_y = Get_CoP_P2y(s, s13, s23, cos_theta, phi);

   double p3_y = -p1_y - p2_y ;

   return p3_y;
}

double ThreeBodies_Kinematics::Get_CoP_P3z(double s, double s13, double s23, double cos_theta, double phi)
{
   double p1_z = Get_CoP_P1z(s, s13, s23, cos_theta, phi);
   double p2_z = Get_CoP_P2z(s, s13, s23, cos_theta, phi);

   double p3_z = -p1_z - p2_z ;

   return p3_z;
}



