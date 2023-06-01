//This is a header file for calculating the 3-bodies' kinematics
#ifndef THREEBODIES_KINEMATICS_HPP
#define THREEBODIES_KINEMATICS_HPP

#include "Kinematics.hpp"

class ThreeBodies_Kinematics
{
public :

   // Some background information
   // s_12 + s_23 + s_13 = s + m1^2 + m2^2 +m3^2
   // s_ij is the invariant mass of the pair of particles i and j
   // That is s_ij = (4p_i+4p_j)^2, where 4p is the 4-momentum

   ThreeBodies_Kinematics();
   ThreeBodies_Kinematics(double mass1, double mass2, double mass3);
   void SetM1M2M3(double mass1, double mass2, double mass3);

   double Get_s12(double s, double s13, double s23);

   double GetMass1() const;
   double GetMass2() const;
   double GetMass3() const;

   double Get_CoP_Energy1(double s, double s13, double s23);
   double Get_CoP_Energy2(double s, double s13, double s23);
   double Get_CoP_Energy3(double s, double s13, double s23);

   double Get_CoP_Momentum1(double s, double s13, double s23);
   double Get_CoP_Momentum2(double s, double s13, double s23);
   double Get_CoP_Momentum3(double s, double s13, double s23);

   double Get_CoP_P1x(double s, double s13, double s23, double cos_theta, double phi);
   double Get_CoP_P2x(double s, double s13, double s23, double cos_theta, double phi);
   double Get_CoP_P3x(double s, double s13, double s23, double cos_theta, double phi);

   double Get_CoP_P1y(double s, double s13, double s23, double cos_theta, double phi);
   double Get_CoP_P2y(double s, double s13, double s23, double cos_theta, double phi);
   double Get_CoP_P3y(double s, double s13, double s23, double cos_theta, double phi);

   double Get_CoP_P1z(double s, double s13, double s23, double cos_theta, double phi);
   double Get_CoP_P2z(double s, double s13, double s23, double cos_theta, double phi);
   double Get_CoP_P3z(double s, double s13, double s23, double cos_theta, double phi);

   double p1_3Dot_p2(double s, double s13, double s23);
   double p1_Cross_p2(double s, double s13, double s23);
   double p1_4Dot_p2(double s, double s13, double s23);

private :
 
   double m_mass[3];
};

#endif
