//including headers
#include "Dalitz_Region.hpp"
#include "Kinematics.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>

Dalitz_Region::Dalitz_Region()
{
   m_Dalitz[0] = 0.0;
   m_Dalitz[1] = 0.0;
   m_Dalitz[2] = 0.0;
   m_Dalitz[3] = 0.0;
}

void Dalitz_Region::SetSMaxM1M2M3(double sMax, double mass1, double mass2, double mass3)
{
   assert(mass1 >= 0.0);
   assert(mass2 >= 0.0);
   assert(mass3 >= 0.0);
   assert(sMax >= 0.0);
   
   m_Dalitz[0] = sMax;
   m_Dalitz[1] = mass1;
   m_Dalitz[2] = mass2;
   m_Dalitz[3] = mass3;
}


double Dalitz_Region::GetDRMass1() const
{
   double mass1 = m_Dalitz[1];
   return mass1;
}

double Dalitz_Region::GetDRMass2() const
{
   double mass2 = m_Dalitz[2];
   return mass2;
}

double Dalitz_Region::GetDRMass3() const
{
   double mass3 = m_Dalitz[3];
   return mass3;
}

double Dalitz_Region::GetSMax() const
{
   double Smax = m_Dalitz[0];
   return Smax;
}

double Dalitz_Region::GetSMin()
{
   double Smin = pow( (m_Dalitz[1] + m_Dalitz[2] + m_Dalitz[3]) , 2.0);
   return Smin;
}

double Dalitz_Region::Get_Sij_Min(int i, int j)
{
  double Sij_min = pow( (m_Dalitz[i] + m_Dalitz[j]) , 2.0);
  return Sij_min;
}

double Dalitz_Region::Get_Sij_Max(int i, int j)
{
  double Sij_max;
  //Sij_max = ((m_Dalitz[0])^0.5 - m_Dalitz[k])^2;
  //But we can't directly get k
  double sMax = m_Dalitz[0];
  double Dalitz_k = m_Dalitz[1] + m_Dalitz[2] + m_Dalitz[3] - m_Dalitz[i] - m_Dalitz[j];
  Sij_max = pow( (sqrt( sMax ) - Dalitz_k) , 2.0);
  return Sij_max;
}

double Dalitz_Region::Get_S23_Dalitz_Min_wrt_S13(double s, double s13)
{
   assert(s13 > 0.0);
   assert(s > 0.0);

   double m1 = m_Dalitz[1];
   double m2 = m_Dalitz[2];
   double m3 = m_Dalitz[3];

   Kinematics* kine = new Kinematics;

   double E2_13rest;
   E2_13rest = kine->M2ab_Energy_a_wrt_b(sqrt(s), m2, sqrt(s13));
   //E2_13rest = (s - s13 - pow(m2,2.0) ) / (2.0*sqrt( s13 ) );

   double E3_13rest;
   E3_13rest = kine->M2ab_CoP_Energy_a(sqrt(s13), m3, m1);
   //E3_13rest = (s13 - pow(m1,2.0) + pow(m3,2.0) ) / (2.0*sqrt( s13 ) );

   double p2_13rest;
   double p3_13rest;
   p2_13rest = sqrt( pow(E2_13rest,2.0) - pow(m2,2.0) );
   p3_13rest = sqrt( pow(E3_13rest,2.0) - pow(m3,2.0) );

   double s23_Dmin; 

   if( (E2_13rest >= m2) && (E3_13rest >= m3) )
   {
      s23_Dmin = pow( (E2_13rest + E3_13rest) , 2.0) - pow( (p2_13rest + p3_13rest) , 2.0) ;
   }
   else
   {
      s23_Dmin = 0.0;
   }
  
   delete kine;

   return s23_Dmin; 
}

double Dalitz_Region::Get_S23_Dalitz_Max_wrt_S13(double s, double s13)
{
   assert(s13 > 0.0);
   assert(s > 0.0);

   double m1 = m_Dalitz[1];
   double m2 = m_Dalitz[2];
   double m3 = m_Dalitz[3];

   Kinematics* kine = new Kinematics;

   double E2_13rest;
   E2_13rest = kine->M2ab_Energy_a_wrt_b(sqrt(s), m2, sqrt(s13));
   //E2_13rest = (s - s13 - pow(m2,2.0) ) / (2.0*sqrt( s13 ) );

   double E3_13rest;
   E3_13rest = kine->M2ab_CoP_Energy_a(sqrt(s13), m3, m1);
   //E3_13rest = (s13 - pow(m1,2.0) + pow(m3,2.0) ) / (2.0*sqrt( s13 ) );


   double p2_13rest;
   double p3_13rest;
   p2_13rest = sqrt( pow(E2_13rest,2.0) - pow(m2,2.0) );
   p3_13rest = sqrt( pow(E3_13rest,2.0) - pow(m3,2.0) );

   double s23_Dmax; 

   if( (E2_13rest >= m2) && (E3_13rest >= m3) )
   {
      s23_Dmax = pow( (E2_13rest + E3_13rest) , 2.0) - pow( (p2_13rest - p3_13rest) , 2.0) ;
   }
   else
   {
      s23_Dmax = 0.0;
   }

   delete kine;

   return s23_Dmax; 
}
