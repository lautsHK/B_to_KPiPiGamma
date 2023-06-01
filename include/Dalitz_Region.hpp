//This is a header file for setting the Dalitz region
#ifndef DALITZ_REGION_HPP
#define DALITZ_REGION_HPP

#include "Kinematics.hpp"

class Dalitz_Region
{
public :

   // Some background information
   // s_12 + s_23 + s_13 = s + m1^2 + m2^2 +m3^2
   // s_ij is the invariant mass of the pair of particles i and j
   // That is s_ij = (4p_i+4p_j)^2, where 4p is the 4-momentum

   Dalitz_Region();
   void SetSMaxM1M2M3(double sMax, double mass1, double mass2, double mass3); //constructor of m1, m2, m3, SqrtS
   
   double GetDRMass1() const;
   double GetDRMass2() const;
   double GetDRMass3() const;
   double GetSMax() const; //The maximum enerfy is set by the user
   double GetSMin(); //The minimum energy to decay to 3 rest particles

   double Get_Sij_Min(int i, int j);
   double Get_Sij_Max(int i, int j);

   double Get_S23_Dalitz_Min_wrt_S13(double s, double s13);
   double Get_S23_Dalitz_Max_wrt_S13(double s, double s13);


private :
 
   double m_Dalitz[4];
};

#endif
