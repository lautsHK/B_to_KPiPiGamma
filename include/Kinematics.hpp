//This is a header file for calculating energies and momentums for 2 bodies decay, both center-of-momentum frame and decayed particle frame

#ifndef KINEMATICS_HPP 
#define KINEMATICS_HPP 

class Kinematics
{
public :
  
   double M2ab_CoP_Momentum(double mass_M, double mass_a, double mass_b);
   double M2ab_CoP_Energy_a(double mass_M, double mass_a, double mass_b);

   double M2ab_Momentum_a_wrt_b(double mass_M, double mass_a, double mass_b);
   double M2ab_Energy_a_wrt_b(double mass_M, double mass_a, double mass_b);
   

private :


};

#endif
