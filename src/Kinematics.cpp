//including headers
#include "Kinematics.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>

double Kinematics::M2ab_CoP_Momentum(double mass_M, double mass_a, double mass_b)
{
   double M2 = pow(mass_M, 2.0);

   double sum_of_decay_masses = mass_a + mass_b;
   assert(mass_M >= sum_of_decay_masses);
   double sum2 = pow(sum_of_decay_masses, 2.0);

   double diff_of_decay_masses = mass_a - mass_b;
   double diff2 = pow(diff_of_decay_masses, 2.0);

   double p;
   p = sqrt( (M2-sum2)*(M2-diff2) ) / (2.0*mass_M);

   return p;
}


double Kinematics::M2ab_CoP_Energy_a(double mass_M, double mass_a, double mass_b)
{
   double M2 = pow(mass_M, 2.0);
   double Ma2 = pow(mass_a, 2.0);
   double Mb2 = pow(mass_b, 2.0);

   double Ea;
   Ea = (M2 + Ma2 - Mb2) / (2.0*mass_M);

   return Ea;
}

double Kinematics::M2ab_Momentum_a_wrt_b(double mass_M, double mass_a, double mass_b)
{
   double M2 = pow(mass_M, 2.0);

   double sum_of_decay_masses = mass_a + mass_b;
   assert(mass_M >= sum_of_decay_masses);
   double sum2 = pow(sum_of_decay_masses, 2.0);

   double diff_of_decay_masses = mass_a - mass_b;
   double diff2 = pow(diff_of_decay_masses, 2.0);

   double pa;
   pa = sqrt( (M2-sum2)*(M2-diff2) ) / (2.0*mass_b);

   return pa;
}

double Kinematics::M2ab_Energy_a_wrt_b(double mass_M, double mass_a, double mass_b)
{
   double M2 = pow(mass_M, 2.0);
   double Ma2 = pow(mass_a, 2.0);
   double Mb2 = pow(mass_b, 2.0);

   double Ea;
   Ea = (M2 - Ma2 - Mb2) / (2.0*mass_b);

   return Ea;
}



