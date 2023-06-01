//This is a header file for setting the Dalitz region
#ifndef KRES_TO_VP_HPP
#define KRES_TO_VP_HPP

#include "Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp"
#include "Kinematics.hpp"
#include "Particle.hpp"

#include <cmath>

class Kres_to_VP
{
public :
   //Actually K1270 and K1400 can bd combined into one function. But the argument would be difficult to identify. Espeicially there will be at least 2 more K_{res}

   //For K*pi decay, int decay = 0
   //For rhoK decay, int decay = 1

   double Get_theta_K() const;

   //the "mass" below is for sqrt(s)
   //and the "mass_V" is for sqrt(s_ij)
   //They are not the rest mass 

   double MatrixElement_S_K1270(double mass, double mass_V, int decay); 
   double MatrixElement_S_K1270(double mass, double mass_V, int decay, double damp_fac); //Without f2 = 3.0 assuption

   double MatrixElement_S_K1400(double mass, double mass_V, int decay); 
   double MatrixElement_S_K1400(double mass, double mass_V, int decay, double damp_fac);

   double MatrixElement_D_K1270(double mass, double mass_V, int decay); 
   double MatrixElement_D_K1270(double mass, double mass_V, int decay, double damp_fac); 

   double MatrixElement_D_K1400(double mass, double mass_V, int decay); 
   double MatrixElement_D_K1400(double mass, double mass_V, int decay, double damp_fac); 

private :

   double theta_K = M_PI/3.0;

};

#endif
