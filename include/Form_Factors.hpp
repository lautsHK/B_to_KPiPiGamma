//This is a header file for calculaating Form Factos, which will have contribution to the S-wave and D-wave
#ifndef FORM_FACTORS_HPP
#define FORM_FACTORS_HPP

#include "Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp"
#include "Quark_Pair_Creation_Model/Kres_to_VP.hpp"
#include "Kinematics.hpp"
#include "Particle.hpp"

#include <cmath>

class Form_Factors
{
public :

   double A_K1270(double s, double s_ij, int decay); 
   double B_K1270(double s, double s_ij, int decay); 

   double A_K1400(double s, double s_ij, int decay); 
   double B_K1400(double s, double s_ij, int decay);    

private :

};

#endif
