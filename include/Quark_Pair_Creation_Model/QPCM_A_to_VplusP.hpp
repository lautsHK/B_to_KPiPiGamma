//This is a header file for setting the Dalitz region
#ifndef QPCM_A_TO_VPLUSP_HPP
#define QPCM_A_TO_VPLUSP_HPP

#include "Kinematics.hpp"

class QPCM_A_to_VplusP
{
public :

  double Get_gammaPQC() const;
  double Get_f2() const;
  double Get_R_SHO() const;

  double I0(double mass_A, double mass_V, double mass_P);
  double I0(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P);

  double I1(double mass_A, double mass_V, double mass_P);
  double I1(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P);

  //needa check the meaning of f2
  double MSpqc(double mass_A, double mass_V, double mass_P); //S-wave
  double MSpqc(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P, double damp_fac); //S-wave

  double MDpqc(double mass_A, double mass_V, double mass_P); //D-wave
  double MDpqc(double mass_A, double mass_V, double mass_P, double radius_A, double radius_V, double radius_P, double damp_fac); //D-wave
  

private :

  double gammaPQC = 4.0; //pair creation constant, phenomenological dimesnionless
  double f2 = 3.0; //damping factor
  double R_SHO = 2.5; //Harmonic Oscillator radius of meson wave function
  /*
 * double Delta_E = energy shift between ground and first excited state.
 * Delta_E = 2./(mass_quark * R_SHO^2.0);
 * Hence R_SHO can be calculated 
 * With modification according to Godfrey and Isgur model,
 * for most L = 0,1 cases, R_SHO = 2.5 is a good approximation.
 * For more information, see Chapter 3.2.2 in Andrey Tayduganov's doctor thesis. 
 */
 
};

#endif
