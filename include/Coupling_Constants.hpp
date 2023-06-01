//This is a header file for setting coupling constants
#ifndef COUPLING_CONSTANTS_HPP
#define COUPLING_CONSTANTS_HPP

#include "Kinematics.hpp"
#include "Particle.hpp"

class Coupling_Constants
{
public :

   // Some background information
   // These numbers are copied from a mathematic file
   // No calculation will be provided here

   Coupling_Constants(int charge);

   double Get_Coupling_Kst_to_K_Pi() const;
   double Get_Coupling_Rho_to_Pi_Pi() const;
   double Get_Coupling_K0st_to_K_Pi() const;
   double Get_Coupling_K2st_to_Kst_Pi() const;
   double Get_Coupling_K2st_to_rho_K() const;
   double Get_Coupling_KstP_to_Kst_Pi() const;
   double Get_Coupling_KstP_to_rho_K() const;
   double Get_Coupling_KstPP_to_Kst_Pi() const;
   double Get_Coupling_KstPP_to_rho_K() const;
   double Get_Coupling_K1_to_kappa2KPi_Pi() const;

   //Clebsch - Gordan coefficients
   double Get_CG_Kst_Pi(); 
   double Get_CG_Rho_K();

   //0 or 1, depends on wherther the kaon is charged 
   double Get_C13(); //There will not be s13 contribution for charged kaon

private :
 
   double m_g_Kst_K_Pi = 5.67689; //Kst is K*(892)
   double m_g_rho_Pi_Pi = 5.97994;
   double m_g_K0st_K_Pi = 5.28597; //K0st is K*_{0}(1430). This number, actually hasn't been used in the this MC generator.
   double m_g_K2st_Kst_Pi = 20.; //K2st is K*_{2}(1430), denoted as K*2
   double m_g_K2st_rho_K = 20.; //K2st is K*_{2}(1430)
   double m_g_KstP_Kst_Pi = 9.0; //KstP is K*_{1}(1410), denoted as K'*.
   double m_g_KstP_rho_K = 0.1; //KstP is K*_{1}(1410)
   double m_g_KstPP_Kst_Pi = 9.0; //KstP is K*_{1}(1680), denoted as K''*.
   double m_g_KstPP_rho_K = 0.1; //KstP is K*_{1}(1680)
   double m_g_K1_kappa2KPi_Pi = 43.066; //K1 is only for K1(1270), and kappa is \kappa(800), also called K*_{0}(700)

   int m_charge;
};
#endif
