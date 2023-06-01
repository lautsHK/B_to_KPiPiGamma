//including headers
#include "Coupling_Constants.hpp"
#include "Kinematics.hpp"
#include "Particle.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>

Coupling_Constants::Coupling_Constants(int charge)
{
   m_charge = charge;
}

double Coupling_Constants::Get_Coupling_Kst_to_K_Pi() const
{
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   double mass_Kst_sq = pow(mass_Kst, 2.);
 
   Kinematics* kine = new Kinematics;
   double P = kine->M2ab_CoP_Momentum(mass_Kst, mass_pi , mass_K);
   double P3 = pow(P, 3.);
  
   double C;
   C = sqrt( 3. * 2. * M_PI * mass_Kst_sq * width_Kst / P3);

   delete Pion;
   delete Kaon;
   delete KStar;
   delete kine;

   return C;
}

double Coupling_Constants::Get_Coupling_Rho_to_Pi_Pi() const
{
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();
   double mass_rho_sq = pow(mass_rho, 2.);
 
   Kinematics* kine = new Kinematics;
   double P = kine->M2ab_CoP_Momentum(mass_rho, mass_pi , mass_pi);
   double P3 = pow(P, 3.);

   double C;
   C = sqrt( 3. * 2. * M_PI * mass_rho_sq * width_rho / P3);

   return C;
}

double Coupling_Constants::Get_Coupling_K0st_to_K_Pi() const
{
   double C;
   C = m_g_K0st_K_Pi;

   return C;
}

double Coupling_Constants::Get_Coupling_K2st_to_Kst_Pi() const
{
   double C;
   C = m_g_K2st_Kst_Pi;

   return C;
}

double Coupling_Constants::Get_Coupling_K2st_to_rho_K() const
{
   double C;
   C = m_g_K2st_rho_K;

   return C;
}

double Coupling_Constants::Get_Coupling_KstP_to_Kst_Pi() const
{
   double C;
   C = m_g_KstP_Kst_Pi;

   return C;
}

double Coupling_Constants::Get_Coupling_KstP_to_rho_K() const
{
   double C;
   C = m_g_KstP_rho_K;

   return C;
}

double Coupling_Constants::Get_Coupling_KstPP_to_Kst_Pi() const
{
   double C;
   C = m_g_KstPP_Kst_Pi;

   return C;
}

double Coupling_Constants::Get_Coupling_KstPP_to_rho_K() const
{
   double C;
   C = m_g_KstPP_rho_K;

   return C;
}

double Coupling_Constants::Get_Coupling_K1_to_kappa2KPi_Pi() const
{
   double C;
   C = m_g_K1_kappa2KPi_Pi;

   return C;
}

double Coupling_Constants::Get_CG_Kst_Pi()
{
   double CG;

   if (m_charge == 0)
   {
      CG = sqrt(2.)/3.;
   }

   if (m_charge == 1)
   {
      CG = -2./3.;
   }

   return CG;
}

double Coupling_Constants::Get_CG_Rho_K()
{
   double CG;

   if (m_charge == 0)
   {
      CG = 1./sqrt(3.);
   }

   if (m_charge == 1)
   {
      CG = -1./sqrt(6.);
   }

   return CG;
}

double Coupling_Constants::Get_C13()
{
   double C13;

   if (m_charge == 0)
   {
      C13 = 1.;
   }

   if (m_charge == 1)
   {
      C13 = 0.;
   }

   return C13;
}
