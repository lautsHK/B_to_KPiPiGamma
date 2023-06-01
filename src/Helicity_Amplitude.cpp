//include headers
#include "Helicity_Amplitude.hpp"
#include "Kinematics.hpp"
#include "ThreeBodies_Kinematics.hpp"
#include "Form_Factors.hpp"
#include "Particle.hpp"
#include "Coupling_Constants.hpp"
#include "Functions.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <complex>

std::complex<double> Helicity_Amplitude::C1_1270(double s, double s13, double s23, int k_charge)
{
   //***Set up the PiPiKaon decay***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.
   double E1 = PiPiK->Get_CoP_Energy1(s, s13, s23); 
   double E2 = PiPiK->Get_CoP_Energy2(s, s13, s23); 
   double p1p2 = PiPiK->p1_4Dot_p2(s, s13, s23);

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst);
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Form Factors***
   Form_Factors* F = new Form_Factors;

   double A_Kst_s23;
   A_Kst_s23 = F->A_K1270(s, s23, 0); //last interger is decay mode.
   double A_Kst_s13;
   A_Kst_s13 = F->A_K1270(s, s13, 0); //last interger is decay mode.
   double B_Kst_s23;
   B_Kst_s23 = F->B_K1270(s, s23, 0); //last interger is decay mode.

   double A_rho_s12;
   A_rho_s12 = F->A_K1270(s, s12, 1); //last interger is decay mode.
   double B_rho_s12;
   B_rho_s12 = F->B_K1270(s, s12, 1); //last interger is decay mode.

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   //***Calculation of C1_1270***
   //The K* Part
   std::complex<double> C1_Kst;

   /*breaking down the details*/
   double AK_BW23_C = 1. + (m2_sq-m3_sq)/mKst_sq; 
   double BK_BW23_C = -1. * AK_BW23_C * ( sqrt(s)*E1 - m1_sq ) + 2.*p1p2;  

   std::complex<double> C1_Kst_23;
   C1_Kst_23 = ( A_Kst_s23 * AK_BW23_C + B_Kst_s23 * BK_BW23_C ) * BW_Kst_s23;  

   std::complex<double> C1_Kst_13;
   C1_Kst_13 = -2. * A_Kst_s13 * BW_Kst_s13;

   C1_Kst = CG_KstPi * gKstKPi * (C1_Kst_23 + C1_Kst_13 * c13);

   //The Rho Part
   std::complex<double> C1_rho;

   double BR_BW12_C = -1. * sqrt(s) * (E1 - E2);

   C1_rho = CG_rhoK * gRhoPiPi * ( A_rho_s12 + B_rho_s12 * BR_BW12_C) * BW_rho_s12;

   //Combine
   std::complex<double> C1;
   C1 = C1_Kst + C1_rho;


   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete Rho;
   delete FN;
   delete F; 
   delete CCs;

   return C1;
}


std::complex<double> Helicity_Amplitude::C2_1270(double s, double s13, double s23, int k_charge)
{
   //***Set up the PiPiKaon decay***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.
   double E1 = PiPiK->Get_CoP_Energy1(s, s13, s23); 
   double E2 = PiPiK->Get_CoP_Energy2(s, s13, s23); 
   double p1p2 = PiPiK->p1_4Dot_p2(s, s13, s23);

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Form Factors***
   Form_Factors* F = new Form_Factors;

   double A_Kst_s23;
   A_Kst_s23 = F->A_K1270(s, s23, 0); //last interger is decay mode.
   double A_Kst_s13;
   A_Kst_s13 = F->A_K1270(s, s13, 0); //last interger is decay mode.
   double B_Kst_s13;
   B_Kst_s13 = F->B_K1270(s, s13, 0); //last interger is decay mode.

   double A_rho_s12;
   A_rho_s12 = F->A_K1270(s, s12, 1); //last interger is decay mode.
   double B_rho_s12;
   B_rho_s12 = F->B_K1270(s, s12, 1); //last interger is decay mode.

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   //***Calculation of C2_1270***
   //The K* Part
   std::complex<double> C2_Kst;

   /*breaking down the details*/
   double AK_BW13_C = 1. + (m1_sq-m3_sq)/mKst_sq; 
   double BK_BW13_C = -1. * AK_BW13_C * ( sqrt(s)*E2 - m2_sq ) + 2.*p1p2;  

   std::complex<double> C2_Kst_13;
   C2_Kst_13 = ( A_Kst_s13 * AK_BW13_C + B_Kst_s13 * BK_BW13_C ) * BW_Kst_s13;  

   std::complex<double> C2_Kst_23;
   C2_Kst_23 = -2. * A_Kst_s23 * BW_Kst_s23;

   C2_Kst = CG_KstPi * gKstKPi * (C2_Kst_23 + C2_Kst_13 * c13);

   //The Rho Part
   std::complex<double> C2_rho;

   double BR_BW12_C = sqrt(s) * (E1 - E2);

   C2_rho = CG_rhoK * gRhoPiPi * ( A_rho_s12 + B_rho_s12 * BR_BW12_C) * BW_rho_s12;

   //Combine
   std::complex<double> C2;
   C2 = C2_Kst + C2_rho;


   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete Rho;
   delete FN;
   delete F; 
   delete CCs;

   return C2;
}


std::complex<double> Helicity_Amplitude::C1_1400(double s, double s13, double s23, int k_charge)
{
   //***Set up the PiPiKaon decay***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.
   double E1 = PiPiK->Get_CoP_Energy1(s, s13, s23); 
   double E2 = PiPiK->Get_CoP_Energy2(s, s13, s23); 
   double p1p2 = PiPiK->p1_4Dot_p2(s, s13, s23);

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Form Factors***
   Form_Factors* F = new Form_Factors;

   double A_Kst_s23;
   A_Kst_s23 = F->A_K1400(s, s23, 0); //last interger is decay mode.
   double A_Kst_s13;
   A_Kst_s13 = F->A_K1400(s, s13, 0); //last interger is decay mode.
   double B_Kst_s23;
   B_Kst_s23 = F->B_K1400(s, s23, 0); //last interger is decay mode.

   double A_rho_s12;
   A_rho_s12 = F->A_K1400(s, s12, 1); //last interger is decay mode.
   double B_rho_s12;
   B_rho_s12 = F->B_K1400(s, s12, 1); //last interger is decay mode.

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   //***Calculation of C1_1400***
   //The K* Part
   std::complex<double> C1_Kst;

   /*breaking down the details*/
   double AK_BW23_C = 1. + (m2_sq-m3_sq)/mKst_sq; 
   double BK_BW23_C = -1. * AK_BW23_C * ( sqrt(s)*E1 - m1_sq ) + 2.*p1p2;  

   std::complex<double> C1_Kst_23;
   C1_Kst_23 = ( A_Kst_s23 * AK_BW23_C + B_Kst_s23 * BK_BW23_C ) * BW_Kst_s23;  

   std::complex<double> C1_Kst_13;
   C1_Kst_13 = -2. * A_Kst_s13 * BW_Kst_s13;

   C1_Kst = CG_KstPi * gKstKPi * (C1_Kst_23 + C1_Kst_13 * c13);

   //The Rho Part
   std::complex<double> C1_rho;

   double BR_BW12_C = -1. * sqrt(s) * (E1 - E2);

   C1_rho = CG_rhoK * gRhoPiPi * ( A_rho_s12 + B_rho_s12 * BR_BW12_C) * BW_rho_s12;

   //Combine
   std::complex<double> C1;
   C1 = C1_Kst + C1_rho;


   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete Rho;
   delete FN;
   delete F; 
   delete CCs;

   return C1;
}


std::complex<double> Helicity_Amplitude::C2_1400(double s, double s13, double s23, int k_charge)
{
   //***Set up the PiPiKaon decay***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.
   double E1 = PiPiK->Get_CoP_Energy1(s, s13, s23); 
   double E2 = PiPiK->Get_CoP_Energy2(s, s13, s23); 
   double p1p2 = PiPiK->p1_4Dot_p2(s, s13, s23);

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Form Factors***
   Form_Factors* F = new Form_Factors;

   double A_Kst_s23;
   A_Kst_s23 = F->A_K1400(s, s23, 0); //last interger is decay mode.
   double A_Kst_s13;
   A_Kst_s13 = F->A_K1400(s, s13, 0); //last interger is decay mode.
   double B_Kst_s13;
   B_Kst_s13 = F->B_K1400(s, s13, 0); //last interger is decay mode.

   double A_rho_s12;
   A_rho_s12 = F->A_K1400(s, s12, 1); //last interger is decay mode.
   double B_rho_s12;
   B_rho_s12 = F->B_K1400(s, s12, 1); //last interger is decay mode.

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   //***Calculation of C2_1400***
   //The K* Part
   std::complex<double> C2_Kst;

   /*breaking down the details*/
   double AK_BW13_C = 1. + (m1_sq-m3_sq)/mKst_sq; 
   double BK_BW13_C = -1. * AK_BW13_C * ( sqrt(s)*E2 - m2_sq ) + 2.*p1p2;  

   std::complex<double> C2_Kst_13;
   C2_Kst_13 = ( A_Kst_s13 * AK_BW13_C + B_Kst_s13 * BK_BW13_C ) * BW_Kst_s13;  

   std::complex<double> C2_Kst_23;
   C2_Kst_23 = -2. * A_Kst_s23 * BW_Kst_s23;

   C2_Kst = CG_KstPi * gKstKPi * (C2_Kst_23 + C2_Kst_13 * c13);
    
   //The Rho Part
   std::complex<double> C2_rho;

   double BR_BW12_C = sqrt(s) * (E1 - E2);

   C2_rho = CG_rhoK * gRhoPiPi * ( A_rho_s12 + B_rho_s12 * BR_BW12_C) * BW_rho_s12;

   //Combine
   std::complex<double> C2;
   C2 = C2_Kst + C2_rho;


   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete Rho;
   delete FN;
   delete F; 
   delete CCs;

   return C2;
}

std::complex<double> Helicity_Amplitude::C3_1410( double s, double s13, double s23, int k_charge)
{
   //***Set up the PiPiKaon decay for s12 calculation***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics to calculate s12***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   //double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst);
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double gKstP_to_Kst_Pi = CCs->Get_Coupling_KstP_to_Kst_Pi();
   double gKstP_to_rho_K = CCs->Get_Coupling_KstP_to_rho_K();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   //***Calculation of C3_1410***
   //The K* Part
   Particle* KStarPrime = new Particle(100313); //100313 is the PDGID of K'*0 (100323 for K'*+)
   double mass_KstP = KStarPrime->Get_Particle_Mass();

   std::complex<double> Im_i(0.0 , 1.0);

   std::complex<double> C3_Kst;
   C3_Kst = CG_KstPi * (2. * Im_i * mass_KstP * gKstP_to_Kst_Pi * gKstKPi) * (BW_Kst_s23 + BW_Kst_s13*c13); 

   std::complex<double> C3_rho;
   C3_rho = CG_rhoK * (2. * Im_i * mass_KstP * gKstP_to_rho_K * gRhoPiPi) * (BW_rho_s12); 

   //Combine
   std::complex<double> C3;
   C3 = C3_Kst + C3_rho;


   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete KStarPrime;
   delete Rho;
   delete FN;
   delete CCs;

   return C3;
}

std::complex<double> Helicity_Amplitude::C3_1680( double s, double s13, double s23, int k_charge)
{
   //***Set up the PiPiKaon decay for s12 calculation***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics to calculate s12***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   //double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst);
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double gKstP_to_Kst_Pi = CCs->Get_Coupling_KstP_to_Kst_Pi(); //Value of KstPP should be used here. 
   double gKstP_to_rho_K = CCs->Get_Coupling_KstP_to_rho_K(); //Value of KstPP hould be used here.
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   //***Calculation of C3_1680***
   //The K* Part
   Particle* KStarPrimePrime = new Particle(30313); //30313 is the PDGID of K''*0 (30323 for K''*+)
   double mass_KstPP = KStarPrimePrime->Get_Particle_Mass();

   std::complex<double> Im_i(0.0 , 1.0);

   std::complex<double> C3_Kst;
   C3_Kst = CG_KstPi * (2. * Im_i * mass_KstPP * gKstP_to_Kst_Pi * gKstKPi) * (BW_Kst_s23 + BW_Kst_s13*c13); 

   std::complex<double> C3_rho;
   C3_rho = CG_rhoK * (2. * Im_i * mass_KstPP * gKstP_to_rho_K * gRhoPiPi) * (BW_rho_s12); 

   //Combine
   std::complex<double> C3;
   C3 = C3_Kst + C3_rho;


   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete KStarPrimePrime;
   delete Rho;
   delete FN;
   delete CCs;

   return C3;
}

std::complex<double> Helicity_Amplitude::C1_1430( double s, double s13, double s23, double cos_theta, double phi1, double phi2, int k_charge)
{
   //***Set up the PiPiKaon decay for s12 calculation***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.
   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23);
   double sin_delta = sin( phi2 - phi1 );
   double sin_theta = sqrt(1. - pow(cos_theta, 2.)); //theta's range is from 0 to M_PI. So we don't need to worry about the sign of sin_theta

   double p1p2_sinD_cosT = p1 * p2 * sin_delta * cos_theta; 
   double p1_sinT_cosP1 = p1 * sin_theta * cos(phi1); 
   double p2_sinT_cosP2 = p2 * sin_theta * cos(phi2); 

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   //double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst);
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double gK2st_to_Kst_Pi = CCs->Get_Coupling_K2st_to_Kst_Pi();
   double gK2st_to_rho_K = CCs->Get_Coupling_K2st_to_rho_K();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   std::complex<double> Im_i(0.0 , 1.0);
   Particle* K2Star = new Particle(315); //315 is the PDGID of K*2 (325 for K*2+)
   double mass_K2st = K2Star->Get_Particle_Mass();
   std::complex<double> Coupling_Kst;
   Coupling_Kst = 2. * Im_i * mass_K2st * gK2st_to_Kst_Pi * gKstKPi;
   std::complex<double> Coupling_rho;
   Coupling_rho = 2. * Im_i * mass_K2st * gK2st_to_rho_K * gRhoPiPi;

   //***Calculation of C1_1430***
   //The K2*(1430) Part

   std::complex<double> C1_Kst;
   C1_Kst = CG_KstPi * (-1.*Coupling_Kst) * p1p2_sinD_cosT * BW_Kst_s23; 

   std::complex<double> C1_rho;
   C1_rho = CG_rhoK * Coupling_rho * p1p2_sinD_cosT * BW_rho_s12;

   std::complex<double> C1;
   C1 = C1_Kst + C1_rho;

   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete K2Star;
   delete Rho;
   delete FN;
   delete CCs;

   return C1;
}


std::complex<double> Helicity_Amplitude::C2_1430( double s, double s13, double s23, double cos_theta, double phi1, double phi2, int k_charge)
{
   //***Set up the PiPiKaon decay for s12 calculation***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.
   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23);
   double sin_delta = sin( phi2 - phi1 );
   double sin_theta = sqrt(1. - pow(cos_theta, 2.));

   double p1p2_sinD_cosT = p1 * p2 * sin_delta * cos_theta; 
   double p1_sinT_cosP1 = p1 * sin_theta * cos(phi1); 
   double p2_sinT_cosP2 = p2 * sin_theta * cos(phi2); 

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   //double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst);
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double gK2st_to_Kst_Pi = CCs->Get_Coupling_K2st_to_Kst_Pi();
   double gK2st_to_rho_K = CCs->Get_Coupling_K2st_to_rho_K();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   std::complex<double> Im_i(0.0 , 1.0);
   Particle* K2Star = new Particle(315); //315 is the PDGID of K*2 (325 for K*2+)
   double mass_K2st = K2Star->Get_Particle_Mass();
   std::complex<double> Coupling_Kst;
   Coupling_Kst = 2. * Im_i * mass_K2st * gK2st_to_Kst_Pi * gKstKPi;
   std::complex<double> Coupling_rho;
   Coupling_rho = 2. * Im_i * mass_K2st * gK2st_to_rho_K * gRhoPiPi;

   //***Calculation of C2_1430***
   //The K2*(1430) Part

   std::complex<double> C2_Kst;
   C2_Kst = CG_KstPi * Coupling_Kst * p1p2_sinD_cosT * BW_Kst_s13 * c13; 

   std::complex<double> C2_rho;
   C2_rho = CG_rhoK * (-1.*Coupling_rho) * p1p2_sinD_cosT * BW_rho_s12;

   std::complex<double> C2;
   C2 = C2_Kst + C2_rho;

   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete K2Star;
   delete Rho;
   delete FN;
   delete CCs;

   return C2;
}


std::complex<double> Helicity_Amplitude::C3_1430( double s, double s13, double s23, double cos_theta, double phi1, double phi2, int k_charge)
{
   //***Set up the PiPiKaon decay for s12 calculation***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   double m1 = mass_pi; 
   double m2 = mass_pi; 
   double m3 = mass_K; 

   double m1_sq = pow(m1, 2.);
   double m2_sq = pow(m2, 2.);
   double m3_sq = pow(m3, 2.);

   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(m1, m2, m3);

   double s12 = PiPiK->Get_s12(s, s13, s23); //We will use it for rho-meson.
   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23);
   double sin_delta = sin( phi2 - phi1 );
   double sin_theta = sqrt(1. - pow(cos_theta, 2.));

   double p1p2_sinD_cosT = p1 * p2 * sin_delta * cos_theta; 
   double p1_sinT_cosP1 = p1 * sin_theta * cos(phi1); 
   double p2_sinT_cosP2 = p2 * sin_theta * cos(phi2); 

   //***Breit-Wigner***
   Particle* KStar = new Particle(313); //313 is the PDGID of K*0 (323 for K*+)
   double mass_Kst = KStar->Get_Particle_Mass();
   double width_Kst = KStar->Get_Particle_Width();
   //double mKst_sq = pow(mass_Kst, 2.);

   Particle* Rho = new Particle(113); //113 is the PDGID of rho0 (213 for rho+)
   double mass_rho = Rho->Get_Particle_Mass();
   double width_rho = Rho->Get_Particle_Width();

   Functions* FN = new Functions; 
   std::complex<double> BW_Kst_s23 = FN->BWPropagator(s23 , mass_Kst , width_Kst); //(s, mass, width)
   std::complex<double> BW_Kst_s13 = FN->BWPropagator(s13 , mass_Kst , width_Kst);
   std::complex<double> BW_rho_s12 = FN->BWPropagator(s12 , mass_rho , width_rho);

   //***Coupling Constants***
   Coupling_Constants* CCs = new Coupling_Constants(k_charge); //0 for neutral kaon, 1 for charged kaon
   double gKstKPi = CCs->Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs->Get_Coupling_Rho_to_Pi_Pi();
   double gK2st_to_Kst_Pi = CCs->Get_Coupling_K2st_to_Kst_Pi();
   double gK2st_to_rho_K = CCs->Get_Coupling_K2st_to_rho_K();
   double CG_KstPi = CCs->Get_CG_Kst_Pi();
   double CG_rhoK = CCs->Get_CG_Rho_K();
   double c13 = CCs->Get_C13();

   std::complex<double> Im_i(0.0 , 1.0);
   Particle* K2Star = new Particle(315); //315 is the PDGID of K*2 (325 for K*2+)
   double mass_K2st = K2Star->Get_Particle_Mass();
   std::complex<double> Coupling_Kst;
   Coupling_Kst = 2. * Im_i * mass_K2st * gK2st_to_Kst_Pi * gKstKPi;
   std::complex<double> Coupling_rho;
   Coupling_rho = 2. * Im_i * mass_K2st * gK2st_to_rho_K * gRhoPiPi;

   //***Calculation of C3_1430***
   //The K2*(1430) Part

   std::complex<double> C3_Kst;
   C3_Kst = CG_KstPi * Coupling_Kst * ( (p2_sinT_cosP2 * BW_Kst_s13 * c13) + (p1_sinT_cosP1 * BW_Kst_s23) ); 

   std::complex<double> C3_rho;
   C3_rho = CG_rhoK * (-1.*Coupling_rho) * (p1_sinT_cosP1 + p2_sinT_cosP2) * BW_rho_s12;

   std::complex<double> C3;
   C3 = C3_Kst + C3_rho;

   delete Pion;
   delete Kaon;
   delete PiPiK;
   delete KStar;
   delete K2Star;
   delete Rho;
   delete FN;
   delete CCs;

   return C3;
}
/////////////////////////////////////////////////////////////////////////
//                                                                     //
//         The function below is for validation purpose only.          //
//        They contain only K1(1270) and K1(1400) information.         //
//      For 5 resonances, please go to the Matrix Element method.      //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

double Helicity_Amplitude::J2_1270( double s, double s13, double s23, int k_charge)
{
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   std::complex<double> C1;
   C1 = C1_1270(s, s13, s23, k_charge);
   double C1_sq;
   C1_sq = pow( abs(C1), 2.0 );

   std::complex<double> C2;
   C2 = C2_1270(s, s13, s23, k_charge);
   double C2_sq;
   C2_sq = pow( abs(C2), 2.0 );

   std::complex<double> C1C2_z = C1 * std::conj(C2);
   std::complex<double> C2C1_z = C2 * std::conj(C1);

   double C1C2 = std::real(C1C2_z);
   double C2C1 = std::real(C2C1_z);

   double J2;
   J2 = C1_sq * p1_sq + C2_sq * p2_sq - (C1C2 + C2C1) * p1p2;  

   delete Pion;
   delete Kaon;
   delete PiPiK;

   return J2;
}

double Helicity_Amplitude::J2_1400( double s, double s13, double s23, int k_charge)
{
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   std::complex<double> C1;
   C1 = C1_1400(s, s13, s23, k_charge);
   double C1_sq;
   C1_sq = pow( abs(C1), 2.0 );

   std::complex<double> C2;
   C2 = C2_1400(s, s13, s23, k_charge);
   double C2_sq;
   C2_sq = pow( abs(C2), 2.0 );

   std::complex<double> C1C2_z = C1 * std::conj(C2);
   std::complex<double> C2C1_z = C2 * std::conj(C1);

   double C1C2 = std::real(C1C2_z);
   double C2C1 = std::real(C2C1_z);

   double J2;
   J2 = C1_sq * p1_sq + C2_sq * p2_sq - (C1C2 + C2C1) * p1p2;  

   delete PiPiK;
   delete Pion;
   delete Kaon;

   return J2;
}

double Helicity_Amplitude::J2_Total( double s, double s13, double s23, int k_charge, double T1K1270, double T1K1400)
{
   //***Particles***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   Particle* K1270 = new Particle(10313); //311 is the PDGID of K0
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //311 is the PDGID of K0
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   //***Breit-Wigner for K1270 and K1400***
   Functions* FN = new Functions; 
   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)

   //***C1_res and C2_res***
   //The C1_res part
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = C1_1400(s, s13, s23, k_charge);
 
   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400;
   
   double C1res_sq;
   C1res_sq = pow( abs(C1_res), 2.0 );
   
   //The C2_res part
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = C2_1400(s, s13, s23, k_charge);
 
   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400;

   double C2res_sq;
   C2res_sq = pow( abs(C2_res), 2.0 );

   //The combined one
   std::complex<double> C1C2_z = C1_res * std::conj(C2_res);
   std::complex<double> C2C1_z = C2_res * std::conj(C1_res);

   double C1C2 = std::real(C1C2_z);
   double C2C1 = std::real(C2C1_z);


   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1 = PiPiK->Get_CoP_Momentum1(s, s13, s23);
   double p1_sq = pow(p1, 2.);
 
   double p2 = PiPiK->Get_CoP_Momentum2(s, s13, s23); 
   double p2_sq = pow(p2, 2.);
   double p1p2 = PiPiK->p1_3Dot_p2(s, s13, s23);

   double J2_res;
   J2_res = C1res_sq * p1_sq + C2res_sq * p2_sq - p1p2 * (C1C2 + C2C1);

   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete FN;

   return J2_res;
}


double Helicity_Amplitude::JcrossJC_1270( double s, double s13, double s23, int k_charge)
{
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1_x_p2 = PiPiK->p1_Cross_p2(s, s13, s23);

   std::complex<double> C1;
   C1 = C1_1270(s, s13, s23, k_charge);

   std::complex<double> C2;
   C2 = C2_1270(s, s13, s23, k_charge);

   std::complex<double> C1C2_z = C1 * std::conj(C2);
   std::complex<double> C2C1_z = C2 * std::conj(C1);

   double C1C2 = std::imag(C1C2_z);
   double C2C1 = std::imag(C2C1_z);

   double JxJC;
   JxJC = p1_x_p2 * ( -C1C2 + C2C1 );  

   delete Pion;
   delete Kaon;
   delete PiPiK;

   return JxJC;
}


double Helicity_Amplitude::JcrossJC_1400( double s, double s13, double s23, int k_charge)
{
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1_x_p2 = PiPiK->p1_Cross_p2(s, s13, s23);

   std::complex<double> C1;
   C1 = C1_1400(s, s13, s23, k_charge);

   std::complex<double> C2;
   C2 = C2_1400(s, s13, s23, k_charge);

   std::complex<double> C1C2_z = C1 * std::conj(C2);
   std::complex<double> C2C1_z = C2 * std::conj(C1);

   double C1C2 = std::imag(C1C2_z);
   double C2C1 = std::imag(C2C1_z);

   double JxJC;
   JxJC = p1_x_p2 * ( -C1C2 + C2C1 );  

   delete Pion;
   delete Kaon;
   delete PiPiK;

   return JxJC;
}


double Helicity_Amplitude::JcrossJC_Total( double s, double s13, double s23, int k_charge, double T1K1270, double T1K1400)
{
   //***Particles***
   Particle* Pion = new Particle(111); //111 is the PDGID of pi0
   double mass_pi = Pion->Get_Particle_Mass();

   Particle* Kaon = new Particle(311); //311 is the PDGID of K0
   double mass_K = Kaon->Get_Particle_Mass();

   Particle* K1270 = new Particle(10313); //311 is the PDGID of K0
   double mass_K1270 = K1270->Get_Particle_Mass();
   double width_K1270 = K1270->Get_Particle_Width();

   Particle* K1400 = new Particle(20313); //311 is the PDGID of K0
   double mass_K1400 = K1400->Get_Particle_Mass();
   double width_K1400 = K1400->Get_Particle_Width();

   //***Breit-Wigner for K1270 and K1400***
   Functions* FN = new Functions; 
   std::complex<double> BW_K1400 = FN->BWPropagator(s, mass_K1400 , width_K1400); //(s, mass, width)
   std::complex<double> BW_K1270 = FN->BWPropagator(s, mass_K1270 , width_K1270); //(s, mass, width)

   //***C1_res and C2_res***
   //The C1_res part
   std::complex<double> C1_res;

   std::complex<double> C1_K1270;
   C1_K1270 = C1_1270(s, s13, s23, k_charge);

   std::complex<double> C1_K1400;
   C1_K1400 = C1_1400(s, s13, s23, k_charge);
 
   C1_res = T1K1270 * BW_K1270 * C1_K1270 + T1K1400 * BW_K1400 * C1_K1400;
   
   
   //The C2_res part
   std::complex<double> C2_res;

   std::complex<double> C2_K1270;
   C2_K1270 = C2_1270(s, s13, s23, k_charge);

   std::complex<double> C2_K1400;
   C2_K1400 = C2_1400(s, s13, s23, k_charge);
 
   C2_res = T1K1270 * BW_K1270 * C2_K1270 + T1K1400 * BW_K1400 * C2_K1400;


   //The combined one
   std::complex<double> C1C2_z = C1_res * std::conj(C2_res);
   std::complex<double> C2C1_z = C2_res * std::conj(C1_res);

   double C1C2 = std::imag(C1C2_z);
   double C2C1 = std::imag(C2C1_z);


   //***Kinematics***
   ThreeBodies_Kinematics* PiPiK = new ThreeBodies_Kinematics(mass_pi, mass_pi, mass_K);

   double p1_x_p2 = PiPiK->p1_Cross_p2(s, s13, s23);

   double JxJC_res;
   JxJC_res = p1_x_p2 * ( -C1C2 + C2C1 );  

   delete Pion;
   delete Kaon;
   delete K1270;
   delete K1400;
   delete FN;

   return JxJC_res;
}

