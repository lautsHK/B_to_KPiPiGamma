//include the header files
#include "RandomGen.hpp"
#include "Dalitz_Region.hpp"
#include "Functions.hpp"
#include "Kinematics.hpp"
#include "ThreeBodies_Kinematics.hpp"
#include "Coupling_Constants.hpp"
#include "Particle.hpp"
#include "Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp"
#include "Quark_Pair_Creation_Model/Kres_to_VP.hpp"
#include "Form_Factors.hpp"
#include "Helicity_Amplitude.hpp"
#include "Matrix_Elements.hpp"
#include "Transform.hpp"

//include the c++ libraries
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <complex>
#include <vector>
#include <algorithm>


int main ( int argc, char* argv[] )
{
   //Test 1
   Coupling_Constants CCs(1); //1 means charged Kaon

   double gKstKPi = CCs.Get_Coupling_Kst_to_K_Pi();
   double gRhoPiPi = CCs.Get_Coupling_Rho_to_Pi_Pi();
 
   std::cout << "\n" << "Using the value same as the Mathematica file,\n" << "gKstKPi = " << gKstKPi << "\n";
   std::cout << "gRhoPiPi = " << gRhoPiPi << "\n\n"; 

   //Test 2
   QPCM_A_to_VplusP Kres;
   std::cout << "QPCM values : \n";
   double value_I0 = Kres.I0(1.27, 0.9, 0.14); //K1 to Kst + Pi
   double value_I1 = Kres.I1(1.27, 0.9, 0.14);
 
   std::cout << "I0(M=1.27, MV=0.9, MP=0.14, RA=RV=RP=2.5) = " << value_I0 << "\n"; 
   std::cout << "I1(M=1.27, MV=0.9, MP=0.14, RA=RV=RP=2.5) = " << value_I1 << "\n"; 
   
   double _MSpqc = Kres.MSpqc(1.27, 0.9, 0.14); 
   double _MDpqc = Kres.MDpqc(1.27, 0.9, 0.14); 

   std::cout << "MSpqc(M=1.27, MV=0.9, MP=0.14) = " << _MSpqc << "\n"; 
   std::cout << "MDpqc(M=1.27, MV=0.9, MP=0.14) = " << _MDpqc << "\n\n"; 
   
   //Test 3
   Kres_to_VP K1270_to_VP;
   std::cout << "Matrix Elements : \n";
   std::cout << "For Kres is K1270 : \n"; 
  
   double MSpqcKst = K1270_to_VP.MatrixElement_S_K1270(1.27, 0.9, 0);
   std::cout << "MSpqcKst = " << MSpqcKst << "\n";

   double MSpqcrho = K1270_to_VP.MatrixElement_S_K1270(1.27, 0.7, 1);
   std::cout << "MSpqcrho = " << MSpqcrho << "\n";

   double MDpqcKst = K1270_to_VP.MatrixElement_D_K1270(1.27, 0.9, 0);
   std::cout << "MDpqcKst = " << MDpqcKst << "\n";

   double MDpqcrho = K1270_to_VP.MatrixElement_D_K1270(1.27, 0.7, 1);
   std::cout << "MDpqcrho = " << MDpqcrho << "\n";

   //Test 4
   Kres_to_VP K1400_to_VP;
   std::cout << "For Kres is K1400 : \n"; 
  
   double MSpqcKst2 = K1400_to_VP.MatrixElement_S_K1400(1.40, 0.9, 0);
   std::cout << "MSpqcKst = " << MSpqcKst2 << "\n";

   double MSpqcrho2 = K1400_to_VP.MatrixElement_S_K1400(1.40, 0.7, 1);
   std::cout << "MSpqcrho = " << MSpqcrho2 << "\n";

   double MDpqcKst2 = K1400_to_VP.MatrixElement_D_K1400(1.40, 0.9, 0);
   std::cout << "MDpqcKst = " << MDpqcKst2 << "\n";

   double MDpqcrho2 = K1400_to_VP.MatrixElement_D_K1400(1.40, 0.7, 1);
   std::cout << "MDpqcrho = " << MDpqcrho2 << "\n\n";

   //Test 5
   Form_Factors FFs;
   std::cout << "Form Factors : \n";

   std::cout << "For Kres is K1270 : \n";
 
   double AKst_1270 = FFs.A_K1270(pow(1.27,2.), pow(0.9,2.), 0);
   std::cout << "AKst = " << AKst_1270 << "\n";

   double Arho_1270 = FFs.A_K1270(pow(1.27,2.), pow(0.7,2.), 1);
   std::cout << "Arho = " << Arho_1270 << "\n";

   double BKst_1270 = FFs.B_K1270(pow(1.27,2.), pow(0.9,2.), 0);
   std::cout << "BKst = " << BKst_1270 << "\n";

   double Brho_1270 = FFs.B_K1270(pow(1.27,2.), pow(0.7,2.), 1);
   std::cout << "Brho = " << Brho_1270 << "\n";

   std::cout << "For Kres is K1400 : \n";
 
   double AKst_1400 = FFs.A_K1400(pow(1.40,2.), pow(0.9,2.), 0);
   std::cout << "AKst = " << AKst_1400 << "\n";

   double Arho_1400 = FFs.A_K1400(pow(1.40,2.), pow(0.7,2.), 1);
   std::cout << "Arho = " << Arho_1400 << "\n";

   double BKst_1400 = FFs.B_K1400(pow(1.40,2.), pow(0.9,2.), 0);
   std::cout << "BKst = " << BKst_1400 << "\n";

   double Brho_1400 = FFs.B_K1400(pow(1.40,2.), pow(0.7,2.), 1);
   std::cout << "Brho = " << Brho_1400 << "\n\n";

   //Test 6
   Helicity_Amplitude HAs;
   std::cout << "Helicity Amplitudes : \n";

   //double Helicity_Amplitude::J2_1270( double s, double s13, double s23, int k_charge)

   double _J2_1270 = HAs.J2_1270(pow(1.27,2.), pow(0.9,2.), pow(0.9,2.), 1); // 1 means charged Kaon
   std::cout << "J2(1270) = " << _J2_1270 << "\n";

   double _JxJC_1270 = HAs.JcrossJC_1270(pow(1.27,2.), pow(0.9,2.), pow(0.9,2.), 1); // 1 means charged Kaon
   std::cout << "JxJC(1270) = " << _JxJC_1270 << "\n";

   double _J2_1400 = HAs.J2_1400(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 1); // 1 means charged Kaon
   std::cout << "J2(1400) = " << _J2_1400 << "\n";

   double _JxJC_1400 = HAs.JcrossJC_1400(pow(1.27,2.), pow(0.9,2.), pow(0.9,2.), 1); // 1 means charged Kaon
   std::cout << "JxJC(1400) = " << _JxJC_1400 << "\n";

   double _J2_total = HAs.J2_Total(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 1, 1.36603, -0.366025); // 1 means charged Kaon
   std::cout << "J2(1270+1400) = " << _J2_total << "\n";

   double _J2_total_0 = HAs.J2_Total(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0, 1.36603, -0.366025); // 1 means charged Kaon
   std::cout << "J2(1270+1400) of neutral Kaon = " << _J2_total_0 << "\n";

   double _JxJC_total = HAs.JcrossJC_Total(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 1, 1.36603, -0.366025); // 1 means charged Kaon
   std::cout << "JxJC(1270+1400) = " << _JxJC_total << "\n\n";

   //Test 7
   Matrix_Elements MEs;
   std::cout << "Matrix Elements : \n";

   std::complex<double> _Jx;
   _Jx = MEs.Jx_Two_Resonances(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025);
   std::cout << "Jx(1270+1400) = " << _Jx << "\n";

   std::complex<double> _Jx5;
   _Jx5 = MEs.Jx_All_Resonances(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025, 0.0, 0.0, 0.0);
   std::cout << "Jx(1270+1400+3Res) = " << _Jx5 << "\n";

   std::complex<double> _Jy;
   _Jy = MEs.Jy_Two_Resonances(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025);
   std::cout << "Jy(1270+1400) = " << _Jy << "\n";

   std::complex<double> _Jy5;
   _Jy5 = MEs.Jy_All_Resonances(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025, 0.0, 0.0, 0.0);
   std::cout << "Jy(1270+1400+3Res) = " << _Jy5 << "\n";

   std::complex<double> _Jz;
   _Jz = MEs.Jz_Two_Resonances(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025);
   std::cout << "Jz(1270+1400) = " << _Jz << "\n";

   std::complex<double> _Jz5;
   _Jz5 = MEs.Jz_All_Resonances(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025, 0.0, 0.0, 0.0);
   std::cout << "Jz(1270+1400+3Res) = " << _Jz5 << "\n\n";

   std::complex<double> J2_z = _Jx * std::conj(_Jx) + _Jy * std::conj(_Jy) + _Jz * std::conj(_Jz);
   std::cout << "J2(1270+1400) = " << J2_z << ", with cos_theta = phi1 = 0.5\n";

   double _J2;
   _J2 = MEs.J2(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025, 0.0, 0.0, 0.0); 
   std::cout << "J2(1270+1400) using ME = " << _J2 << ", with cos_theta = phi1 = 0.5\n\n";

   double _nJxJC;
   _nJxJC = MEs.Im_N_dot_J_cross_Jc(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025, 0.0, 0.0, 0.0); 
   std::cout << "JxJC(1270+1400) using ME = " << _nJxJC << ", with cos_theta = phi1 = 0.5\n\n";

   double _MEL2;
   _MEL2 = MEs.ME2_LH(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025, 0.0, 0.0, 0.0); 
   std::cout << "Left Polarized Matrix ELement = " << _MEL2 << ", with cos_theta = phi1 = 0.5\n\n";

   double _MER2;
   _MER2 = MEs.ME2_RH(pow(1.4,2.), pow(0.9,2.), pow(0.9,2.), 0.5, 0.5, 1, 1.36603, -0.366025, 0.0, 0.0, 0.0); 
   std::cout << "Right Polarized Matrix ELement = " << _MER2 << ", with cos_theta = phi1 = 0.5\n\n";

   return 0;
} 
