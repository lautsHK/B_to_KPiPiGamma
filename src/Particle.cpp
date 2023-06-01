//including headers
#include "Particle.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cassert>
#include <map>
#include <string> 

Particle::Particle(int PDGID)
{
   m_PDGID = PDGID;

   //numbers are from PDG 2022
   m_Name_Map[111] = "pi0";
   //m_Mass_Map[111] = 0.1349768; //GeV
   m_Mass_Map[111] = 0.14; //GeV, MI value
   m_Width_Map[111] = 7.73*pow(10.0, -9.0); //GeV
   m_Charge_Map[111] = 0; 
   
   m_Name_Map[211] = "pi+";
   m_Mass_Map[211] = 0.13957039; //GeV
   m_Width_Map[211] = 0.0; //GeV
   m_Charge_Map[211] = 1; 

   m_Name_Map[113] = "rho0";
   //m_Mass_Map[113] = 0.77526; //GeV
   m_Mass_Map[113] = 0.775; //GeV, MI value
   //m_Width_Map[113] = 0.1491; //GeV
   m_Width_Map[113] = 0.149; //GeV, MI value
   m_Charge_Map[113] = 0; 
   
   m_Name_Map[213] = "rho+";
   m_Mass_Map[213] = 0.77526; //GeV
   m_Width_Map[213] = 0.1491; //GeV
   m_Charge_Map[213] = 1; 
   
   m_Name_Map[223] = "omega";
   m_Mass_Map[223] = 0.78266; //GeV
   m_Width_Map[223] = 0.00868; //GeV
   m_Charge_Map[223] = 0; 
   
   m_Name_Map[311] = "K0";
   //m_Mass_Map[311] = 0.497611; //GeV
   m_Mass_Map[311] = 0.498; //GeV, MI value
   m_Width_Map[311] = 0.0; //GeV
   m_Charge_Map[311] = 0; 

   m_Name_Map[321] = "K+";
   m_Mass_Map[321] = 0.493677; //GeV
   m_Width_Map[321] = 0.0; //GeV
   m_Charge_Map[321] = 1; 

   m_Name_Map[313] = "K*0"; // K*(892)
   //m_Mass_Map[313] = 0.89555; //GeV
   m_Mass_Map[313] = 0.892; //GeV, MI value
   //m_Width_Map[313] = 0.0473; //GeV
   m_Width_Map[313] = 0.05; //GeV, MI value
   m_Charge_Map[313] = 0; 

   m_Name_Map[323] = "K*+"; // K*(892)
   m_Mass_Map[323] = 0.89167; //GeV
   m_Width_Map[323] = 0.0514; //GeV
   m_Charge_Map[323] = 1; 
   
   m_Name_Map[9000311] = "kappa"; // kappa, K0*(700), also called K0*(800)
   m_Mass_Map[9000311] = 0.89555; //GeV
   m_Width_Map[9000311] = 0.468; //GeV
   m_Charge_Map[9000311] = 0; 

   m_Name_Map[10313] = "K_10"; // K1(1270)
   //m_Mass_Map[10313] = 1.253; //GeV
   m_Mass_Map[10313] = 1.272; //GeV, MI value
   m_Width_Map[10313] = 0.09; //GeV
   m_Charge_Map[10313] = 0; 

   m_Name_Map[10323] = "K_1+"; // K1(1270)+
   m_Mass_Map[10323] = 1.253; //GeV
   m_Width_Map[10323] = 0.09; //GeV
   m_Charge_Map[10323] = 1; 

   m_Name_Map[20313] = "K'_10"; // K1(1400)
   m_Mass_Map[20313] = 1.403; //GeV
   m_Width_Map[20313] = 0.174; //GeV
   m_Charge_Map[20313] = 0; 

   m_Name_Map[20323] = "K'_1+"; // K1(1400)+
   m_Mass_Map[20323] = 1.403; //GeV
   m_Width_Map[20323] = 0.174; //GeV
   m_Charge_Map[20323] = 1; 

   m_Name_Map[100313] = "K'*0"; // K*(1410)
   m_Mass_Map[100313] = 1.414; //GeV, MI value
   m_Width_Map[100313] = 0.232; //GeV
   m_Charge_Map[100313] = 0; 

   m_Name_Map[100323] = "K'*+"; // K*(1410)+
   m_Mass_Map[100323] = 1.414; //GeV
   m_Width_Map[100323] = 0.232; //GeV
   m_Charge_Map[100323] = 1; 

   m_Name_Map[10311] = "K_0*0"; // K*0(1430)
   m_Mass_Map[10311] = 1.425; //GeV
   m_Width_Map[10311] = 0.27; //GeV
   m_Charge_Map[10311] = 0; 

   m_Name_Map[10321] = "K_0*+"; // K*0(1430)+
   m_Mass_Map[10321] = 1.425; //GeV
   m_Width_Map[10321] = 0.27; //GeV
   m_Charge_Map[10321] = 1; 

   m_Name_Map[315] = "K_2*0"; // K*2(1430)
   m_Mass_Map[315] = 1.4324; //GeV
   m_Width_Map[315] = 0.109; //GeV
   m_Charge_Map[315] = 0; 

   m_Name_Map[325] = "K_2*+"; // K*2(1430)+
   m_Mass_Map[325] = 1.4273; //GeV
   m_Width_Map[325] = 0.100; //GeV
   m_Charge_Map[325] = 1; 

   m_Name_Map[30313] = "K''*0"; // K*(1680)
   m_Mass_Map[30313] = 1.718; //GeV
   m_Width_Map[30313] = 0.322; //GeV
   m_Charge_Map[30313] = 0; 

   m_Name_Map[30323] = "K''*+"; // K*(1680)+
   m_Mass_Map[30323] = 1.718; //GeV
   m_Width_Map[30323] = 0.322; //GeV
   m_Charge_Map[30323] = 1; 

   m_Name_Map[511] = "B0";
   m_Mass_Map[511] = 5.27966; //GeV
   m_Width_Map[511] = 4.3*pow(10.0, -16.0); //GeV
   m_Charge_Map[511] = 0; 
   
   m_Name_Map[521] = "B+";
   m_Mass_Map[521] = 5.27934; //GeV
   m_Width_Map[521] = 4.3*pow(10.0, -16.0); //GeV
   m_Charge_Map[521] = 0; 
   
}

std::string Particle::Get_Particle_Name()
{
   std::string name;
   name = m_Name_Map[m_PDGID];

   return name;
}

double Particle::Get_Particle_Mass()
{
   double mass;
   mass = m_Mass_Map[m_PDGID];

   return mass;
}

double Particle::Get_Particle_Width()
{
   double width;
   width = m_Width_Map[m_PDGID];

   return width;
}

int Particle::Get_Particle_Charge()
{
   int charge;
   charge = m_Charge_Map[m_PDGID];

   return charge;
}
