//This is a header file for storing particle information
#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <map>
#include <string> 

class Particle
{
public :

   // Some background information
   // These numbers are copied from PDG in 2022
   Particle(int PDGID);

   std::string Get_Particle_Name();
   
   double Get_Particle_Mass();

   double Get_Particle_Width();

   int Get_Particle_Charge();

private :

   int m_PDGID;

   std::map<int, std::string> m_Name_Map;
   std::map<int, double> m_Mass_Map;
   std::map<int, double> m_Width_Map;  
   std::map<int, int> m_Charge_Map;  
};

#endif
