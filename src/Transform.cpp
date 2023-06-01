//including headers
#include "Transform.hpp"

//including c++ libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

double Transform::Euler_Rotation_x( double angle1, double Cos_angle2, double angle3, double vector_x, double vector_y, double vector_z )
{
   double rotated_x;

   double c1 = cos(angle1); 
   double c2 = Cos_angle2; 
   double c3 = cos(angle3); 

   double s1 = sin(angle1); 
   double s2 = sqrt(1. - pow(c2, 2.)); 
   double s3 = sin(angle3); 

   double R_xx = c1*c3 - c2*s3*s1;
   double R_xy = c1*s3 + c2*c1*s1;
   double R_xz = s1*s2;

   rotated_x = R_xx * vector_x + R_xy * vector_y + R_xz * vector_z;

   return rotated_x;
}

double Transform::Euler_Rotation_y( double angle1, double Cos_angle2, double angle3, double vector_x, double vector_y, double vector_z )
{
   double rotated_y;

   double c1 = cos(angle1); 
   double c2 = Cos_angle2; 
   double c3 = cos(angle3); 

   double s1 = sin(angle1); 
   double s2 = sqrt(1. - pow(c2, 2.)); 
   double s3 = sin(angle3); 

   double R_yx = -s1*c3 - c2*s3*c1;
   double R_yy = -s1*s3 + c2*c3*c1;
   double R_yz = c1*s2;

   rotated_y = R_yx * vector_x + R_yy * vector_y + R_yz * vector_z;

   return rotated_y;
}

double Transform::Euler_Rotation_z( double angle1, double Cos_angle2, double angle3, double vector_x, double vector_y, double vector_z )
{
   double rotated_z;

   double c1 = cos(angle1); 
   double c2 = Cos_angle2; 
   double c3 = cos(angle3); 

   double s1 = sin(angle1); 
   double s2 = sqrt(1. - pow(c2, 2.)); 
   double s3 = sin(angle3); 

   double R_zx = s2*s3;
   double R_zy = -s2*c3;
   double R_zz = c2;

   rotated_z = R_zx * vector_x + R_zy * vector_y + R_zz * vector_z;

   return rotated_z;
}

double Transform::Lorentz_Boost_E( double beta, double energy, double momentum )
{
   double gamma = 1./ sqrt( 1. - pow(beta, 2.0) );

   double transformed_E = gamma * energy + (-gamma*beta) * momentum;

   return transformed_E;
}

double Transform::Lorentz_Boost_p( double beta, double energy, double momentum )
{
   double gamma = 1./ sqrt( 1. - pow(beta, 2.0) );

   double transformed_p = (-gamma*beta) * energy + gamma * momentum;

   return transformed_p;
}
