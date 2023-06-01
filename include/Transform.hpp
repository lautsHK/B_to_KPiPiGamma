//This is a header file for rotating and Lorentz boost

#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

class Transform
{
public:

   //Euler Rotaion 
   double Euler_Rotation_x( double angle1, double Cos_angle2, double angle3, double vector_x, double vector_y, double vector_z );
   double Euler_Rotation_y( double angle1, double Cos_angle2, double angle3, double vector_x, double vector_y, double vector_z );
   double Euler_Rotation_z( double angle1, double Cos_angle2, double angle3, double vector_x, double vector_y, double vector_z );

   //Lorentz Transformation
   double Lorentz_Boost_E( double beta, double energy, double momentum );
   double Lorentz_Boost_p( double beta, double energy, double momentum );

private:
};

#endif
