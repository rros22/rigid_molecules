#include "quaternions.hpp"
#include <cmath>

//setter functions
void quaternion::set_angle(double theta)
{
  this->q0 = cos(theta/2);
}

void quaternion::set_axis(std::array<double, 3> &q, double theta)
{
  (this->q)[0] = sin(theta/2)*q[0];
  (this->q)[1] = sin(theta/2)*q[1];
  (this->q)[2] = sin(theta/2)*q[2];
}

std::array<double, 3> quaternion::normalise(std::array<double, 3> &q){

  double norm = sqrt(pow(q[0], 2) + pow(q[1], 2) + pow(q[2], 2));
  return {q[0]/norm, q[1]/norm, q[2]/norm};
}

void quaternion::set_quaternion(double theta, std::array<double, 3> &q)
{
  //normalise direction
  std::array<double, 3> u = normalise(q);

  //set members
  set_axis(u, theta);
  set_angle(theta);
}

//getter functions
double quaternion::get_angle()
{
  return 2*acos(this->q0);
}

const std::array<double, 3>& quaternion::get_axis()
{
  return this->q;
}

std::array<double, 4> quaternion::get_parameters()
{
  return {q0, q[0], q[1], q[2]};
}

double quaternion::get_norm()
{
  return sqrt(pow(q0, 2) + pow(q[0], 2) + pow(q[1], 2) + pow(q[2], 2));
}



/*

  Rotation matrix class

*/

void rot_matrix::set_elements(quaternion Q){

  double q0 = Q.get_parameters()[0];
  double q1 = Q.get_parameters()[1];
  double q2 = Q.get_parameters()[2];
  double q3 = Q.get_parameters()[3];

  //set elements
  elements[0] = pow(q0, 2) + pow(q1, 2) - pow(q2, 2) - pow(q3, 2);
  elements[1] = 2*(q1*q2 - q0*q3);
  elements[2] = 2*(q1*q3 + q0*q2);
  elements[3] = 2*(q1*q2 + q0*q3);
  elements[4] = pow(q0, 2) - pow(q1, 2) + pow(q2, 2) - pow(q3, 2);
  elements[5] = 2*(q2*q3 - q0*q1);
  elements[6] = 2*(q1*q3 - q0*q2);
  elements[7] = 2*(q2*q3 + q0*q1);
  elements[8] = pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2);

}

rot_matrix::rot_matrix(quaternion Q){

  set_elements(Q);
}

//matrix-vector multiplication
std::array<double, 3> rot_matrix::matrix_vector(std::array<double, 3> a)
{

  std::array<double, 3> result;

  for (int i = 0; i < 3; i++){

    result[i] = elements[3*i]*a[0] + elements[1 + 3*i]*a[1] + elements[2 + 3*i]*a[2];
  }

  return result;
}
