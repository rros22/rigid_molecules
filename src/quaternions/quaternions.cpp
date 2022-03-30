#include "quaternions.hpp"
#include <cmath>


void quaternion::set_angle(double theta){

  this->q0 = cos(theta/2);
}

void quaternion::set_axis(std::array<double, 3> &q, double theta){

  (this->q)[0] = sin(theta/2)*q[0];
  (this->q)[1] = sin(theta/2)*q[1];
  (this->q)[2] = sin(theta/2)*q[2];
}

std::array<double, 3> quaternion::normalise(std::array<double, 3> &q){

  double norm = sqrt(pow(q[0], 2) + pow(q[1], 2) + pow(q[2], 2));
  return {q[0]/norm, q[1]/norm, q[1]/norm};
}

void quaternion::set_quaternion(double theta, std::array<double, 3> &q){

  //normalise direction
  std::array<double, 3> u = normalise(q);

  //set members
  set_axis(u, theta);
  set_angle(theta);
}

double quaternion::get_angle(){

  return 2*acos(this->q0);
}

std::array<double, 3> quaternion::get_axis(){

  return this->q;
}
