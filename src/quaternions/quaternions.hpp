#ifndef QUATERNIONS_HPP
#define QUATERNIONS_HPP

#include <array>

struct quaternion
{
  double q0;   //scalar part
  double q[3]; //vector part

  //member functions
  quaternion(){};
  quaternion(double theta, std::array<double, 3> &q);

  void set_angle(double theta);
  void set_axis(std::array<double, 3> &q, double theta);
  std::array<double, 3> normalise(std::array<double, 3> &q);
  //setter functions
  void set_quaternion(double theta, std::array<double, 3> &q);
  //getter functions
  double get_angle();
  double get_norm();
  //calculate orientation. Compact the matrix struct definition
  void transform_vector(double* input, double* offset, double* output);

};



#endif
