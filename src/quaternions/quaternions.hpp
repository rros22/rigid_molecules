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

};

struct rot_matrix
{

  double elements[9];

  //create matrix with quaternion defined components
  rot_matrix(quaternion Q);
  //set matrix elements with quaternion
  void set_elements(quaternion Q);
  //return array o input rotated coordinates
  std::array<double, 3> matrix_vector(std::array<double, 3> a);
  //overload operaator '*' to do the same as the above
  std::array<double, 3> operator * (const std::array<double, 3>& coordinate)
  {
    std::array<double, 3> result;

    for (int i = 0; i < 3; i++){

      result[i] = elements[3*i]*coordinate[0] + elements[1 + 3*i]*coordinate[1] + elements[2 + 3*i]*coordinate[2];
    }

    return result;
  }
};

#endif
