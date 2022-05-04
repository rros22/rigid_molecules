#ifndef QUATERNIONS_HPP
#define QUATERNIONS_HPP

#include <array>

class quaternion
{
private:
  //scalar part
  double q0;
  //vector part
  std::array<double, 3> q;
  //constructor called functions
  void set_angle(double theta);
  void set_axis(std::array<double, 3> &q, double theta);
  std::array<double, 3> normalise(std::array<double, 3> &q);

public:
  //setter functions
  void set_quaternion(double theta, std::array<double, 3> &q);
  //getter functions
  double get_angle();
  const std::array<double, 3>& get_axis();
  std::array<double, 4> get_parameters();
  double get_norm();

};

class rot_matrix
{
private:
  //store matrix elements
  std::array<double, 9> elements;

public:
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
