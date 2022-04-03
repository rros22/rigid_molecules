#ifndef QUATERNIONS_HPP
#define QUATERNIONS_HPP

#include <array>

class quaternion{

private:
  //scalar part
  double q0;

  //vector part
  std::array<double, 3> q;

  void set_angle(double theta);
  void set_axis(std::array<double, 3> &q, double theta);

  std::array<double, 3> normalise(std::array<double, 3> &q);

public:

  void set_quaternion(double theta, std::array<double, 3> &q);

  double get_angle();
  std::array<double, 3> get_axis();

  std::array<double, 4> get_parameters();

  double get_norm();
  


};

class rot_matrix{

private:

  std::array<double, 9> elements;

  void set_elements(quaternion Q);

public:

  rot_matrix(quaternion Q);

  std::array<double, 3> matrix_vector(std::array<double, 3> a);


};

#endif
