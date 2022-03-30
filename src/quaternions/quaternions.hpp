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

};

#endif
