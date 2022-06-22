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

double quaternion::get_norm()
{
  return sqrt(pow(q0, 2) + pow(q[0], 2) + pow(q[1], 2) + pow(q[2], 2));
}

void quaternion::transform_vector(double* input, double* offset, double* output)
{
    //compiler will inline
    double rotation_matrix[9];
    //set matrix elements
    rotation_matrix[0] = pow(q0, 2) + pow(q[0], 2) - pow(q[1], 2) - pow(q[2], 2);
    rotation_matrix[1] = 2*(q[0]*q[1] - q0*q[2]);
    rotation_matrix[2] = 2*(q[0]*q[2] + q0*q[1]);
    rotation_matrix[3] = 2*(q[0]*q[1] + q0*q[2]);
    rotation_matrix[4] = pow(q0, 2) - pow(q[0], 2) + pow(q[1], 2) - pow(q[2], 2);
    rotation_matrix[5] = 2*(q[1]*q[2] - q0*q[0]);
    rotation_matrix[6] = 2*(q[0]*q[2] - q0*q[1]);
    rotation_matrix[7] = 2*(q[1]*q[2] + q0*q[0]);
    rotation_matrix[8] = pow(q0, 2) - pow(q[0], 2) - pow(q[1], 2) + pow(q[2], 2);
    //compute rotated vector
    for (int i = 0; i < 3; i++)
    {
      output[i] = rotation_matrix[3*i]*input[0] + rotation_matrix[1 + 3*i]*input[1] +
                  rotation_matrix[2 + 3*i]*input[2];

      output[i] += offset[i];  
    }

}
