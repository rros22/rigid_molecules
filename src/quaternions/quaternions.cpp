#include "quaternions.hpp"
#include <cmath>

//setter functions
void quaternion::set_angle(double theta)
{
  this->q0 = cos(theta/2);
}

void quaternion::set_axis(double theta, double q[3])
{
  (this->q)[0] = sin(theta/2)*q[0];
  (this->q)[1] = sin(theta/2)*q[1];
  (this->q)[2] = sin(theta/2)*q[2];
}

void quaternion::normalise(double q[3])
{
  double norm = sqrt(pow(q[0], 2) + pow(q[1], 2) + pow(q[2], 2));
  q[0] /= norm;
  q[1] /= norm;
  q[2] /= norm;
}

void quaternion::set_quaternion(double theta, double q[3])
{
  //normalise direction
  normalise(q);
  //set members
  set_axis(theta, q);
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

double quaternion::get_axis(unsigned i)
{
  if (i == 0 || i == 1 || i == 2)
  {
    return q[i]/sqrt(pow(q[0], 2) + pow(q[1], 2) + pow(q[2], 2));
  }

  else
  {
    return 0;
  }
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

void quaternion::transform_vector_invert(double* input, double* offset, double* output)
{
    //compiler will inline
    double rotation_matrix[9];
    //set matrix elements
    rotation_matrix[0] = pow(q0, 2) + pow(q[0], 2) - pow(q[1], 2) - pow(q[2], 2);
    rotation_matrix[1] = 2*(q[0]*q[1] + q0*q[2]);
    rotation_matrix[2] = 2*(q[0]*q[2] - q0*q[1]);
    rotation_matrix[3] = 2*(q[0]*q[1] - q0*q[2]);
    rotation_matrix[4] = pow(q0, 2) - pow(q[0], 2) + pow(q[1], 2) - pow(q[2], 2);
    rotation_matrix[5] = 2*(q[1]*q[2] + q0*q[0]);
    rotation_matrix[6] = 2*(q[0]*q[2] + q0*q[1]);
    rotation_matrix[7] = 2*(q[1]*q[2] - q0*q[0]);
    rotation_matrix[8] = pow(q0, 2) - pow(q[0], 2) - pow(q[1], 2) + pow(q[2], 2);
    //compute rotated vector
    for (int i = 0; i < 3; i++)
    {
      input[i] -= offset[i];

      output[i] = rotation_matrix[3*i]*input[0] + rotation_matrix[1 + 3*i]*input[1] +
                  rotation_matrix[2 + 3*i]*input[2];
    }
}

//operators
quaternion quaternion::operator+(quaternion quat)
{
  quaternion result;
  result.q0 = this->q0 + quat.q0;
  (result.q)[0] = (this->q)[0] + (quat.q)[0];
  (result.q)[1] = (this->q)[1] + (quat.q)[1];
  (result.q)[2] = (this->q)[2] + (quat.q)[2];
  return result;
}

quaternion quaternion::operator-(quaternion quat)
{
  quaternion result;
  result.q0 = this->q0 - quat.q0;
  (result.q)[0] = (this->q)[0] - (quat.q)[0];
  (result.q)[1] = (this->q)[1] - (quat.q)[1];
  (result.q)[2] = (this->q)[2] - (quat.q)[2];
  return result;
}

quaternion quaternion::multiply(double x)
{
  quaternion result;
  result.q0 = this->q0 * x;
  (result.q)[0] = (this->q)[0] * x;
  (result.q)[1] = (this->q)[1] * x;
  (result.q)[2] = (this->q)[2] * x;
  return result;
}
