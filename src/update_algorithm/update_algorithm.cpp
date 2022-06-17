#include "update_algorithm.hpp"
#include <cmath>

void coulombic_force(site_positions* site_a, site_positions* site_b, double q_a, double q_b, site_forces* forces)
{
  //force constant
  double k_e = 8.988E9;
  //local variables for simplification, compiler will inline
  double x_a = site_a->x;
  double y_a = site_a->y;
  double z_a = site_a->z;

  double x_b = site_b->x;
  double y_b = site_b->y;
  double z_b = site_b->z;

  //forces
  double f_x = - k_e*q_a*q_b*(x_b - x_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 3/2);
  double f_y = - k_e*q_a*q_b*(y_b - y_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 3/2);
  double f_z = - k_e*q_a*q_b*(z_b - z_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 3/2);

  forces->fx = f_x;
  forces->fy = f_y;
  forces->fz = f_z;
}

void lj_force(site_positions* site_a, site_positions* site_b, double sigma, double epsilon, site_forces* forces)
{
  //local variables for simplification, compiler will inline
  double x_a = site_a->x;
  double y_a = site_a->y;
  double z_a = site_a->z;

  double x_b = site_b->x;
  double y_b = site_b->y;
  double z_b = site_b->z;

  //forces
  double f_x = - 4*epsilon*((12*pow(sigma, 12)*(x_b - x_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 7)) - (6*pow(sigma, 6)*(x_b - x_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 4)));
  double f_y = - 4*epsilon*((12*pow(sigma, 12)*(y_b - y_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 7)) - (6*pow(sigma, 6)*(y_b - y_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 4)));
  double f_z = - 4*epsilon*((12*pow(sigma, 12)*(z_b - z_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 7)) - (6*pow(sigma, 6)*(z_b - z_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 4)));

  forces->fx = f_x;
  forces->fy = f_y;
  forces->fz = f_z;
}
