#include "sites.hpp"

/*

Site

*/

void site::set_symbol(std::string symbol){

  this->symbol = symbol;
}

void site::set_local_coordinates(std::array<double, 3> X_l){

  this->X_l = X_l;
}

void site::set_global_coordinates(std::array<double, 3> X){

  this->X = X;
}

std::string site::get_symbol(){

  return symbol;
}

std::array<double, 3> site::get_local_coordinates(){

  return this->X_l;
}

std::array<double, 3> site::get_global_coordinates(){

  return this->X;
}

double site::get_mass(){

  return this->m;
}

std::array<double, 3> site::get_forces(){

  return F;
}

void site::add_forces(std::array<double, 3> F){

  this->F[0] += F[0];
  this->F[1] += F[1];
  this->F[2] += F[2];

}

void site::reset_forces(){

  F[0] = 0;
  F[1] = 0;
  F[2] = 0;
}

/*

Charge site

*/
double charge::get_charge(){

  return q;
}

void charge::set_parameters(double q, double m){

  this->q = q;
  this->m = m;
}

void charge::calculate_forces(charge *charge){

  //define dielectric constant
  double k_e = 8.988E9;

  //obtain charge values
  double q_1 = this->q;
  double q_2 = charge->get_charge();

  //obtain global coordinates
  std::array<double, 3> X_1 = this->get_global_coordinates();
  std::array<double, 3> X_2 = charge->get_global_coordinates();

  //arrays to store site forces
  std::array<double, 3> F_1;
  std::array<double, 3> F_2;

  //calculate coulomb force
  F_1[0] = - k_e*q_1*q_2*(X_2[0] - X_1[0])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 3/2);
  F_1[1] = - k_e*q_1*q_2*(X_2[1] - X_1[1])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 3/2);
  F_1[2] = - k_e*q_1*q_2*(X_2[2] - X_1[2])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 3/2);

  F_2[0] = - F_1[0];
  F_2[1] = - F_1[1];
  F_2[2] = - F_1[2];

  //add forces to both charges
  this->add_forces(F_1);
  charge->add_forces(F_2);

}

/*

Lennard jones site

*/

void lj_site::set_parameters(double sigma, double epsilon, double m){

  this->sigma = sigma;
  this->epsilon = epsilon;
  this->m = m;
}

void lj_site::calculate_forces(lj_site *lj_site){

  //define lennard jones parameters
  double epsilon =  this->epsilon;
  double sigma = this->sigma;

  //obtain global coordinates
  std::array<double, 3> X_1 = this->get_global_coordinates();
  std::array<double, 3> X_2 = lj_site->get_global_coordinates();

  //arrays to store site forces
  std::array<double, 3> F_1;
  std::array<double, 3> F_2;

  //calculate lennard jones force
  F_1[0] = - 4*epsilon*((12*pow(sigma, 12)*(X_2[0] - X_1[0])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 7)) - (6*pow(sigma, 6)*(X_2[0] - X_1[0])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 4)));
  F_1[1] = - 4*epsilon*((12*pow(sigma, 12)*(X_2[1] - X_1[1])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 7)) - (6*pow(sigma, 6)*(X_2[1] - X_1[1])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 4)));
  F_1[2] = - 4*epsilon*((12*pow(sigma, 12)*(X_2[2] - X_1[2])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 7)) - (6*pow(sigma, 6)*(X_2[2] - X_1[2])/pow(pow(X_2[0] - X_1[0], 2) + pow(X_2[1] - X_1[1], 2) + pow(X_2[2] - X_1[2], 2), 4)));

  F_2[0] = - F_1[0];
  F_2[1] = - F_1[1];
  F_2[2] = - F_1[2];

  //add forces to both charges
  this->add_forces(F_1);
  lj_site->add_forces(F_2);
}
