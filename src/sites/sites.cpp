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

/*

Charge site

*/

void charge::set_parameters(double q, double m){

  this->q = q;
  this->m = m;
}

/*

Lennard jones site

*/

void lj_site::set_parameters(double sigma, double epsilon, double m){

  this->sigma = sigma;
  this->epsilon = epsilon;
  this->m = m;
}
