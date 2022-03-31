#ifndef SITES_HPP
#define SITES_HPP

#include <iostream>
#include <array>

class site{

protected:

  //global position
  std::array<double, 3> X;

  //local position (X - X(CoM))
  std::array<double, 3> X_l;

  //velocity
  std::array<double, 3> V;

  //forces
  std::array<double, 3> F;

  //name
  std::string symbol;

public:

  void set_symbol(std::string symbol);
  void set_local_coordinates(std::array<double, 3> X_l);
  void set_global_coordinates(std::array<double, 3> X);

  std::string get_symbol();
  std::array<double, 3> get_local_coordinates();
  std::array<double, 3> get_global_coordinates();

};

class atom: public site{

  //mass
  double m;

};

class charge: public site{

  //charge
  double q;

  //mass
  double m;

public:

  void set_parameters(double q, double m);

};

class lj_site: public site{

  //parameters
  double sigma;
  double epsilon;

  //mass
  double m;

public:

  void set_parameters(double sigma, double epsilon, double m);

};

#endif
