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

  //mass
  double m;

  //name
  std::string symbol;

public:

  void set_symbol(std::string symbol);
  void set_local_coordinates(std::array<double, 3> X_l);
  void set_global_coordinates(std::array<double, 3> X);
  virtual void set_forces() = 0;

  std::string get_symbol();
  std::array<double, 3> get_local_coordinates();
  std::array<double, 3> get_global_coordinates();
  double get_mass();

};

class atom: public site{

};

class charge: public site{

  //charge
  double q;

public:

  void set_parameters(double q, double m);
  void set_forces(){}

};

class lj_site: public site{

  //parameters
  double sigma;
  double epsilon;

public:

  void set_parameters(double sigma, double epsilon, double m);
  void set_forces(){}
};

#endif
