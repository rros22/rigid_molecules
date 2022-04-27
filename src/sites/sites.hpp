#ifndef SITES_HPP
#define SITES_HPP

#include <iostream>
#include <array>
#include <memory>

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

  //private methods to be accessed by constructor or member methods


public:

  site();

  //constructor initialise m
  site(double m): m(m) {}

  void set_symbol(std::string symbol);

  void set_local_coordinates(std::array<double, 3> X_l);
  void set_global_coordinates(std::array<double, 3> X);

  std::string get_symbol();
  std::array<double, 3> get_local_coordinates();
  std::array<double, 3> get_global_coordinates();
  double get_mass();
  std::array<double, 3> get_forces();

  void add_forces(std::array<double, 3> F);
  void reset_forces();

  virtual void calculate_forces(std::shared_ptr<site> site) = 0;

};

class atom: public site{

};

class charge: public site{

  //charge
  double q;

public:

  charge();

  charge(double q, double m): site(m), q(q) {}

  double get_charge();

  void set_parameters(double q, double m);

  void calculate_forces(std::shared_ptr<site> site){}

  //void calculate_forces(charge *charge);

};

class lj_site: public site{

  //parameters
  double sigma;
  double epsilon;

public:

  lj_site();

  lj_site(double sigma, double epsilon, double m): site(m), sigma(sigma), epsilon(epsilon) {}

  void set_parameters(double sigma, double epsilon, double m);
  void set_forces(){}

  void calculate_forces(std::shared_ptr<site> site){}

  //void calculate_forces(lj_site *lj_site);
};

#endif
