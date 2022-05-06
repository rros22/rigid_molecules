#ifndef SITES_HPP
#define SITES_HPP

#include <iostream>
#include <array>
#include <memory>


//kernel to any site class
class site
{
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
  //constructor
  site();
  //constructor initialise m
  site(double m): m(m) {}

  //default setters
  void set_global_coordinates(std::array<double, 3> X);
  void set_local_coordinates(std::array<double, 3> X_l);
  void set_velocity(std::array<double, 3> V);
  void set_force(std::array<double, 3> F);
  void set_mass(double m);
  void set_symbol(std::string symbol);

  //default getters
  const std::array<double, 3>& get_global_coordinates();
  const std::array<double, 3>& get_local_coordinates();
  const std::array<double, 3>& get_velocity();
  const std::array<double, 3>& get_force();
  double get_mass();
  const std::string& get_symbol();


  void add_forces(std::array<double, 3> F);
  void reset_forces();

};

/*

Derived classes

*/


class charge
{
private:
  //charge
  double q = 0;

public:
  //contains a site
  site kernel;
  //default constructor
  charge();
  //initialiser list constructor
  charge(double q, double m);

  //default setters
  void set_charge(double q);

  //default getters
  double get_charge();

  void set_parameters(double q, double m);
  void calculate_forces(std::shared_ptr<site> site){}

  //void calculate_forces(charge *charge);
};

class lj_site
{
private:
  //parameters
  double sigma;
  double epsilon;

public:
  //constains a site
  site kernel;
  //default constructor
  lj_site();
  
  //initialise constructor
  lj_site(double sigma, double epsilon, double m);

  //setter functions
  void set_sigma(double sigma);
  void set_epsilon(double epsilon);

  //getter functions
  double get_sigma();
  double get_epsilon();

  void set_parameters(double sigma, double epsilon, double m);
  void set_forces(){}

  void calculate_forces(std::shared_ptr<site> site){}

  //void calculate_forces(lj_site *lj_site);
};

#endif
