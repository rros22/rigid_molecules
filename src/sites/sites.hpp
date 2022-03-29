#ifndef SITES_HPP
#define SITES_HPP

#include <iostream>
#include <array>

class site{

public:

  //position
  std::array<double, 3> X;

  //velocity
  std::array<double, 3> V;

  //forces
  std::array<double, 3> F;

};

class atom: public site{

  //mass
  double m;

};

class charge: public site{

  //charge
  double q;

};

#endif
