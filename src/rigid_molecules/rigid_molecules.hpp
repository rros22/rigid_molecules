#ifndef RIGID_MOLECULES_HPP
#define RIGID_MOLECULES_HPP

#include <array>

#include "../sites/sites.hpp"

class rigid_molecule{

  double CoM;

  //moments of inertia, wrt body frame (Ixx, Iyy, Izz)
  std::array<double, 3> I;

};

class water: public rigid_molecule{


};

#endif
