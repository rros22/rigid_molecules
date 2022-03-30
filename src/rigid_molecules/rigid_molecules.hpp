#ifndef RIGID_MOLECULES_HPP
#define RIGID_MOLECULES_HPP

#include <array>

#include "../sites/sites.hpp"
#include "../quaternions/quaternions.hpp"

class rigid_molecule{

protected:

  //global coordinates of the center of mass
  std::array<double, 3> CoM;

  //body frame orientation wrt global frame
  quaternion Q;

  //moments of inertia, wrt body frame (Ixx, Iyy, Izz)
  std::array<double, 3> I;

public:

  void set_position(std::array<double, 3> &CoM);
  void set_orientation(quaternion &Q);

};

class water: public rigid_molecule{

private:

  //interaction sites (atoms and dummy charges)
  lj_site O1;
  charge H1;
  charge H2;

  charge q1;

  //local coordinates of sites
  std::array<double, 3> O1_X;
  std::array<double, 3> H1_X;
  std::array<double, 3> H2_X;
  std::array<double, 3> q1_X;

public:

  water();

  void set_global_coordinates();



};

#endif
