#ifndef RIGID_MOLECULES_HPP
#define RIGID_MOLECULES_HPP

#include <array>
#include <vector>
#include <iostream>
#include <memory>

#include "../sites/sites.hpp"
#include "../quaternions/quaternions.hpp"

//abstract base class

class rigid_molecule
{
private:
  //global coordinates of the center of mass
  std::array<double, 3> CoM;
  //body frame orientation wrt global frame
  quaternion Q;
  //mass of molecule
  double m;
  //moments of inertia, wrt body frame (Ixx, Iyy, Izz)
  std::array<double, 3> I;
  //center of mass force, in global coordinates
  std::array<double, 3> F;
  //torque, about principal axes i.e. local coordinates
  std::array<double, 3> T;
  //name
  std::string symbol;

public:
  //setter functions
  void set_mass(double m);
  void set_mom_inertia(std::array<double, 3> I);
  void set_position(std::array<double, 3> &CoM);
  void set_orientation(quaternion &Q);
  void set_CoM_force(std::array<double, 3> F);
  void set_torque();
  void set_symbol(std::string symbol);

  //getter functions
  const std::array<double, 3>& get_position();
  const quaternion& get_orientation();
  double get_mass();
  const std::array<double, 3>& get_mom_inertia();
  const std::array<double, 3>& get_CoM_force();
  const std::array<double, 3>& get_torque();
  const std::string& get_symbol();

  void reset_CoM_force();

};

class H2O
{
private:
  lj_site O1;
  charge H1;
  charge H2;
  charge q;

public:
  //common to any molecule
  rigid_molecule kernel;

  H2O();

  ~H2O() {}

  H2O(std::array<double, 3> CoM, quaternion Q);

  //site coordinates
  void set_global_coordinates();
  //forces
  void set_forces(std::shared_ptr<rigid_molecule>){};

};


#endif
