#ifndef RIGID_MOLECULES_HPP
#define RIGID_MOLECULES_HPP

#include <array>
#include <vector>
#include <iostream>
#include <memory>

#include "../sites/sites.hpp"
#include "../quaternions/quaternions.hpp"

//abstract base class

class rigid_molecule{

protected:

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

  //list of interaction sites to interaction sites
  std::vector<std::shared_ptr<site>> sites;

  //name
  std::string symbol;

  //set molecule mass (cannot vary)
  void set_mass();

  //private methods to be accessed by constructor or member methods
  void set_symbol(std::string symbol);


public:

  virtual ~rigid_molecule(){}

  void set_position(std::array<double, 3> &CoM);
  void set_orientation(quaternion &Q);
  double get_mass();
  std::string get_symbol();

  void set_global_coordinates();
  std::vector<std::shared_ptr<site>> return_sites_list();
  std::array<double, 3> return_coordinates_site(int i);
  void reset_CoM_force();
  void set_CoM_force();

  //methods to be implemented by derived classes
  virtual void set_forces(rigid_molecule *molecule) = 0;

};

class H2O: public rigid_molecule{

  //interaction sites allocated by constructor on sites vector (atoms and dummy charges)

public:

  H2O();

  ~H2O() {}

  H2O(std::array<double, 3> CoM, quaternion Q);

  //forces
  void set_forces(rigid_molecule *molecule){};

};

class N2: public rigid_molecule{



public:

  N2();

  ~N2() {};

  N2(std::array<double, 3> CoM, quaternion Q);

  //forces
  void set_forces(rigid_molecule *molecule){}

};

#endif
