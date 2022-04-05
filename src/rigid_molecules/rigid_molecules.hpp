#ifndef RIGID_MOLECULES_HPP
#define RIGID_MOLECULES_HPP

#include <array>
#include <vector>
#include <iostream>

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

  //list of pointers to interaction sites
  std::vector<site*> sites;

  //set molecule mass (cannot vary)
  void set_mass();


public:

  virtual ~rigid_molecule(){}

  void set_position(std::array<double, 3> &CoM);
  void set_orientation(quaternion &Q);
  double get_mass();

  void set_global_coordinates();
  std::vector<site*> return_sites_list();
  std::array<double, 3> return_coordinates_site(int i);

  //methods to be implemented by derived classes


};

class H2O: public rigid_molecule{

private:

  //interaction sites (atoms and dummy charges)
  lj_site O_1;
  charge H_1;
  charge H_2;

  charge q_1;

public:

  H2O();

  ~H2O() {}

  H2O(std::array<double, 3> CoM, quaternion Q);

};

class N2: public rigid_molecule{

  //interaction sites
  lj_site N_1;
  lj_site N_2;

  charge q_1;

  //local coordinates of sites
  std::array<double, 3> N1_X;
  std::array<double, 3> N2_X;

public:

  N2();

  ~N2() {};

  N2(std::array<double, 3> CoM, quaternion Q);

};

#endif
