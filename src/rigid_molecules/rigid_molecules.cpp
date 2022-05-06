#include "rigid_molecules.hpp"
#include <cmath>

#define k_b 1.38064852E-23 * 1E20 * 1/1.66054 * 1E27

/*

Rigid molecule

*/

//setter functions
void rigid_molecule::set_position(std::array<double, 3> &CoM)
{
  this->CoM = CoM;
}

void rigid_molecule::set_orientation(quaternion &Q)
{
  this->Q = Q;
}

void rigid_molecule::set_mass(double m)
{
  //reset mass
  this->m = m;
}

void rigid_molecule::set_mom_inertia(std::array<double, 3> I)
{
  this->I = I;
}

void rigid_molecule::set_symbol(std::string symbol)
{
  this->symbol = symbol;
}

//getter functions
double rigid_molecule::get_mass()
{
  return this->m;
}

const std::string& rigid_molecule::get_symbol()
{
  return symbol;
}


void rigid_molecule::set_CoM_force(std::array<double, 3> F)
{
    this->F[0] += F[0];
    this->F[1] += F[1];
    this->F[2] += F[2];
}

void rigid_molecule::reset_CoM_force()
{
  F[0] = 0;
  F[1] = 0;
  F[2] = 0;
}

/*

H2O

*/

H2O::H2O(){

  //set molecule symbol
  kernel.set_symbol("H2O");

  //set site names
  O1.kernel.set_symbol("O");
  H1.kernel.set_symbol("H");
  H2.kernel.set_symbol("H");
  q.kernel.set_symbol("X");

  //set interaction site parameters
  O1.set_parameters(3.15365, 78*k_b, 15.999);
  H1.set_parameters(0.52, 1.00784);
  H2.set_parameters(0.52, 1.00784);

  q.set_parameters(-1.04, 0);

  //set site coordinates
  O1.kernel.set_local_coordinates({0, -0.065555, 0});
  H1.kernel.set_local_coordinates({-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  H2.kernel.set_local_coordinates({0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  q.kernel.set_local_coordinates({0, 0.15 - 0.065555, 0});

  //define molecular mass
  kernel.set_mass(O1.kernel.get_mass() + H1.kernel.get_mass() + H2.kernel.get_mass());

  //define molecular moments of inertia
  std::array<double, 3> I = {0.614812, 1.154932, 1.78001};

  kernel.set_mom_inertia(I);

}

H2O::H2O(std::array<double, 3> CoM, quaternion Q): H2O(){

    kernel.set_position(CoM);
    kernel.set_orientation(Q);

}

void H2O::set_global_coordinates(){

  //rotate local coordinates to molecule orientation
  rot_matrix rotation(kernel.get_orientation());

  //store global coordinates, and intermediate states of calculation
  std::array<double, 3> global;
  std::array<double, 3> CoM;

  //O1
  global = rotation * O1.kernel.get_local_coordinates();
  CoM = kernel.get_position();

  for (int j = 0; j < 3; j++)
  {
    global[j] += CoM[j];
  }

  O1.kernel.set_global_coordinates(global);

  //H1
  global = rotation * H1.kernel.get_local_coordinates();
  CoM = kernel.get_position();

  for (int j = 0; j < 3; j++)
  {
    global[j] += CoM[j];
  }

  H1.kernel.set_global_coordinates(global);

  //H2
  global = rotation * H2.kernel.get_local_coordinates();
  CoM = kernel.get_position();

  for (int j = 0; j < 3; j++)
  {
    global[j] +=CoM[j];
  }

  H2.kernel.set_global_coordinates(global);

  //q
  global = rotation * q.kernel.get_local_coordinates();
  CoM = kernel.get_position();

  for (int j = 0; j < 3; j++)
  {
    global[j] += CoM[j];
  }

  q.kernel.set_global_coordinates(global);

}

/*
void H2O::set_forces(rigid_molecule *molecule){

  //get list of interaction sites
  std::vector<site*> sites = molecule->return_sites_list();

  //reset forces
  for (int i = 0; i < sites.size(); i++){

    sites[i]->reset_forces();
  }

  //select between different child classes
  if (molecule->get_symbol() == "H2O"){

    //H2O *h2o = (H2O *) molecule;

    //forces on oxygen atom
    sites[0]->calculate_forces((lj_site *) sites[0]);

    //forces on hydrogen atom 1
    sites[1]->calculate_forces((charge *) sites[1]);
    sites[1]->calculate_forces((charge *) sites[2]);
    sites[1]->calculate_forces((charge *) sites[3]);

    //forces on hydrogen atom 2
    sites[2]->calculate_forces((charge *) sites[1]);
    sites[2]->calculate_forces((charge *) sites[2]);
    sites[2]->calculate_forces((charge *) sites[3]);

    //forces on dummy charge
    sites[3]->calculate_forces((charge *) sites[1]);
    sites[3]->calculate_forces((charge *) sites[2]);
    sites[3]->calculate_forces((charge *) sites[3]);

    //calculate force on CoM
    set_CoM_force();

  }

  else if (molecule->get_symbol() == "N2"){

    //forces on oxygen atom
    sites[0]->calculate_forces((lj_site *) sites[0]);
    sites[0]->calculate_forces((lj_site *) sites[1]);

    //calculate force on CoM
    set_CoM_force();

  }

}
*/
