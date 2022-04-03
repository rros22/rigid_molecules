#include "rigid_molecules.hpp"
#include <cmath>

#define k_b 1.38064852E-23 * 1E20 * 1/1.66054 * 1E27

/*

Rigid molecule

*/

void rigid_molecule::set_position(std::array<double, 3> &CoM){

  this->CoM = CoM;
}

void rigid_molecule::set_orientation(quaternion &Q){

  this->Q = Q;
}

void rigid_molecule::set_mass(){

  //reset mass
  this->m = 0;

  //loop through site list
  for (int i = 0; i < sites.size(); i++){

    this->m += sites[i]->get_mass();
  }
}

double rigid_molecule::get_mass(){

  return this->m;
}

/*

Water

*/

void water::set_global_coordinates(){

  //oxygen atom
  //rotate local coordinates to molecule orientation
  rot_matrix rotation(this->Q);
  std::array<double, 3> global = rotation.matrix_vector(O1.get_local_coordinates());

  //add center of mass position
  for (int i = 0; i < 3; i++){

    global[i] = global[i] + CoM[i];
  }

  O1.set_global_coordinates(global);

  //hydrogen atom 1
  //rotate local coordinates to molecule orientation
  global = rotation.matrix_vector(H1.get_local_coordinates());

  //add center of mass position
  for (int i = 0; i < 3; i++){

    global[i] = global[i] + CoM[i];
  }

  H1.set_global_coordinates(global);

  //hydrogen atom 2
  //rotate local coordinates to molecule orientation
  global = rotation.matrix_vector(H2.get_local_coordinates());

  //add center of mass position
  for (int i = 0; i < 3; i++){

    global[i] = global[i] + CoM[i];
  }

  H2.set_global_coordinates(global);

  //dummy charge
  //rotate local coordinates to molecule orientation
  global = rotation.matrix_vector(q1.get_local_coordinates());

  //add center of mass position
  for (int i = 0; i < 3; i++){

    global[i] = global[i] + CoM[i];
  }

  q1.set_global_coordinates(global);


}

water::water(){

  //define interaction site pointer array
  sites.push_back(&O1);
  sites.push_back(&H1);
  sites.push_back(&H2);
  sites.push_back(&q1);

  //set site names
  O1.set_symbol("O");
  H1.set_symbol("H");
  H2.set_symbol("H");
  q1.set_symbol("X");

  //set interaction site parameters
  O1.set_parameters(3.15365, 78*k_b, 15.999);
  H1.set_parameters(0.52, 1.00784);
  H2.set_parameters(0.52, 1.00784);

  q1.set_parameters(-1.04, 0);

  //define site coordinates

  O1_X = {0, -0.065555, 0};
  H1_X = {-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};
  H2_X = {0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};

  q1_X = {0, 0.15 - 0.065555, 0};

  //set site coordinates
  O1.set_local_coordinates(O1_X);
  H1.set_local_coordinates(H1_X);
  H2.set_local_coordinates(H2_X);
  q1.set_local_coordinates(q1_X);

  //define molecular mass
  set_mass();

  //define molecular moments of inertia
  I[0] = 0.614812;
  I[1] = 1.154932;
  I[2] = 1.78001;



}

water::water(std::array<double, 3> CoM, quaternion Q): water(){

    set_position(CoM);
    set_orientation(Q);

    set_global_coordinates();
}

std::array<double, 3> water::return_coordinates_site(int i){

  return sites[i]->get_global_coordinates();

}

std::vector<site*> water::return_sites_list(){

  return sites;
}
