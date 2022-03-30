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

  I[0] = 0.614812;
  I[1] = 1.154932;
  I[2] = 1.78001;

  O1.set_parameters(3.15365, 78*k_b, 15.999);
  H1.set_parameters(0.52, 1.00784);
  H2.set_parameters(0.52, 1.00784);

  q1.set_parameters(-1.04, 0);

  //local site coordinates
  O1_X = {0, 0.65555, 0};
  H1_X = {-0.9572*sin(52.26), 0.9572*cos(52.26) - 0.65555, 0};
  H2_X = {0.9572*sin(52.26), 0.9572*cos(52.26) - 0.65555, 0};

  q1_X = {0, 0.15, 0};

  //set site coordinates
  O1.set_local_coordinates(O1_X);
  H1.set_local_coordinates(H1_X);
  H2.set_local_coordinates(H2_X);
  q1.set_local_coordinates(q1_X);

}
