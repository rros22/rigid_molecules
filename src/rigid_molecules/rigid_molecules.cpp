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

void rigid_molecule::set_global_coordinates(){

  //rotate local coordinates to molecule orientation
  rot_matrix rotation(this->Q);

  //store global coordinates, and intermediate states of calculation
  std::array<double, 3> global;

  //loop through all sites in molecules
  for (int i = 0; i < sites.size(); i++){

    //tranform local coordinates according to molecule orientation
    global = rotation.matrix_vector(sites[i]->get_local_coordinates());

    //offset by center of mass position
    for (int j = 0; j < 3; j++){

      global[j] += CoM[j];
    }

    //set global coordinate of site
    sites[i]->set_global_coordinates(global);

  }
}

std::vector<site*> rigid_molecule::return_sites_list(){

  return sites;
}

std::array<double, 3> rigid_molecule::return_coordinates_site(int i){

  if (i >= 0 && i < sites.size()){

    return sites[i]->get_global_coordinates();

  }

  else{

    std::cout << "Index i is out of bounds: cannot return site." << std::endl;

    return {0, 0, 0};
  }

}

/*

H2O

*/

H2O::H2O(){

  //define interaction site pointer array
  sites.push_back(&O_1);
  sites.push_back(&H_1);
  sites.push_back(&H_2);
  sites.push_back(&q_1);

  //set site names
  O_1.set_symbol("O");
  H_1.set_symbol("H");
  H_2.set_symbol("H");
  q_1.set_symbol("X");

  //set interaction site parameters
  O_1.set_parameters(3.15365, 78*k_b, 15.999);
  H_1.set_parameters(0.52, 1.00784);
  H_2.set_parameters(0.52, 1.00784);

  q_1.set_parameters(-1.04, 0);

  //set site coordinates
  O_1.set_local_coordinates({0, -0.065555, 0});
  H_1.set_local_coordinates({-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  H_2.set_local_coordinates({0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  q_1.set_local_coordinates({0, 0.15 - 0.065555, 0});

  //define molecular mass
  set_mass();

  //define molecular moments of inertia
  I[0] = 0.614812;
  I[1] = 1.154932;
  I[2] = 1.78001;



}

H2O::H2O(std::array<double, 3> CoM, quaternion Q): H2O(){

    set_position(CoM);
    set_orientation(Q);

    set_global_coordinates();
}

/*

N2

*/

N2::N2(){

  //define interaction site pointer array
  sites.push_back(&N_1);
  sites.push_back(&N_2);

  //for debugging purposes
  sites.push_back(&q_1);

  //set site names
  N_1.set_symbol("N");
  N_2.set_symbol("N");

  //set interaction site parameters
  N_1.set_parameters(3.17, 78*k_b, 14.0067);
  N_2.set_parameters(3.17, 78*k_b, 14.0067);

  //set site coordinates
  N_1.set_local_coordinates({0, -0.545, 0});
  N_2.set_local_coordinates({0, 0.545, 0});

  //define molecular mass
  set_mass();

  //define molecular moments of inertia
  I[0] = 8.32068;
  I[1] = 0;
  I[2] = 8.32068;
}
