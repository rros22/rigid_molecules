#include "rigid_molecules.hpp"
#include <cmath>

#define k_b 1.38064852E-23 * 1E20 * 1/1.66054 * 1E27

/*

Rigid molecule

*/
void rigid_molecule::set_symbol(std::string symbol){

  this->symbol = symbol;
}

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

std::string rigid_molecule::get_symbol(){

  return symbol;
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

std::vector<std::shared_ptr<site>> rigid_molecule::return_sites_list(){

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

void rigid_molecule::reset_CoM_force(){

  F[0] = 0;
  F[1] = 0;
  F[2] = 0;
}

void rigid_molecule::set_CoM_force(){

  //reset previous CoM force
  reset_CoM_force();

  //temp array to store contributions of sites
  std::array<double, 3> forces;

  //add contributions from every site
  for (int i = 0; i < sites.size(); i++){

    //obtain forces from each site
    forces = sites[i]->get_forces();

    //add (Fx, Fy, Fz) contributions from every site (vector addition)
    F[0] += forces[0];
    F[1] += forces[1];
    F[2] += forces[2];

  }

}

/*

H2O

*/

H2O::H2O(){

  //set molecule symbol
  set_symbol("H2O");

  //define interaction sites
  sites.emplace_back(new lj_site(3.15365, 78*k_b, 15.999));
  sites.emplace_back(new charge(0.52, 1.00784));
  sites.emplace_back(new charge(0.52, 1.00784));
  sites.emplace_back(new charge(-1.04, 0));

  //set site names
  sites[0]->set_symbol("O");
  sites[1]->set_symbol("H");
  sites[2]->set_symbol("H");
  sites[3]->set_symbol("X");

  //set interaction site parameters
  //sites[0]->set_parameters(3.15365, 78*k_b, 15.999);
  //sites[1]->set_parameters(0.52, 1.00784);
  //sites[2]->set_parameters(0.52, 1.00784);

  //sites[3]->set_parameters(-1.04, 0);

  //set site coordinates
  sites[0]->set_local_coordinates({0, -0.065555, 0});
  sites[1]->set_local_coordinates({-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  sites[2]->set_local_coordinates({0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  sites[3]->set_local_coordinates({0, 0.15 - 0.065555, 0});

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

/*

N2

*/

N2::N2(){

  //set molecule symbol
  set_symbol("N2");

  //define interaction site pointer array
  sites.emplace_back(new lj_site(3.17, 78*k_b, 14.0067));
  sites.emplace_back(new lj_site(3.17, 78*k_b, 14.0067));

  //for debugging purposes
  sites.emplace_back(new charge(0, 0));

  //set site names
  sites[0]->set_symbol("N");
  sites[1]->set_symbol("N");

  //set interaction site parameters
  //sites[0]->set_parameters(3.17, 78*k_b, 14.0067);
  //sites[1]->set_parameters(3.17, 78*k_b, 14.0067);

  //set site coordinates
  sites[0]->set_local_coordinates({0, -0.545, 0});
  sites[1]->set_local_coordinates({0, 0.545, 0});

  //define molecular mass
  set_mass();

  //define molecular moments of inertia
  I[0] = 8.32068;
  I[1] = 0;
  I[2] = 8.32068;
}

N2::N2(std::array<double, 3> CoM, quaternion Q): N2(){

    set_position(CoM);
    set_orientation(Q);

    set_global_coordinates();
}
