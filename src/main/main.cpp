#include <iostream>
#include <vector>
#include <array>
#include <cmath>

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../box/box.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"

int main(){

  //Box* domain = new Box(100, 100, 100, 3*3*3);



  //box_pdb(domain, "results/domain.pdb");

  //std::cout << domain->get_molecules().size();

  //delete domain;

  //set position and orientation of molecule
  std::array<double, 3> position = {0, 0, 0};

  quaternion orientation;
  std::array<double, 3> axis = {0, 0, 1};
  double angle = 90*M_PI/180;
  orientation.set_quaternion(angle, axis);

  rigid_molecule* molecule = new water(position, orientation);


  int k = 1;

  molecule_pdb(molecule, k, "results/domain.pdb");

  delete molecule;

  return 0;
}
