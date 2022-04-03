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

  //delete domain;

  rigid_molecule* molecule = new water();

  std::cout << molecule->get_mass() << std::endl;

  return 0;
}
