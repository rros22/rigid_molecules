#include <iostream>
#include <vector>
#include <array>
#include <cmath>

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../box/box.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"

int main(){

  Box domain(100, 100, 100, 3*3*3);



  box_pdb(domain, "results/domain.pdb");


  return 0;
}
