#include <iostream>
#include <vector>

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../box/box.hpp"
#include "../post_processing/post_processing.hpp"

int main(){

  Box* domain = new Box(100, 100, 100, 3*3*3);



  box_pdb(domain, "results/domain.pdb");

  std::cout << domain->get_molecules().size();

  delete domain;
  
  return 0;
}
