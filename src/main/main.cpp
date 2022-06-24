#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"
#include "../update_algorithm/update_algorithm.hpp"


int main(int argc, char* argv[])
{
  if (argv[1] == NULL) {return 1;}

  unsigned molecule_no = std::stoi(argv[1]);
  std::string path = "results/results.pdb";
  double domain[6] = {-5, 5, -5, 5, -5, 5};
  h2o_buffer water_molecules(molecule_no);
  water_molecules.initialise(domain);
  water_buffer_pdb(&water_molecules, path);
  int i = 0;
  while(1)
  /*
  {
    verlet_integrate(&water_molecules, 1E-15);
    water_buffer_pdb(&water_molecules, path);
    std::cout << i++ << std::endl;
    if (i > 1000) {break;}
  }
  */
  return 0;
}
