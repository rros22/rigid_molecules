#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"
#include "../update_algorithm/update_algorithm.hpp"
#include "../debug/debug.hpp"


int main(int argc, char* argv[])
{
  if (argv[1] == NULL)
  {
    std::cout << "Program executed without providing enough arguments"
    << std::endl << "Number of molecules unspecified" << std::endl;
    return 1;
  }

  std::string path = "results/results.pdb";
  std::string debug_path = "debug/results.csv";

  unsigned molecule_no = std::stoi(argv[1]);

  h2o_buffer water_molecules(1);
  water_molecules.debug();
  water_buffer_pdb(&water_molecules, path);
  csv_forces(water_molecules, debug_path);

  int i = 0;
  int iterations = 1;
  while(i < iterations)

  {
    verlet_integrate(&water_molecules, 1E-15);
    csv_forces(water_molecules, debug_path);
    //water_buffer_pdb(&water_molecules, path);
    i++;
    std::cout << i << std::endl;
  }

  return 0;
}
