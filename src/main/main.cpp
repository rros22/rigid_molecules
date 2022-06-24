#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"


int main(int argc, char* argv[])
{
  unsigned molecule_no = std::stoi(argv[1]);
  std::string path = "results/results.pdb";
  double domain[6] = {-50, 50, -50, 50, -50, 50};
  h2o_buffer water_molecules(molecule_no);
  water_molecules.initialise(domain);
  water_buffer_pdb(&water_molecules, path);

  return 0;
}
