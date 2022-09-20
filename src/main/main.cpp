#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"
#include "../update_algorithm/update_algorithm.hpp"
#include "../parse/parse.hpp"
#include "../debug/debug.hpp"


int main(int argc, char* argv[])
{
  unsigned molecule_no;
  if (parse_molecules(argc, argv, molecule_no)){return 1;}

  std::string path = "results/results.pdb";
  std::string debug_path = "debug/results.csv";

  //allocate memory
  h2o_buffer water_molecules;
  water_molecules.debug();

  //write molecular coordinates
  water_buffer_pdb(&water_molecules, path);
  //write molecular and site forces
  csv_forces(water_molecules, debug_path);
  //array to store torques
  double torque[3];

  //define number of iterations
  int i = 0;
  int iterations = 10000;
  while(i < iterations)
  {
    verlet_integrate(&water_molecules, 0.001);
    quaternion_verlet_integrate(&water_molecules, 0.001);
    water_buffer_pdb(&water_molecules, path);
    i++;
  }

  return 0;
}
