#ifndef BOX_HPP
#define BOX_HPP

#include <array>
#include <vector>

#include "../rigid_molecules/rigid_molecules.hpp"

class Box{

private:
  //dimensions
  std::array<double, 3> X;

  //molecules
  std::vector<rigid_molecule> molecules;

  void set_dimensions(double x, double y, double z);

public:

  Box(double x, double y, double z, int molecule_no);

  rigid_molecule* get_molecule(int i);

  std::vector<rigid_molecule> get_molecules();

};

#endif
