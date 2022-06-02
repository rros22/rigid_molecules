#ifndef BOX_HPP
#define BOX_HPP

#include <array>
#include <vector>
#include <memory>

#include "../rigid_molecules/rigid_molecules.hpp"

class Box{

private:
  //dimensions
  std::array<double, 3> X;
  //molecules (pointers, due to being an abstract class)
  std::vector<std::shared_ptr<rigid_molecule>> molecules;

  void set_dimensions(double x, double y, double z);

public:
  //create box with specified dimensions
  Box(double x, double y, double z, int molecule_no);

  //getter functions
  const std::array<double, 3>& get_dimensions();
  std::vector<std::shared_ptr<rigid_molecule>>& get_molecules();
  std::shared_ptr<rigid_molecule> get_molecule(int i);
//test upload
};

#endif
