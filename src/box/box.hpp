#ifndef BOX_HPP
#define BOX_HPP

#include <array>
#include <vector>

#include "../rigid_molecules/rigid_molecules.hpp"

class Box{

private:
  //dimensions
  std::array<double, 3> X;

  //molecules (pointers, due to being an abstract class)
  std::vector<rigid_molecule*> molecules;

  void set_dimensions(double x, double y, double z);

public:

  Box(double x, double y, double z, int molecule_no);

  //explicit destructor, object contains rerferences to heap allocated members
  ~Box();

  rigid_molecule* get_molecule(int i);

  std::vector<rigid_molecule*> get_molecules();

};

#endif
