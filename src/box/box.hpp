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
  std::vector<H2O> water_molecules;
  std::vector<N2> nitrogen2_molecules;

  //list of pointers to molecules
  std::vector<rigid_molecule*> molecules;

  void set_dimensions(double x, double y, double z);

public:

  Box(double x, double y, double z, int water_no, int nitrogen2_no);

  //explicit destructor, object contains rerferences to heap allocated members
  ~Box();

  rigid_molecule* get_molecule(int i);

  std::vector<rigid_molecule*> get_molecules();

};

#endif
