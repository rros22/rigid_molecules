#include "box.hpp"
#include <cmath>

/*

Box

*/

void Box::set_dimensions(double x, double y, double z)
{
  this->X = {x, y, z};
}

Box::Box(double x, double y, double z, int molecule_no)
{
  //set box dimensions
  set_dimensions(x, y, z);
  //add molecules to vector
  for (int i = 0; i < molecule_no/2; i++)
  {
    molecules.push_back(std::make_shared<H2O> ());
  }

  for (int i = molecule_no/2; i < molecule_no; i++){

    molecules.push_back(std::make_shared<N2> ());
  }

  //initialise molecule positions randomly

  int side = cbrt(molecules.size());

  //define quaternion (0, 0, 1), 90 deg, converted to radians
  std::array<double, 3> axis = {1, 0, 0};
  quaternion q;
  q.set_quaternion(90*M_PI/180, axis);

  //define position array
  std::array<double, 3> position;

  for (int i = 0; i < side; i++){

    for (int j = 0; j < side; j++){

      for (int k = 0; k < side; k++){

        //position
        position = {x/side + i*x/side, y/side + j*y/side, z/side + k*z/side};
        molecules[i + side*(j + side*k)]->set_position(position);

        //orientation
        molecules[i + side*(j + side*k)]->set_orientation(q);
      }

    }

  }

  //calculate site positions from molecule positions and orientations
  for (int i = 0; i < molecule_no; i++){

    molecules[i]->set_global_coordinates();
  }


}

const std::array<double, 3>& Box::get_dimensions()
{
  return X;
}

std::vector<std::shared_ptr<rigid_molecule>>& Box::get_molecules()
{
  return molecules;
}

std::shared_ptr<rigid_molecule> Box::get_molecule(int i)
{
  return molecules[i];
}
