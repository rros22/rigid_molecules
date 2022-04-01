#include "box.hpp"

#include <cmath>

/*

Box

*/

void Box::set_dimensions(double x, double y, double z){

  this->X = {x, y, z};
}

Box::Box(double x, double y, double z, int molecule_no){

  //set box dimensions
  set_dimensions(x, y, z);

  //add molecules to vector
  for (int i = 0; i < molecule_no; i++){

    molecules.push_back(new water());
  }

  //initialise molecule positions randomly

  int side = cbrt(molecules.size());

  //define quaternion (1, 1, 1), 45 deg
  std::array<double, 3> axis = {1, 1, 1};
  quaternion q;
  q.set_quaternion(45, axis);

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

Box::~Box(){

  //free memory from all elements of vector
  for (auto p : molecules){

    delete p;
  }

}

std::vector<rigid_molecule*> Box::get_molecules(){

  return molecules;
}
