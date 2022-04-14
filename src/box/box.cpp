#include "box.hpp"

#include <cmath>

/*

Box

*/

void Box::set_dimensions(double x, double y, double z){

  this->X = {x, y, z};
}

Box::Box(double x, double y, double z, int water_no, int nitrogen2_no){

  //set box dimensions
  set_dimensions(x, y, z);

  //add molecules
  water_molecules.assign(water_no, H2O());
  nitrogen2_molecules.assign(nitrogen2_no, N2());

  //allocate vector of molecules
  for (int i = 0; i < water_no; i++){

    molecules.push_back(&(water_molecules[i]));
  }

  for (int i = 0; i < nitrogen2_no; i++){

    molecules.push_back(&(nitrogen2_molecules[i]));
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
  for (int i = 0; i < molecules.size(); i++){

    molecules[i]->set_global_coordinates();
  }


}

Box::~Box(){

}

std::vector<rigid_molecule*> Box::get_molecules(){

  return molecules;
}
