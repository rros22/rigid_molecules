#include "rigid_molecules.hpp"
#include "../quaternions/quaternions.hpp"
#include <cmath>


/*
  Water molecule buffer
*/
h2o_buffer::h2o_buffer(int i)
{
  std::cout << "Default constructor" << std::endl;
}

h2o_buffer::h2o_buffer(unsigned n)
{
  allocate(n); //allocale memory for buffer
}

void h2o_buffer::allocate(unsigned n)
{
  //set parameters
  this->n = n;
  //test for allocation
  if (buffer != NULL)
  {
    free(buffer);
  }
  //calculate memory required
  size_t buffer_size = n*(sizeof(water_site_positions) + sizeof(water_site_forces) +
                          3*sizeof(lin_dyn_x) + sizeof(quaternion));
  //allocate
  buffer = malloc(buffer_size);
  //initialise pointers to different segments of the buffer
  water_site_pos = (water_site_positions *) buffer;
  water_site_fr = (water_site_forces *)(water_site_pos + n); //4 sites per molecule
  x_lin_dyn = (lin_dyn_x *)(water_site_fr + n);
  y_lin_dyn = (lin_dyn_y *)(x_lin_dyn + n);
  z_lin_dyn = (lin_dyn_z *)(y_lin_dyn + n);
  orientations = (quaternion *)(z_lin_dyn + n);

}

void h2o_buffer::initialise(double xyz[6])
{
  double x1 = xyz[0];
  double x2 = xyz[1];
  double y1 = xyz[2];
  double y2 = xyz[3];
  double z1 = xyz[4];
  double z2 = xyz[5];

  int side = cbrt(n);

  double angle = 90*M_PI/180;
  double axis[3] = {0, 0, 1};

  for (int i = 0; i < side; i++)
  {
    for (int j = 0; j < side; j++)
    {
      for (int k = 0; k < side; k++)
      {
      //position
      x_lin_dyn[i + side*(j + side*k)].com_x = x1 + (i + 1)*(x2-x1)/side;
      y_lin_dyn[i + side*(j + side*k)].com_y = y1 + (j + 1)*(y2-y1)/side;
      z_lin_dyn[i + side*(j + side*k)].com_z = z1 + (k + 1)*(z2-z1)/side;
      //orientation
      orientations[i + side*(j + side*k)].set_quaternion(angle, axis);
      }
    }
  }

  site_global_coordiantes();

}

void h2o_buffer::site_global_coordiantes()
{
  //site local coordinates
  double O_L[3] = {0, -0.065555, 0};
  double H1_L[3] = {-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};
  double H2_L[3] = {0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};
  double q1_L[3] = {0, 0.15 - 0.065555, 0};

  //array to store coordinates
  double global_coordinates[3];
  double offset[3];

  //loop through all molecules
  for (int i = 0; i < n; i++)
  {
    offset[0] = x_lin_dyn[i].com_x;
    offset[1] = y_lin_dyn[i].com_y;
    offset[2] = z_lin_dyn[i].com_z;

    //O
    orientations[i].transform_vector(O_L, offset, global_coordinates);
    water_site_pos[i].O.x = global_coordinates[0];
    water_site_pos[i].O.y = global_coordinates[1];
    water_site_pos[i].O.z = global_coordinates[2];
    //H1
    orientations[i].transform_vector(H1_L, offset, global_coordinates);
    water_site_pos[i].H1.x = global_coordinates[0];
    water_site_pos[i].H1.y = global_coordinates[1];
    water_site_pos[i].H1.z = global_coordinates[2];
    //H2
    orientations[i].transform_vector(H2_L, offset, global_coordinates);
    water_site_pos[i].H2.x = global_coordinates[0];
    water_site_pos[i].H2.y = global_coordinates[1];
    water_site_pos[i].H2.z = global_coordinates[2];
    //q1
    orientations[i].transform_vector(q1_L, offset, global_coordinates);
    water_site_pos[i].q1.x = global_coordinates[0];
    water_site_pos[i].q1.y = global_coordinates[1];
    water_site_pos[i].q1.z = global_coordinates[2];
  }
}

//debug
void h2o_buffer::debug()
{
  allocate(2);

  double x1 = -5;
  double x2 = 5;
  double y1 = -5;
  double y2 = 5;
  double z1 = -5;
  double z2 = 5;

  double angle = -90*M_PI/180;
  double angle_2 = 0;
  double axis[3] = {0.7, 0, 1};
  double axis_2[3] = {0, 3.4, 1};

  x_lin_dyn[0].com_x = 3;
  y_lin_dyn[0].com_y = -2;
  z_lin_dyn[0].com_z = 0.5;

  x_lin_dyn[1].com_x = 4;
  y_lin_dyn[1].com_y = 3.2;
  z_lin_dyn[1].com_z = -2.6;

  orientations[0].set_quaternion(angle, axis);
  orientations[1].set_quaternion(angle_2, axis_2);


  site_global_coordiantes();

}
