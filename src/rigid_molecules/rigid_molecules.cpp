#include "rigid_molecules.hpp"
#include <cmath>


/*
  Water molecule buffer
*/

h2o_buffer::h2o_buffer(unsigned n)
{
  allocate(n); //allocale memory for buffer
}

void h2o_buffer::allocate(unsigned n)
{
  if (buffer != NULL)
  {
    free(buffer);
  }
  else
  {
    //calculate memory required
    size_t buffer_size = n*(4*sizeof(site_positions) + sizeof(site_forces) +
                            3*sizeof(lin_dyn_x) + sizeof(quaternion));
    //allocate
    buffer = malloc(buffer_size);
    //initialise pointers to different segments of the buffer
    site_pos = (site_positions *) buffer;
    site_fr = (site_forces *)(site_pos + 4*n); //4 sites per molecule
    x_lin_dyn = (lin_dyn_x *)(site_fr + n);
    y_lin_dyn = (lin_dyn_y *)(x_lin_dyn + n);
    z_lin_dyn = (lin_dyn_z *)(y_lin_dyn + n);
    orientations = (quaternion *)(z_lin_dyn + n);
  }
}



/*set site coordinates
  sites[0]->set_local_coordinates({0, -0.065555, 0});
  sites[1]->set_local_coordinates({-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  sites[2]->set_local_coordinates({0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0});
  sites[3]->set_local_coordinates({0, 0.15 - 0.065555, 0});
*/
