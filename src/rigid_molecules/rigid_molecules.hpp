#ifndef RIGID_MOLECULES_HPP
#define RIGID_MOLECULES_HPP

#include <iostream>
#include "../quaternions/quaternions.hpp"

/*
  Molecular data is organised to optimise sequential access pattern
  Each molecule memory buffer contains different components, which are
  the data inputs to algorithm subroutines.
*/

/*
  molecules_buffer[N] = [ Comp A[x*N] | Comp B[y*N] | ... | Comp X[l*N] ]
*/

struct site_positions
{
  double x;
  double y;
  double z;
};

//linear dynamics components
struct lin_dyn_x
{
  double com_x;
  double com_u;
  double com_Fx;
};

struct lin_dyn_y
{
  double com_y;
  double com_v;
  double com_Fy;
};

struct lin_dyn_z
{
  double com_z;
  double com_w;
  double com_Fz;
};

//memory buffers for each molecule type

/*
  Water molecule:

  [ site_positions[n] | site_forces_x[n] | site_forces_y[n] | site_forces_z[n] |
    lin_dyn_x[n] | lin_dyn_y[n] | lin_dyn_z[n] ]
*/
struct h2o_buffer
{
  unsigned n;   //number of water molecules
  void* buffer = NULL; // memory for the water molecules

  site_positions* site_pos; //array of site positions both charges and lj
  double* site_forces_x;          //array of x site forces
  double* site_forces_y;          //array of y site forces
  double* site_forces_z;          //array of z site forces
  lin_dyn_x* x_lin_dyn;           //array of x linear dynamics components
  lin_dyn_y* y_lin_dyn;           //array of y linear dynamics components
  lin_dyn_z* z_lin_dyn;           //array of z linear dynamics components
  quaternion* orientations;       //array of molecule orientations

  h2o_buffer(){};
  h2o_buffer(unsigned n);
  void allocate(unsigned n);

};


#endif
