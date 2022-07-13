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

struct water_site_positions
{
  site_positions O;
  site_positions H1;
  site_positions H2;
  site_positions q1;
};

//site forces
struct site_forces
{
  double fx;
  double fy;
  double fz;
};


struct water_site_forces
{
  site_forces O;
  site_forces H1;
  site_forces H2;
  site_forces q1;
};

//linear dynamics components
struct lin_dyn_x
{
  double com_x;
  double com_u;
  double com_Fx;
  double com_Fx_n;
};

struct lin_dyn_y
{
  double com_y;
  double com_v;
  double com_Fy;
  double com_Fy_n;
};

struct lin_dyn_z
{
  double com_z;
  double com_w;
  double com_Fz;
  double com_Fz_n;
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

  water_site_positions* water_site_pos; //array of site positions both charges and lj (in struct)
  water_site_forces* water_site_fr;     //array of site forces (in struct)
  lin_dyn_x* x_lin_dyn;     //array of x linear dynamics components
  lin_dyn_y* y_lin_dyn;     //array of y linear dynamics components
  lin_dyn_z* z_lin_dyn;     //array of z linear dynamics components
  quaternion* orientations; //array of molecule orientations

  h2o_buffer(int i);
  h2o_buffer(unsigned n);
  void allocate(unsigned n);
  //initialise particles positions
  void site_global_coordiantes();
  void initialise(double xyz[6]);
  void debug();


};


#endif
