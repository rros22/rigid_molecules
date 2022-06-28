#include "debug.hpp"
//#include <fstream>

void csv_site_force(site_forces* forces, std::ofstream& file)
{
  file << forces->fx << ',';
  file << forces->fy << ',';
  file << forces->fz << std::endl;
}

void csv_CoM_force(lin_dyn_x* x_lin_dyn, lin_dyn_y* y_lin_dyn, lin_dyn_z* z_lin_dyn, std::ofstream& file)
{
  file << x_lin_dyn->com_Fx << ',';
  file << y_lin_dyn->com_Fy << ',';
  file << z_lin_dyn->com_Fz << std::endl;
}

void csv_forces(h2o_buffer& molecules, std::string path)
{
  std::ofstream file;
  file.open(path, std::fstream::app);
  //print CoM forces
  int molecule_no = molecules.n;
  lin_dyn_x* x_lin_dyn = molecules.x_lin_dyn;
  lin_dyn_y* y_lin_dyn = molecules.y_lin_dyn;
  lin_dyn_z* z_lin_dyn = molecules.z_lin_dyn;

  for (int i = 0; i < molecule_no; i++)
  {
    csv_CoM_force(x_lin_dyn, y_lin_dyn, z_lin_dyn, file);
  }
}
