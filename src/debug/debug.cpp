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
  file << "Center of Mass forces" << std::endl << std::endl;
  int molecule_no = molecules.n;
  lin_dyn_x* x_lin_dyn = molecules.x_lin_dyn;
  lin_dyn_y* y_lin_dyn = molecules.y_lin_dyn;
  lin_dyn_z* z_lin_dyn = molecules.z_lin_dyn;
  water_site_forces* water_site_fr = molecules.water_site_fr;

  for (int i = 0; i < molecule_no; i++)
  {
    file << i + 1 << ". " ;
    csv_CoM_force(x_lin_dyn, y_lin_dyn, z_lin_dyn, file);
  }
  file << std::endl;
  //print site forces
  file << "Site forces" << std::endl << std::endl;
  for (int i = 0; i < molecule_no; i++)
  {
    file << i + 1 << "_" << "O. " ;
    csv_site_force(&water_site_fr[i].O, file);
    file << i + 1 << "_" << "H1. " ;
    csv_site_force(&water_site_fr[i].H1, file);
    file << i + 1 << "_" << "H2. " ;
    csv_site_force(&water_site_fr[i].H2, file);
    file << i + 1 << "_" << "q1. " ;
    csv_site_force(&water_site_fr[i].q1, file);
  }
  file << std::endl;

}
