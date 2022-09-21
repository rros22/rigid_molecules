#include "debug.hpp"
#include <cmath>
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
    csv_CoM_force(&x_lin_dyn[i], &y_lin_dyn[i], &z_lin_dyn[i], file);
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

void csv_orientations(h2o_buffer& molecules, std::string path, double t)
{
  std::ofstream file;
  file.open(path, std::fstream::app);
  //print CoM forces
  unsigned molecule_no = molecules.n;
  quaternion_orientation* orientations = molecules.orientations;

  file << "Time: " << t << std::endl;
  for (int i = 0; i < molecule_no; i++)
  {
    file << "MOLECULE " << i + 1 << ": " << " Angle: " << orientations[i].current.get_angle()*180/M_PI;
    file << " | Axis: (" << orientations[i].current.get_axis(0) << ", "
         << orientations[i].current.get_axis(1) << ", "
         << orientations[i].current.get_axis(2) << ") " << std::endl;
    file << "MOLECULE " << i + 1 << ": " << " Angle: " << orientations[i].next.get_angle()*180/M_PI;
    file << " | Axis: (" << orientations[i].next.get_axis(0) << ", "
         << orientations[i].next.get_axis(1) << ", "
         << orientations[i].next.get_axis(2) << ") " << std::endl;
  }
    file << std::endl;

}
