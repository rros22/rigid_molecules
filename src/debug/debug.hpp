#ifndef DEBUG_HPP
#define DEBUG_HPP

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"
#include "../update_algorithm/update_algorithm.hpp"

void csv_site_force(site_forces* forces, std::ofstream& file);
void csv_CoM_force(lin_dyn_x* x_lin_dyn, lin_dyn_y* y_lin_dyn, lin_dyn_z* z_lin_dyn, std::ofstream& file);
void csv_forces(h2o_buffer& molecules, std::string path);
void csv_orientations(h2o_buffer& molecules, std::string path, double t);

#endif
