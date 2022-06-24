#ifndef DEBUG_HPP
#define DEBUG_HPP

#include "../rigid_molecules/rigid_molecules.hpp"
#include "../post_processing/post_processing.hpp"
#include "../quaternions/quaternions.hpp"
#include "../update_algorithm/update_algorithm.hpp"

void csv_site_force(site_forces* forces);
void csv_CoM_force(water_site_forces* forces);
void csv_forces(h2o_buffer& molecules);

#endif
