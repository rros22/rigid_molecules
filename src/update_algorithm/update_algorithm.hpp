#ifndef UPDATE_ALGORITHM_HPP
#define UPDATE_ALGORITHM_HPP

#include <iostream>
#include "../rigid_molecules/rigid_molecules.hpp"

/*
    1. Calculate site positions from molecule positions and orientations
    2. Calculate site forces
    3. Position integration
        3.a. Calculate CoM forces
        3,b. Calculate velocity
        3.b. Verlet integrate position
    4. Orientation integration
        4.a. Calculate torques
        4.b. Calculate quaternion derivatives
        4.c. Verlet integrate
*/

/*
  1. SIte position calculation
*/



/*
  1. Site force calculation
*/

//site forces
void coulombic_force(site_positions* site_a, site_positions* site_b, site_forces* forces, double q_a, double q_b);
void lj_force(site_positions* site_a, site_positions* site_b, site_forces* forces, double sigma, double epsilon);
void accumulate_force_sites(site_forces* site_a, site_forces* site_b, site_forces* forces);
//set site forces molecules
void reset_forces_site(site_forces* site);
void set_forces_sites(h2o_buffer* water_molecules);
//calculate CoM forces
void set_CoM_force(h2o_buffer* water_molecules);
void set_CoM_force_n(h2o_buffer* water_molecules);

/*
  3. Position integration
*/

void next_position(h2o_buffer* water_molecules, double dt);
void next_velocity(h2o_buffer* water_molecules, double dt);
void verlet_integrate(h2o_buffer* water_molecules, double dt);

/*
  4. Orientation integration
*/

void compute_torques(h2o_buffer* water_molecules);
void compute_angular_velocity(quaternion orientation, quaternion velocity, double omega[3]);
void compute_quat_accel(h2o_buffer* water_molecules);
void half_quat_velocity(h2o_buffer* water_molecules, double dt);
void next_orientation(h2o_buffer* water_molecules, double dt);

void cross_product(double vector_a[3], double vector_b[3], double output[3]);
double torque_magnitude(torques* torque);

#endif
