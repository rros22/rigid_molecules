#ifndef UPDATE_ALGORITHM_HPP
#define UPDATE_ALGORITHM_HPP

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

//site forces
void coulombic_force(site_positions* site_a, site_positions* site_b, double q_a, double q_b, site_forces* forces);
void lj_force(site_positions* site_a, site_positions* site_b, double sigma, double epsilon, site_forces* forces);

void set_forces();

#endif
