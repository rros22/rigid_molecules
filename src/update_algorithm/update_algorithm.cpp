#include "update_algorithm.hpp"
#include <cmath>

//#define k_b 1.38064852E-23
//#define u 1.66054E-27

/*
O mass = 15.9994
H mass = 1.008
O charge = -1.1794
H charge = 0.5897
r_0 of OH bond = 0.9572
Theta of HOH angle = 104.52
OM distance = 0.1577
LJ epsilon of O-O = 0.21084
LJ sigma of O-O = 3.1668
LJ sigma , epsilon  of OH, HH = 0.0
Coulomb cutoff = 8.5
*/

void coulombic_force(site_positions* site_a, site_positions* site_b, site_forces* forces, double q_a, double q_b)
{
  //force constant
  double k_e = 1; //8.988E9
  //local variables for simplification, compiler will inline
  double x_a = site_a->x;
  double y_a = site_a->y;
  double z_a = site_a->z;

  double x_b = site_b->x;
  double y_b = site_b->y;
  double z_b = site_b->z;

  //forces
  forces->fx = -k_e*q_a*q_b*(x_b - x_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 3/2);
  forces->fy = -k_e*q_a*q_b*(y_b - y_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 3/2);
  forces->fz = -k_e*q_a*q_b*(z_b - z_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 3/2);

}
void lj_force(site_positions* site_a, site_positions* site_b, site_forces* forces, double sigma, double epsilon)
{
  //local variables for simplification, compiler will inline
  double x_a = site_a->x;
  double y_a = site_a->y;
  double z_a = site_a->z;

  double x_b = site_b->x;
  double y_b = site_b->y;
  double z_b = site_b->z;

  //forces
  forces->fx = -4*epsilon*((12*pow(sigma, 12)*(x_b - x_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 7)) - (6*pow(sigma, 6)*(x_b - x_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 4)));
  forces->fy = -4*epsilon*((12*pow(sigma, 12)*(y_b - y_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 7)) - (6*pow(sigma, 6)*(y_b - y_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 4)));
  forces->fz = -4*epsilon*((12*pow(sigma, 12)*(z_b - z_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 7)) - (6*pow(sigma, 6)*(z_b - z_a)/pow(pow(x_b - x_a, 2) + pow(y_b - y_a, 2) + pow(z_b - z_a, 2), 4)));
}

void accumulate_force_sites(site_forces* site_a, site_forces* site_b, site_forces* forces)
{
  site_a->fx += forces->fx;
  site_a->fy += forces->fy;
  site_a->fz += forces->fz;

  site_b->fx -= forces->fx;
  site_b->fy -= forces->fy;
  site_b->fz -= forces->fz;

}

void set_forces_sites(h2o_buffer* water_molecules)
{
  //physics constants
  double sigma = 3.1589; //angstroms
  double epsilon = 0.1852;

  double q_H = 1;//0.52;
  double q_q = -1.04;
  //inlined by compiler
  unsigned molecule_no = water_molecules->n;
  water_site_positions* water_site_pos = water_molecules->water_site_pos;
  water_site_forces* water_site_fr = water_molecules->water_site_fr;
  //struct to store temp forces results
  site_forces forces;
  //call site forces functions
  for (int i = 0; i < molecule_no; i++)
  {
    for (int j = i + 1; j < molecule_no; j++)
    {
      // O - O interaction
      lj_force(&water_site_pos[i].O, &water_site_pos[j].O, &forces, sigma, epsilon);
      accumulate_force_sites(&water_site_fr[i].O, &water_site_fr[j].O, &forces);

      // H1 - H1 interaction
      coulombic_force(&water_site_pos[i].H1, &water_site_pos[j].H1, &forces, q_H, q_H);
      accumulate_force_sites(&water_site_fr[i].H1, &water_site_fr[j].H1, &forces);
      // H1 - H2 interaction
      coulombic_force(&water_site_pos[i].H1, &water_site_pos[j].H2, &forces, q_H, q_H);
      accumulate_force_sites(&water_site_fr[i].H1, &water_site_fr[j].H2, &forces);
      // H1 - q1 interaction
      coulombic_force(&water_site_pos[i].H1, &water_site_pos[j].q1, &forces, q_H, q_q);
      accumulate_force_sites(&water_site_fr[i].H1, &water_site_fr[j].q1, &forces);

      // H2 - H1 interaction
      coulombic_force(&water_site_pos[i].H2, &water_site_pos[j].H1, &forces, q_H, q_H);
      accumulate_force_sites(&water_site_fr[i].H2, &water_site_fr[j].H1, &forces);
      // H2 - H2 interaction
      coulombic_force(&water_site_pos[i].H2, &water_site_pos[j].H2, &forces, q_H, q_H);
      accumulate_force_sites(&water_site_fr[i].H2, &water_site_fr[j].H2, &forces);
      // H1 - q1 interaction
      coulombic_force(&water_site_pos[i].H2, &water_site_pos[j].q1, &forces, q_H, q_q);
      accumulate_force_sites(&water_site_fr[i].H2, &water_site_fr[j].q1, &forces);

      // q1 - H1 interaction
      coulombic_force(&water_site_pos[i].q1, &water_site_pos[j].H1, &forces, q_q, q_H);
      accumulate_force_sites(&water_site_fr[i].q1, &water_site_fr[j].H1, &forces);
      // q1 - H2 interaction
      coulombic_force(&water_site_pos[i].q1, &water_site_pos[j].H2, &forces, q_q, q_H);
      accumulate_force_sites(&water_site_fr[i].q1, &water_site_fr[j].H2, &forces);
      // q1 - q1 interaction
      coulombic_force(&water_site_pos[i].q1, &water_site_pos[j].q1, &forces, q_q, q_q);
      accumulate_force_sites(&water_site_fr[i].q1, &water_site_fr[j].q1, &forces);

    }
  }
}

void set_CoM_force(h2o_buffer* water_molecules)
{
  unsigned molecule_no = water_molecules->n;
  water_site_forces* water_site_fr = water_molecules->water_site_fr;
  lin_dyn_x* x_lin_dyn = water_molecules->x_lin_dyn;
  lin_dyn_y* y_lin_dyn = water_molecules->y_lin_dyn;
  lin_dyn_z* z_lin_dyn = water_molecules->z_lin_dyn;

  for (int i = 0; i < molecule_no; i++)
  {
    x_lin_dyn[i].com_Fx = water_site_fr[i].O.fx + water_site_fr[i].H1.fx +
                       water_site_fr[i].H2.fx + water_site_fr[i].q1.fx;

    y_lin_dyn[i].com_Fy = water_site_fr[i].O.fy + water_site_fr[i].H1.fy +
                       water_site_fr[i].H2.fy + water_site_fr[i].q1.fy;

    z_lin_dyn[i].com_Fz = water_site_fr[i].O.fz + water_site_fr[i].H1.fz +
                       water_site_fr[i].H2.fz + water_site_fr[i].q1.fz;
  }
}

void set_CoM_force_n(h2o_buffer* water_molecules)
{
  unsigned molecule_no = water_molecules->n;
  water_site_forces* water_site_fr = water_molecules->water_site_fr;
  lin_dyn_x* x_lin_dyn = water_molecules->x_lin_dyn;
  lin_dyn_y* y_lin_dyn = water_molecules->y_lin_dyn;
  lin_dyn_z* z_lin_dyn = water_molecules->z_lin_dyn;

  for (int i = 0; i < molecule_no; i++)
  {
    x_lin_dyn[i].com_Fx_n = water_site_fr[i].O.fx + water_site_fr[i].H1.fx +
                       water_site_fr[i].H2.fx + water_site_fr[i].q1.fx;

    y_lin_dyn[i].com_Fy_n = water_site_fr[i].O.fy + water_site_fr[i].H1.fy +
                       water_site_fr[i].H2.fy + water_site_fr[i].q1.fy;

    z_lin_dyn[i].com_Fz_n = water_site_fr[i].O.fz + water_site_fr[i].H1.fz +
                       water_site_fr[i].H2.fz + water_site_fr[i].q1.fz;
  }
}

/*
  Position integration
*/

void next_position(h2o_buffer* water_molecules, double dt)
{
  double m_water = 18.01468;
  unsigned molecule_no = water_molecules->n;
  lin_dyn_x* x_lin_dyn = water_molecules->x_lin_dyn;
  lin_dyn_y* y_lin_dyn = water_molecules->y_lin_dyn;
  lin_dyn_z* z_lin_dyn = water_molecules->z_lin_dyn;

  //loop through every dimension sequentially
  for (int i = 0; i < molecule_no; i++)
  {
    x_lin_dyn[i].com_x = x_lin_dyn[i].com_x + dt*x_lin_dyn[i].com_u + pow(dt, 2)/2/m_water*x_lin_dyn[i].com_Fx;
  }

  for (int i = 0; i < molecule_no; i++)
  {
    y_lin_dyn[i].com_y = y_lin_dyn[i].com_y + dt*y_lin_dyn[i].com_v + pow(dt, 2)/2/m_water*y_lin_dyn[i].com_Fy;
  }

  for (int i = 0; i < molecule_no; i++)
  {
    z_lin_dyn[i].com_z = z_lin_dyn[i].com_z + dt*z_lin_dyn[i].com_w + pow(dt, 2)/2/m_water*z_lin_dyn[i].com_Fz;
  }
}

void next_velocity(h2o_buffer* water_molecules, double dt)
{
  double m_water = 18.01468;
  unsigned molecule_no = water_molecules->n;
  lin_dyn_x* x_lin_dyn = water_molecules->x_lin_dyn;
  lin_dyn_y* y_lin_dyn = water_molecules->y_lin_dyn;
  lin_dyn_z* z_lin_dyn = water_molecules->z_lin_dyn;

  //loop through every dimension sequentially
  for (int i = 0; i < molecule_no; i++)
  {
    x_lin_dyn[i].com_u = x_lin_dyn[i].com_u + dt/2/m_water*(x_lin_dyn[i].com_Fx + x_lin_dyn[i].com_Fx_n);
  }

  for (int i = 0; i < molecule_no; i++)
  {
    y_lin_dyn[i].com_v = y_lin_dyn[i].com_v + dt/2/m_water*(y_lin_dyn[i].com_Fy + y_lin_dyn[i].com_Fy_n);
  }

  for (int i = 0; i < molecule_no; i++)
  {
    z_lin_dyn[i].com_w = z_lin_dyn[i].com_w + dt/2/m_water*(z_lin_dyn[i].com_Fz + z_lin_dyn[i].com_Fz_n);
  }
}

void verlet_integrate(h2o_buffer* water_molecules, double dt)
{
  set_forces_sites(water_molecules);
  set_CoM_force(water_molecules);
  next_position(water_molecules, dt);

  set_forces_sites(water_molecules);
  set_CoM_force_n(water_molecules);
  next_velocity(water_molecules, dt);
}
