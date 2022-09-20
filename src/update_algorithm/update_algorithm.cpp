#include "update_algorithm.hpp"
#include "../debug/debug.hpp"
#include "../quaternions/quaternions.hpp"
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

void reset_forces_site(site_forces* site)
{
  site->fx = 0;
  site->fy = 0;
  site->fz = 0;
}

void set_forces_sites(h2o_buffer* water_molecules)
{
  //physics constants
  double sigma = 1;//3.1589; //angstroms
  double epsilon = 1;// 0.1852;

  double q_H = 1;//0.52;
  double q_q = 1; //-1.04
  //inlined by compiler
  unsigned molecule_no = water_molecules->n;
  water_site_positions* water_site_pos = water_molecules->water_site_pos;
  water_site_forces* water_site_fr = water_molecules->water_site_fr;
  //struct to store temp forces results
  site_forces forces;
  //reset forces
  for (int i = 0; i < molecule_no; i++)
  {
    reset_forces_site(&water_site_fr[i].O);
    reset_forces_site(&water_site_fr[i].H1);
    reset_forces_site(&water_site_fr[i].H2);
    reset_forces_site(&water_site_fr[i].q1);
  }
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
      // H2 - q1 interaction
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
  double m_water = 1;// 18.01468;
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
  double m_water = 1;// 18.01468;
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
  water_molecules->site_global_coordiantes();

  set_forces_sites(water_molecules);
  set_CoM_force_n(water_molecules);
  next_velocity(water_molecules, dt);
  water_molecules->site_global_coordiantes();
}

/*
  Orientation integration
*/

void compute_torques(h2o_buffer* water_molecules)
{
  unsigned molecule_no = water_molecules->n;
  water_site_forces* water_site_fr = water_molecules->water_site_fr;
  quaternion_orientation* orientations = water_molecules->orientations;
  torques* water_torques = water_molecules->water_torques;
  //need to update to compute torque about each direction separately
  //transformation offset
  double zero[3] = {0, 0, 0};
  //memory to store global site force
  double force_global[3];

  //site local coordinates
  double O_L[3] = {0, -0.065555, 0};
  double H1_L[3] = {-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};
  double H2_L[3] = {0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};
  double q1_L[3] = {0, 0.15 - 0.065555, 0};

  //memory for forces
  double local_O[3];
  double local_H1[3];
  double local_H2[3];
  double local_q1[3];

  //memory for torques
  double torque_O[3];
  double torque_H1[3];
  double torque_H2[3];
  double torque_q1[3];

  for (int i = 0; i < molecule_no; i++)
  {
    //transform to local coordinates

    //O
    force_global[0] = water_site_fr[i].O.fx;
    force_global[1] = water_site_fr[i].O.fy;
    force_global[2] = water_site_fr[i].O.fz;

    orientations[i].current.transform_vector_invert(force_global, zero, local_O);

    //H1
    force_global[0] = water_site_fr[i].H1.fx;
    force_global[1] = water_site_fr[i].H1.fy;
    force_global[2] = water_site_fr[i].H1.fz;

    orientations[i].current.transform_vector_invert(force_global, zero, local_H1);

    //H2
    force_global[0] = water_site_fr[i].H2.fx;
    force_global[1] = water_site_fr[i].H2.fy;
    force_global[2] = water_site_fr[i].H2.fz;

    orientations[i].current.transform_vector_invert(force_global, zero, local_H2);

    //q1
    force_global[0] = water_site_fr[i].q1.fx;
    force_global[1] = water_site_fr[i].q1.fy;
    force_global[2] = water_site_fr[i].q1.fz;

    orientations[i].current.transform_vector_invert(force_global, zero, local_q1);

    //compute torques
    cross_product(O_L, local_O, torque_O);
    cross_product(H1_L, local_H1, torque_H1);
    cross_product(H2_L, local_H2, torque_H2);
    cross_product(q1_L, local_q1, torque_q1);

    //sum torques
    water_torques[i].tx = torque_O[0]+ torque_H1[0] + torque_H2[0] + torque_q1[0];
    water_torques[i].ty = torque_O[1]+ torque_H1[1] + torque_H2[1] + torque_q1[1];
    water_torques[i].tz = torque_O[2]+ torque_H1[2] + torque_H2[2] + torque_q1[2];
  }
}

void compute_next_torques(h2o_buffer* water_molecules)
{
  unsigned molecule_no = water_molecules->n;
  water_site_forces* water_site_fr = water_molecules->water_site_fr;
  quaternion_orientation* orientations = water_molecules->orientations;
  torques* water_torques = water_molecules->water_torques;
  //need to update to compute torque about each direction separately
  //transformation offset
  double zero[3] = {0, 0, 0};
  //memory to store global site force
  double force_global[3];

  //site local coordinates
  double O_L[3] = {0, -0.065555, 0};
  double H1_L[3] = {-0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};
  double H2_L[3] = {0.9572*sin(52.26*M_PI/180), 0.9572*cos(52.26*M_PI/180) - 0.065555, 0};
  double q1_L[3] = {0, 0.15 - 0.065555, 0};

  //memory for forces
  double local_O[3];
  double local_H1[3];
  double local_H2[3];
  double local_q1[3];

  //memory for torques
  double torque_O[3];
  double torque_H1[3];
  double torque_H2[3];
  double torque_q1[3];

  for (int i = 0; i < molecule_no; i++)
  {
    //transform to local coordinates

    //O
    force_global[0] = water_site_fr[i].O.fx;
    force_global[1] = water_site_fr[i].O.fy;
    force_global[2] = water_site_fr[i].O.fz;

    orientations[i].next.transform_vector_invert(force_global, zero, local_O);

    //H1
    force_global[0] = water_site_fr[i].H1.fx;
    force_global[1] = water_site_fr[i].H1.fy;
    force_global[2] = water_site_fr[i].H1.fz;

    orientations[i].next.transform_vector_invert(force_global, zero, local_H1);

    //H2
    force_global[0] = water_site_fr[i].H2.fx;
    force_global[1] = water_site_fr[i].H2.fy;
    force_global[2] = water_site_fr[i].H2.fz;

    orientations[i].next.transform_vector_invert(force_global, zero, local_H2);

    //q1
    force_global[0] = water_site_fr[i].q1.fx;
    force_global[1] = water_site_fr[i].q1.fy;
    force_global[2] = water_site_fr[i].q1.fz;

    orientations[i].next.transform_vector_invert(force_global, zero, local_q1);

    //compute torques
    cross_product(O_L, local_O, torque_O);
    cross_product(H1_L, local_H1, torque_H1);
    cross_product(H2_L, local_H2, torque_H2);
    cross_product(q1_L, local_q1, torque_q1);

    //sum torques
    water_torques[i].tx = torque_O[0]+ torque_H1[0] + torque_H2[0] + torque_q1[0];
    water_torques[i].ty = torque_O[1]+ torque_H1[1] + torque_H2[1] + torque_q1[1];
    water_torques[i].tz = torque_O[2]+ torque_H1[2] + torque_H2[2] + torque_q1[2];
  }
}

void compute_angular_velocity(quaternion orientation, quaternion velocity, double omega[3])
{
  double q[4] = {orientation.q0, orientation.q[0], orientation.q[1], orientation.q[2]};
  double q_dot[4] = {velocity.q0, velocity.q[0], velocity.q[1], velocity.q[2]};

  omega[0] = -q[1]*q_dot[0] + q[0]*q_dot[1] + q[3]*q_dot[2] - q[2]*q_dot[3];
  omega[1] = -q[2]*q_dot[0] - q[3]*q_dot[1] + q[0]*q_dot[2] + q[1]*q_dot[3];
  omega[2] = -q[3]*q_dot[0] + q[3]*q_dot[1] - q[1]*q_dot[2] + q[0]*q_dot[3];
}

void compute_quat_accel(h2o_buffer* water_molecules)
{
  unsigned molecule_no = water_molecules->n;
  quaternion_orientation* orientations = water_molecules->orientations;
  quaternion_derivatives* quat_derivatives = water_molecules->quat_derivatives;
  torques* water_torques = water_molecules->water_torques;

  double I[3] = {1, 1, 1};
  double q[4];
  double omega[3];
  double v[4];
  for (int i = 0; i < molecule_no; i++)
  {
    //compute angular velocity
    compute_angular_velocity(orientations[i].current, quat_derivatives[i].vel, omega);
    //define rhs vector
    v[0] = -2*pow(orientations[i].current.get_norm(), 2);
    v[1] = water_torques[i].tx + omega[1]*omega[2]*(I[1] - I[2])/I[0];
    v[2] = water_torques[i].ty + omega[2]*omega[0]*(I[2] - I[0])/I[1];
    v[3] = water_torques[i].tz + omega[0]*omega[1]*(I[0] - I[2])/I[2];
    //define q
    q[0] = orientations[i].current.q0;
    q[1] = orientations[i].current.q[0];
    q[2] = orientations[i].current.q[1];
    q[3] = orientations[i].current.q[2];
    //compute quaternion accelerations
    (quat_derivatives[i].accel).q0 = 0.5*(q[0]*v[0] - q[1]*v[1] - q[2]*v[2] - q[3]*v[3]);
    (quat_derivatives[i].accel).q[0] = 0.5*(q[1]*v[0] + q[0]*v[1] - q[3]*v[2] + q[2]*v[3]);
    (quat_derivatives[i].accel).q[1] = 0.5*(q[2]*v[0] + q[3]*v[1] + q[0]*v[2] - q[1]*v[3]);
    (quat_derivatives[i].accel).q[2] = 0.5*(q[3]*v[0] - q[2]*v[1] + q[1]*v[2] + q[0]*v[3]);
  }
}

void compute_next_quat_accel(h2o_buffer* water_molecules)
{

}

void compute_quat_velocity(h2o_buffer* water_molecules, double dt)
{
  unsigned molecule_no = water_molecules->n;
  quaternion_derivatives* quat_derivatives = water_molecules->quat_derivatives;
  quaternion_orientation* orientations = water_molecules->orientations;

  for (int i = 0; i < molecule_no; i++)
  {
    quat_derivatives[i].vel = quat_derivatives[i].vel + (quat_derivatives[i].accel).multiply(0.5*dt);
  }
}

void next_orientation(h2o_buffer* water_molecules, double dt)
{
  unsigned molecule_no = water_molecules->n;
  quaternion_orientation* orientations = water_molecules->orientations;
  quaternion_derivatives* quat_derivatives = water_molecules->quat_derivatives;

  for (int i = 0; i < molecule_no; i++)
  {
    orientations[i].next = orientations[i].current + (quat_derivatives[i].vel).multiply(dt);
  }
}

void converge_quat_velocity(h2o_buffer* water_molecules, double dt)
{
  unsigned molecule_no = water_molecules->n;
  quaternion_orientation* orientations = water_molecules->orientations;
  quaternion_derivatives* quat_derivatives = water_molecules->quat_derivatives;
  torques* water_torques = water_molecules->water_torques;
  //angular velocity guess
  double I[3] = {1, 1, 1};
  double omega[3];
  double next_omega[3];
  double v[3];
  double q[3];
  double condition = 1;
  //iterate through every molecule
  for (int i = 0; i < molecule_no; i++)
  {
    //average quaternion orientaions
    double q[4] = {0.5*(orientations[i].current.q0 + orientations[i].next.q0),
                   0.5*(orientations[i].current.q[0] + orientations[i].next.q[0]),
                   0.5*(orientations[i].current.q[1] + orientations[i].next.q[1]),
                   0.5*(orientations[i].current.q[2] + orientations[i].next.q[2])
                  };

    double q_dot[4] = {quat_derivatives[i].vel.q0, quat_derivatives[i].vel.q[0], quat_derivatives[i].vel.q[1], quat_derivatives[i].vel.q[2]};

    omega[0] = -q[1]*q_dot[0] + q[0]*q_dot[1] + q[3]*q_dot[2] - q[2]*q_dot[3];
    omega[1] = -q[2]*q_dot[0] - q[3]*q_dot[1] + q[0]*q_dot[2] + q[1]*q_dot[3];
    omega[2] = -q[3]*q_dot[0] + q[3]*q_dot[1] - q[1]*q_dot[2] + q[0]*q_dot[3];

    while(condition > 1E-3)
    { //reset guess
      omega[0] = next_omega[0];
      omega[1] = next_omega[1];
      omega[2] = next_omega[2];
      //calculate new angular acceleration
      v[0] = -2*pow(orientations[i].next.get_norm(), 2);
      v[1] = water_torques[i].tx + omega[1]*omega[2]*(I[1] - I[2])/I[0];
      v[2] = water_torques[i].ty + omega[2]*omega[0]*(I[2] - I[0])/I[1];
      v[3] = water_torques[i].tz + omega[0]*omega[1]*(I[0] - I[2])/I[2];
      //define q
      q[0] = orientations[i].next.q0;
      q[1] = orientations[i].next.q[0];
      q[2] = orientations[i].next.q[1];
      q[3] = orientations[i].next.q[2];
      //compute quaternion accelerations
      (quat_derivatives[i].accel).q0 = 0.5*(q[0]*v[0] - q[1]*v[1] - q[2]*v[2] - q[3]*v[3]);
      (quat_derivatives[i].accel).q[0] = 0.5*(q[1]*v[0] + q[0]*v[1] - q[3]*v[2] + q[2]*v[3]);
      (quat_derivatives[i].accel).q[1] = 0.5*(q[2]*v[0] + q[3]*v[1] + q[0]*v[2] - q[1]*v[3]);
      (quat_derivatives[i].accel).q[2] = 0.5*(q[3]*v[0] - q[2]*v[1] + q[1]*v[2] + q[0]*v[3]);
      //update quaternion velocity
      quat_derivatives[i].vel = quat_derivatives[i].vel + (quat_derivatives[i].accel).multiply(0.5*dt);
      //calculate angular velocity
      compute_angular_velocity(orientations[i].next, quat_derivatives[i].vel, next_omega);
      //calculate convergence condition
      condition = sqrt(pow(omega[0] - next_omega[0], 2) + pow(omega[1] - next_omega[1], 2) + pow(omega[2] - next_omega[2], 2));
    }
  }
}

void reset_orientation(h2o_buffer* water_molecules)
{
  unsigned molecule_no = water_molecules->n;
  quaternion_orientation* orientations = water_molecules->orientations;

  for (int i = 0; i < molecule_no; i++)
  {
    orientations[i].current.q0 = orientations[i].next.q0;
    orientations[i].current.q[0] = orientations[i].next.q[0];
    orientations[i].current.q[1] = orientations[i].next.q[1];
    orientations[i].current.q[2] = orientations[i].next.q[2];
  }
}

void cross_product(double vector_a[3], double vector_b[3], double output[3])
{
  output[0] = vector_a[1]*vector_b[2] - vector_a[2]*vector_b[1];
  output[1] = vector_a[2]*vector_b[0] - vector_a[0]*vector_b[2];
  output[2] = vector_a[0]*vector_b[1] - vector_a[1]*vector_b[0];
}

double torque_magnitude(torques* torque)
{
  return pow(pow(torque->tx, 2) + pow(torque->ty, 2) + pow(torque->tz, 2), 0.5);
}

void quaternion_verlet_integrate(h2o_buffer* water_molecules, double dt)
{
  //Step 1
  //calculate torques.
  compute_torques(water_molecules);
  //calculate quaternion accel
  compute_quat_accel(water_molecules);
  //caluclate half step quaternion velocity
  compute_quat_velocity(water_molecules, dt);
  //calculate new quaternion orientations
  next_orientation(water_molecules, dt);

  //Step 2
  //re-orientate molecules and calculate new site forces
  water_molecules->next_site_global_coordinates();
  set_forces_sites(water_molecules);
  //calculate torques at new time step
  compute_next_torques(water_molecules);
  //iterate until quaternion velcity at n+1 converges
  converge_quat_velocity(water_molecules, dt);

  //step 3
  //reset starting parameters
  void reset_orientation(h2o_buffer* water_molecules);

}
