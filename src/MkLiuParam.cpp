#include "MkLiuParam.hpp"

void MkLiuParam::Initialize()
{
  //     dim : Dimension of the problem (1, 2 or 3)
  dim = 3;

  //     maxn    : Maximum number of particles
  //     max_interation : Maximum number of interaction pairs
  maxn = 12000;
  max_interaction = 100 * maxn;

  //     Parameters for the computational geometry,
  //     x_maxgeom : Upper limit of allowed x-regime
  //     x_mingeom : Lower limit of allowed x-regime
  //     y_maxgeom : Upper limit of allowed y-regime
  //     y_mingeom : Lower limit of allowed y-regime
  //     z_maxgeom : Upper limit of allowed z-regime
  //     z_mingeom : Lower limit of allowed z-regime

  x_maxgeom = 10.;
  x_mingeom = -10.e0;
  y_maxgeom = 10.e0;
  y_mingeom = -10.e0;
  z_maxgeom = 10.e0;
  z_mingeom = -10.e0;

  //     SPH algorithm for particle approximation (pa_sph)
  //     pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
  //              2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
  pa_sph = 1;

  //     Nearest neighbor particle searching (nnps) method
  //     nnps = 1 : Simplest and direct searching
  //            2 : Sorting grid linked list
  //            3 : Tree algorithm
  nnps = 1;

  //     Smoothing length evolution (sle) algorithm
  //     sle = 0 : Keep unchanged,
  //           1 : h = fac * (m/rho)^(1/dim)
  //           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
  //           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )

  sle = 0;

  //     Smoothing kernel function
  //     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
  //         = 2, Gauss kernel   (Gingold and Monaghan 1981)
  //         = 3, Quintic kernel (Morris 1997)
  skf = 1;

  //     Switches for different senarios

  //     summation_density = .TRUE. : Use density summation model in the code,
  //                        .FALSE.: Use continuiity equation
  //     average_velocity = .TRUE. : Monaghan treatment on average velocity,
  //                       .FALSE.: No average treatment.
  //     config_input = .TRUE. : Load initial configuration data,
  //                   .FALSE.: Generate initial configuration.
  //     virtual_part = .TRUE. : Use vritual particle,
  //                   .FALSE.: No use of vritual particle.
  //     vp_input = .TRUE. : Load virtual particle information,
  //               .FALSE.: Generate virtual particle information.
  //     visc = .true. : Consider viscosity,
  //           .false.: No viscosity.
  //     ex_force =.true. : Consider external force,
  //               .false.: No external force.
  //     visc_artificial = .true. : Consider artificial viscosity,
  //                      .false.: No considering of artificial viscosity.
  //     heat_artificial = .true. : Consider artificial heating,
  //                      .false.: No considering of artificial heating.
  //     self_gravity = .true. : Considering self_gravity,
  //                    .false.: No considering of self_gravity
  //     nor_density =  .true. : Density normalization by using CSPM,
  //                    .false.: No normalization.

  summation_density = true;
  average_velocity = false;
  config_input = false;
  virtual_part = false; //true;
  vp_input = false;
  visc = true;
  ex_force = true;
  visc_artificial = true;
  heat_artificial = false;
  self_gravity = true;
  nor_density = false;

  //     Symmetry of the problem
  //     nsym = 0 : no symmetry,
  //          = 1 : axis symmetry,
  //          = 2 : center symmetry.
  nsym = 0;

  //     Control parameters for output
  //     int_stat = .true. : Pr statistics about SPH particle interactions.
  //                        including virtual particle information.
  //     print_step: Pr Timestep (On Screen)
  //     save_step : Save Timestep    (To Disk File)
  //     moni_particle: The particle number for information monitoring.

  int_stat = true;

  print_step = 500;
  save_step = 5000;
  moni_particle = 1600;

  pi = 3.14159265358979323846;

  //     Simulation cases
  //     shocktube = .true. : carry out shock tube simulation
  //     shearcavity = .true. : carry out shear cavity simulation

  shocktube = true;
  shearcavity = false; //true;
  lowdensity = false;
}
