//--------------------------------------------------------- 
//     Including file for parameters and constants used  
//     in the entire SPH software packages. 
//--------------------------------------------------------- 
 
#ifndef PARAM_H
#define PARAM_H

//     dim : Dimension of the problem (1, 2 or 3) 
int dim  = 3; 
 
//     maxn    : Maximum number of particles 
//     max_interation : Maximum number of interaction pairs 
int maxn = 12000;
int max_interaction = 100 * maxn; 
 
//     Parameters for the computational geometry,   
//     x_maxgeom : Upper limit of allowed x-regime  
//     x_mingeom : Lower limit of allowed x-regime  
//     y_maxgeom : Upper limit of allowed y-regime  
//     y_mingeom : Lower limit of allowed y-regime  
//     z_maxgeom : Upper limit of allowed z-regime  
//     z_mingeom : Lower limit of allowed z-regime  

double x_maxgeom =  10.;
double x_mingeom = -10.e0; 
double y_maxgeom =  10.e0;  
double y_mingeom = -10.e0;
double z_maxgeom =  10.e0;
double z_mingeom = -10.e0;
     
//     SPH algorithm for particle approximation (pa_sph) 
//     pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j)) 
//              2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2) 
int pa_sph = 2; 
 
//     Nearest neighbor particle searching (nnps) method 
//     nnps = 1 : Simplest and direct searching 
//            2 : Sorting grid linked list 
//            3 : Tree algorithm 
int nnps  = 1;
 
//     Smoothing length evolution (sle) algorithm 
//     sle = 0 : Keep unchanged, 
//           1 : h = fac * (m/rho)^(1/dim) 
//           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt) 
//           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )  
 
int sle = 0; 
 
//     Smoothing kernel function  
//     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//         = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//         = 3, Quintic kernel (Morris 1997) 
int skf= 1; 
 
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

bool summation_density  = true; 
bool average_velocity  = true; 
bool config_input  = false; 
bool virtual_part  = true; 
bool vp_input  = false; 
bool visc  = true;
bool ex_force  = true; 
bool visc_artificial  = false; 
bool heat_artificial  = false; 
bool self_gravity  = false;       
bool nor_density  = true;       
 
//     Symmetry of the problem 
//     nsym = 0 : no symmetry, 
//          = 1 : axis symmetry, 
//          = 2 : center symmetry.      
int    nsym = 0; 
 
//     Control parameters for output  
//     int_stat = .true. : Print statistics about SPH particle interactions. 
//                        including virtual particle information. 
//     print_step: Print Timestep (On Screen) 
//     save_step : Save Timestep    (To Disk File) 
//     moni_particle: The particle number for information monitoring. 

bool int_stat = true; 

int print_step = 100; 
int save_step = 500;
int moni_particle = 1600;
            
double pi  = 3.14159265358979323846;
 
//     Simulation cases 
//     shocktube = .true. : carry out shock tube simulation 
//     shearcavity = .true. : carry out shear cavity simulation 

bool shocktube  = false; 
bool shearcavity  = true;

#endif
