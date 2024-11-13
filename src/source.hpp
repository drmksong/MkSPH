//234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#include <math.h>
#include "param.h"
#include <GL/glut.h>
#include "MkDouble.h"
#include "MkInt.h"

#define DIM 3

struct SPH_Data {
public:
  int NTotal;      // number of particles
  int MaxN;        // maximum number of particles to be accomodated 
  int MaxInteraction; // maximum number of interaction 
  int MaxTimeStep; // maximum number of time steps     
  int NIac;        // number of interaction paris
  int ITimeStep;   // current time step
  int NVirt;
  int Dim;
  MkInt CountIac; // number of neighboring particles
  double Dt;       // time increment
  MkDouble X;      // coordinate of all particles
  MkDouble VX;     // velocity of all particles 
  MkDouble AveVel; // average velocity of each particle
  MkInt IType;     // type of particles 1:ideal gas, 2:water
  MkDouble Mass;   // mass of particles
  MkDouble Rho;    // density of particles
  MkDouble Eta;    // dynamic viscosity
  MkDouble P;      // pressure
  MkDouble T;      // temperature
  MkDouble U;      // internal energy of particles
  MkDouble C;      // sound speed of particles
  MkDouble Hsml;   // smoothing lengths of particles 
  MkDouble XGCell;   // x, y and z coordinate of sorging grid cell 

  MkDouble DEDt;   // produced artificial heat, adding to energy equation
  MkDouble DRhoDt; // density change rate of each particle
  MkDouble DVXDt;  // acceleration with respect to x, y and z
  MkDouble TDSDt;  // production of viscous entroyp t * ds/dt
  MkDouble DUDt;   // du/dt

  MkInt Pair_I;    // list of first partner of interaction pair 
  MkInt Pair_J;    // list of second partner of interaction pair
  MkDouble W;      // kernal of for all interaction pairs
  MkDouble DWDX;   // derivative of kernel with respect to x, y and z

public:
  SPH_Data();
  ~SPH_Data();
  void Initialize(int maxn,int dim);
  void Clear();
  void Dump();
  void Draw();
} ;


struct SPH_Grid {
public:
  int Dim;
  int MaxNGridX;
  int MaxNGridY;
  int MaxNGridZ;
  MkInt NGridX;      // number of sorting grid cells in x,y and z direction
  MkDouble MaxGridX; // maximum x, y and z coordinate of grid range
  MkDouble MinGridX; // minimum x, y and z coordinate of grid range
  MkDouble DGeomX;   // x, y and z expansion of grid range
  MkDouble GHsmlX;   // smoothing length measure in cells of the grid
  MkInt Grid;           // DIM dimension grid

public:
  SPH_Grid();
  ~SPH_Grid();
  void Initialize(int dim, int nx, int ny, int nz);
  void Initialize(int dim);
  void Clear();
  void Dump();
};

void art_heat(SPH_Data &data);
void art_visc(SPH_Data &data); 
void av_vel(SPH_Data &data) ;
void sum_density(SPH_Data &data); 
void con_density(SPH_Data &data); 
void direct_find(SPH_Data &data);
void ext_force(SPH_Data &data) ;
void h_upgrade(SPH_Data &data);
void time_integration(SPH_Data &data);
void viscosity(SPH_Data &data);
void single_step(SPH_Data &data); 
void virt_part(SPH_Data &data);
void output(SPH_Data &data) ; 
void link_list(SPH_Data &data);
void input(SPH_Data &data);
void shock_tube(SPH_Data &data) ;
void shear_cavity(SPH_Data &data) ;
void int_force(SPH_Data &data) ;

void init_grid(SPH_Data &data, SPH_Grid &grid);

void p_gas(double rho, double u, double &p, double &c) ;
void p_art_water(double rho, double &p, double &c) ;
void grid_geom(int i,MkDouble x,int ngridx,MkDouble maxgridx,MkDouble mingridx,MkDouble dgeomx,MkInt & xgcell,SPH_Data &data);
void kernel(double r,MkDouble dx,double hsml,double &w, MkDouble &dwdx,SPH_Data &data);   

void drawBox(void);
void display(void);
void init(void);
void update2(void);
int mainfunc(void);

extern  SPH_Data data;








