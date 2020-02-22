//--------------------------------------------------------- 
//     Including file for parameters and constants used  
//     in the entire SPH software packages. 
//--------------------------------------------------------- 
 
#ifndef LiuPARAM_H
#define LiuPARAM_H

class MkLiuParam {
public:
int dim; 
int maxn;
int max_interaction; 
double x_maxgeom;
double x_mingeom; 
double y_maxgeom;  
double y_mingeom;
double z_maxgeom;
double z_mingeom;
int pa_sph; 
int nnps;
int sle; 
int skf; 
bool summation_density; 
bool average_velocity; 
bool config_input; 
bool virtual_part; 
bool vp_input; 
bool visc;
bool ex_force; 
bool visc_artificial; 
bool heat_artificial; 
bool self_gravity;       
bool nor_density;       
int  nsym; 
bool int_stat; 

int print_step; 
int save_step;
int moni_particle;
            
double pi;
bool shocktube; 
bool shearcavity;
bool lowdensity;

public:
  MkLiuParam(){Initialize();}
  ~MkLiuParam(){}
  void Initialize();
};

#endif
