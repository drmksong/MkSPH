//23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123
//--------------------------------------------------------------------------------------------
#ifndef LiuSPH_H
#define LiuSPH_H
//--------------------------------------------------------------------------------------------
#include <c++\iostream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <glut.h>
#include <unistd.h>

#include "rvgs.h"
#include "MkDouble.h"
#include "MkInt.h"
#include "MkLiuParticle.hpp"
#include "MkLiuKernel.hpp"
#include "MkLiuPair.hpp"
#include "MkLiuParam.hpp"

#define DIM 3

class MkLiuSPH {
public:
  std::string FileName;
  std::string xvFileName;
  std::string stateFileName;
  std::string otherFile;
  std::string vp_xvFileName;
  std::string vp_stateFileName;
  std::string vp_otherFile;

  int Dim;  // check
  int NTotal;  //check
  int NVirt;  //check

  int MaxInteraction; // maximum number of interaction check
  int NIac;           // number of interaction pairs   check

  int MaxTimeStep; // maximum number of time steps     check
  int CurrentTimeStep;   // current time step    check
  double Dt;       // time increment    check
  double Rho_Ref;  // reference density to make particle pressure zero
  double Rho_Norm; // normalized density for Norm_Density()

  MkLiuParticles LiuParticles;
  MkLiuKernel LiuKernel;
  MkLiuPairs LiuPairs;
  MkLiuParam LiuParam;

  //  MkInt Pair_I;    // list of first partner of interaction pair     check
  //  MkInt Pair_J;    // list of second partner of interaction pair    check

public:

  MkLiuSPH();
  ~MkLiuSPH();

  bool Initialize();
  bool Initialize(int np, int nvirt, int gi, int gj,int gk); 
  void Clear();

  void SetParticles();
  void SetKernel();
  void SetDt(double dt){Dt = dt;}
  void SetDim(int dim){Dim = dim;}
  void SetMaxTimeStep(int mt){MaxTimeStep = mt;}
  MkLiuParam &GetLiuParam(){return LiuParam;}
  MkLiuParticles &GetLiuParticles(){return LiuParticles;}
  void SetRhoRef(double rho_ref){Rho_Ref = rho_ref;}
  void SetSmoothLen(double sml){for (int i=0;i<LiuParticles.GetSize();i++) LiuParticles[i].SetSmoothLen(sml);}

  void Draw(); 
  void Run();

  void Time_Integration();
  void Single_Step(); 

  void Ext_Force() ;
  void Int_Force() ;
  void Int_Force_MK();
  void Art_Heat();
  void Art_Visc(); 
  void Art_Drag();  //mk
  void Art_Repel(); //mk
  void Av_Vel() ;

  double Calc_Norm();
  void Norm_Density();  // mk normalize function
  void Sum_Density(); 
  void Con_Density(); 

  int Pair_Count();
  void Direct_Find();
  //  void Link_List();
  void H_Upgrade();

  void Input();  //
  void Output() ;  //check
  void Shock_Tube() ;
  void Shear_Cavity() ;
  void LowDensity();  // mk remove later
  void Virt_LD(); // mk remove later

  void Shake();
  void Droplet(double radius, double spacing, double rho);
  void Scale(double scale);
  void Move(double dx, double dy, double dz);

  void Virt_Part();

  void Viscosity();  //check
  double P_ideal_gas(double rho, double u);  //check
  double C_ideal_gas(double rho, double u);  //check
  double P_art_water(double rho);  //check
  double C_art_water(double rho);  //check
};

#endif
