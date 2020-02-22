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
  void SetSmoothLen(double sml){for (int i=0;i<LiuParticles.GetSize();i++) LiuParticles[i].SetSmoothLen(sml);}

  void Draw(); 
  void Run();

  void Time_Integration();
  void Single_Step(); 

  void Ext_Force() ;
  void Int_Force() ;
  void Art_Heat();
  void Art_Visc(); 
  void Av_Vel() ;

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

  void Shake();

  void Virt_Part();

  void Viscosity();  //check
  double P_ideal_gas(double rho, double u);  //check
  double C_ideal_gas(double rho, double u);  //check
  double P_art_water(double rho);  //check
  double C_art_water(double rho);  //check
};

#endif
