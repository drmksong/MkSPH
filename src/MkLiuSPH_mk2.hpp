//23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123
//--------------------------------------------------------------------------------------------
#ifndef LiuSPH_H
#define LiuSPH_H
//--------------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include <unistd.h>

#include "rvgs.h"
#include "MkArray.hpp"
#include "MkLiuParticle.hpp"
#include "MkLiuBound.hpp"
#include "MkLiuKernel.hpp"
#include "MkLiuPair.hpp"
#include "MkLiuParam.hpp"
#include "MkLiuGrid.hpp"

#define DIM 3

class MkLiuSPH
{
public:
  std::string FileName;
  std::string xvFileName;
  std::string stateFileName;
  std::string otherFile;
  std::string vp_xvFileName;
  std::string vp_stateFileName;
  std::string vp_otherFile;

  int Dim;    // check
  int NTotal; //check
  int NVirt;  //check

  int MaxInteraction; // maximum number of interaction check
  int NIac;           // number of interaction pairs   check
  int NBIac;          // number of boundaries-particles interaction pairs

  int MaxTimeStep;     // maximum number of time steps     check
  int CurrentTimeStep; // current time step    check
  double Dt;           // time increment    check
  double Rho_Ref;      // reference density to make particle pressure zero
  double Rho_Norm;     // normalized density for Norm_Density()
  double Rho_Bnd_Norm; // contribution from single particle

  MkLiuParticles LiuParticles;
  MkLiuBoundarys LiuBoundarys;
  MkLiuKernel LiuKernel;
  MkLiuPairs LiuPairs;
  MkLiuPairs LiuBndPairs;
  MkLiuParam LiuParam;
  MkLiuGrids LiuGrids;

  //  MkInt Pair_I;    // list of first partner of interaction pair     check
  //  MkInt Pair_J;    // list of second partner of interaction pair    check

public:
  MkLiuSPH();
  ~MkLiuSPH();

  bool Initialize();
  bool Initialize(int np, int nvirt, int gi, int gj, int gk);
  void Clear();

  void SetParticles();
  void SetKernel();
  void SetDt(double dt) { Dt = dt; }
  void SetDim(int dim) { Dim = dim; }
  void SetMaxTimeStep(int mt) { MaxTimeStep = mt; }
  MkLiuParam &GetLiuParam() { return LiuParam; }
  MkLiuParticles &GetLiuParticles() { return LiuParticles; }
  void SetRhoRef(double rho_ref) { Rho_Ref = rho_ref; }
  void SetSmoothLen(double sml)
  {
    for (int i = 0; i < LiuParticles.GetSize(); i++)
      LiuParticles[i].SetSmoothLen(sml);
  }

  void Draw();
  void Run();

  void Time_Integration();
  void Single_Step();

  void Ext_Force();
  void Int_Force();
  void Int_Bnd_Force();
  void Art_Heat();
  void Art_Visc();
  void Art_Drag();      //mk
  void Art_Repel();     //mk
  void Art_Bnd_Repel(); //mk
  void Av_Vel();

  double Calc_Bnd_Norm();  // mass of liu_boundaries
  void Norm_Bnd_Density(); // mk normalize function
  double Calc_Norm();      // norm density of particles
  void Norm_Density();     // mk normalize function
  void Sum_Density();
  void Con_Density();

  int Pair_Count();
  void Direct_Find();
  int Pair_Grid_Count();
  void Direct_Grid_Find();
  int Pair_Grid_Count_backup();
  void Direct_Grid_Find_backup();
  int Pair_Bnd_Count();
  void Direct_Bnd_Find();
  //  void Link_List();
  void H_Upgrade();

  void Input();  //
  void Output(); //check
  void Shock_Tube();
  void Shear_Cavity();
  void LowDensity(); // mk remove later
  void Virt_LD();    // mk remove later

  void Shake();
  void Droplet(double radius, double spacing, double rho);
  void Scale(double scale);
  void Move(double dx, double dy, double dz);

  void Virt_Part();

  void Viscosity();                         //check
  double P_ideal_gas(double rho, double u); //check
  double C_ideal_gas(double rho, double u); //check
  double P_art_water(double rho);           //check
  double C_art_water(double rho);           //check
};

#endif
