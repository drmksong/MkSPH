//---------------------------------------------------------------------------
#ifndef MkStrainH
#define MkStrainH

#include <stdio.h>
//#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include "MkDouble.h"
#include "MkInt.h"
#include "MkMisc.h"

#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif

#define M_PI 3.14159265358979323846

class MkStrain
{
public:
  int StrainType; // strainrate = 0, strain = 1, spinrate=2, spin=3, plasticstrainrate = 4, plasticstrain = 5
  double FStrain[3][3];
  double FMajorStrain;
  double FDeviatoric[3][3];

public:
  MkStrain();
  MkStrain(double strain[3][3]);
  void Identity();
  void LoadIdentity() { Identity(); }
  void Clear();
  void CalcMajorStrain();
  void CalcDeviatoric();
  void Update(MkStrain &srate, double dt);
  void SetType(int type) { StrainType = type; }

  MkStrain &operator=(const MkStrain &str);
  double &operator()(int i, int j)
  {
    static double Zero;
    Zero = 0;
    if (i < 3 && i >= 0 && j < 3 && j >= 0)
      return FStrain[i][j];
    else
      exit(-11);
  }

  void Out(char *);
  void Out();
};

extern MkStrain NullStrain;
//---------------------------------------------------------------------------
#endif
