//---------------------------------------------------------------------------
#include "MkFailure.h"

//---------------------------------------------------------------------------
// index
//  [0][0] [0][1] [0][2] 
//  [1][0] [1][1] [1][2] 
//  [2][0] [2][1] [2][2] 
//---------------------------------------------------------------------------
MkFailure NullFailure;

MkFailure::MkFailure()
{
  FAlpha = 0;
  FKay = 0;
  FTensionStatus = false;
  FShearStatus = false;
}

MkFailure::MkFailure(double f, double c) // f is degree
{
  double t = tan(f*3.141592654/180.0);
  FAlpha = t / sqrt(9+12*t);
  FKay = 3*c/ sqrt(9+12*t);
}


void MkFailure::Clear()
{
  FAlpha = 0;
  FKay = 0;
  FTensionStatus = false;
  FShearStatus = false;
}

void MkFailure::Set(double f, double c) // f is degree
{
  double t = tan(f*3.141592654/180.0);
  FAlpha = t / sqrt(9+12*t);
  FKay = 3*c / sqrt(9+12*t);
}

bool MkFailure::CheckTension(MkStress &str) // str.CalcInvariant should be done already
{
  double g;
  g = FAlpha*str.FI_1-FKay;
  if (0<=g) FTensionStatus = true;
  else FTensionStatus = false;
  return FTensionStatus;
}

bool MkFailure::CheckShear(MkStress &str) // str.CalcInvariant should be done already
{
  double g;
  g = sqrt(str.FI_2)+FAlpha*str.FI_1-FKay;
  if (0<=g) FShearStatus = true;
  else FShearStatus = false;
  return FShearStatus;
}

void MkFailure::TreatTension(MkStress &str) // str.CalcInvariant should be done already
{
  double ka = FKay/FAlpha;
  str.FStress[0][0] -= 1.0/3.0*(str.FI_1 - ka);
  str.FStress[1][1] -= 1.0/3.0*(str.FI_1 - ka);
  str.FStress[2][2] -= 1.0/3.0*(str.FI_1 - ka);
}

void MkFailure::TreatScaling(MkStress &str) // str.CalcInvariant should be done already
{
  double rn = -FAlpha*str.FI_1+FKay/sqrt(str.FJ_2);
  str.FStress[0][0] = rn*str.FDeviatoric[0][0] + str.FI_1/3.0;  
  str.FStress[1][1] = rn*str.FDeviatoric[1][1] + str.FI_1/3.0;  
  str.FStress[2][2] = rn*str.FDeviatoric[2][2] + str.FI_1/3.0;  
  str.FStress[0][1] = rn*str.FDeviatoric[0][1];  
  str.FStress[1][0] =    str.FStress[0][1]; 
  str.FStress[0][2] = rn*str.FDeviatoric[0][2];  
  str.FStress[2][0] =    str.FStress[2][0]; 
  str.FStress[1][2] = rn*str.FDeviatoric[1][2];  
  str.FStress[2][1] =    str.FStress[2][1]; 
}

MkFailure & MkFailure::operator=(const MkFailure &failure)
{
  FAlpha = failure.FAlpha;
  FKay = failure.FKay;
  FTensionStatus = failure.FTensionStatus;
  FShearStatus = failure.FShearStatus;

  return *this;
}

void MkFailure::Out(char *fname)
{
  FILE *fp;
  char str[256],s[256];

  fp = fopen(fname,"a");
  fclose(fp);
}

void MkFailure::Out()
{
  /*  char str[256],s[256];

  printf("Output of Failure \n");
  sprintf(str,"A-th row [%10.5f %10.5f %10.5f ]\n",0.0,1.0,2.0);
  printf(str);
  for (int i=0;i<3;i++) {
    sprintf(str,"%d-th row [",i);
    for (int j=0;j<3;j++) {
      sprintf(s,"%10.5f ",FFailure[i][j]);
      strcat(str,s);
    }
    sprintf(s,"]\n");
    strcat(str,s);
    printf(str);
    }*/
}

