//---------------------------------------------------------------------------
#include "MkStrain.h"

//---------------------------------------------------------------------------
// index
//  [0][0] [0][1] [0][2] 
//  [1][0] [1][1] [1][2] 
//  [2][0] [2][1] [2][2] 
//---------------------------------------------------------------------------
MkStrain NullStrain;

MkStrain::MkStrain()
{
   for (int i=0;i<3;i++)
       for (int j=0;j<3;j++)
           FStrain[i][j] = 0;
}

MkStrain::MkStrain(double strain[3][3])
{
   for (int i=0;i<3;i++)
       for (int j=0;j<3;j++)
           FStrain[i][j] = strain[i][j];
}

void MkStrain::Identity()
{
   int i;
   for (i=0;i<3;i++)
       for (int j=0;j<3;j++)
           FStrain[i][j] = 0;

   for (i=0;i<3;i++)
       FStrain[i][i] = 1;
}

void MkStrain::Clear()
{
   for (int i=0;i<3;i++)
       for (int j=0;j<3;j++)
           FStrain[i][j] = 0;
}

void MkStrain::CalcMajorStrain()
{
  FMajorStrain = FStrain[0][0]+FStrain[1][1]+FStrain[2][2];
}

void MkStrain::CalcDeviatoric()
{
  CalcMajorStrain();
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      FDeviatoric[i][j] = FStrain[i][j];
      if (i==j) FDeviatoric[i][j] -= FMajorStrain/3.0;
    }
  }
}

void MkStrain::Update(MkStrain &srate, double dt)
{
  if (srate.StrainType == 1 || srate.StrainType == 3 || srate.StrainType == 5) {
    MkDebug("MkStrain::Update trying to update with strain, instead of strainrate\n");
    exit(-11);
  }
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++){
      FStrain[i][j] += dt*srate.FStrain[i][j];
    }
  }
  CalcDeviatoric();  
} 

MkStrain & MkStrain::operator=(const MkStrain &str)
{
  StrainType = str.StrainType;
  FMajorStrain = str.FMajorStrain ;

  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++){
      FStrain[i][j] = str.FStrain[i][j];
    }
  }
  
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++){
      FDeviatoric[i][j] = str.FDeviatoric[i][j];
    }
  }
  
  return (*this);
}

void MkStrain::Out(char *fname)
{
  FILE *fp;
  char str[256],s[256];

  fp = fopen(fname,"a");
  fprintf(fp,"Output of Strain ");
  sprintf(str,"A-th row [%10.3f %10.3f %10.3f %10.3f ]",0.0,1.0,2.0,3.0);
  fprintf(fp,str);
  for (int i=0;i<4;i++) {
    sprintf(str,"%d-th row [",i);
    for (int j=0;j<3;j++) {
      sprintf(s,"%10.3f ",FStrain[i][j]);
      strcat(str,s);
    }
    sprintf(s,"]");

    strcat(str,s);
    fprintf(fp,str);
  }
  fclose(fp);
}

void MkStrain::Out()
{
  char str[256],s[256];

  printf("Output of Strain \n");
  sprintf(str,"A-th row [%12d %12d %12d ]\n",0,1,2);
  printf(str);
  for (int i=0;i<3;i++) {
    sprintf(str,"%d-th row [",i);
    for (int j=0;j<3;j++) {
      sprintf(s,"%12.7f ",FStrain[i][j]);
      strcat(str,s);
    }
    sprintf(s,"]\n");
    strcat(str,s);
    printf(str);
  }
}

