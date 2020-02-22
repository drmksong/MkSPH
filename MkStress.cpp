//---------------------------------------------------------------------------
#include "MkStress.h"

//---------------------------------------------------------------------------
// index
//  [0][0] [0][1] [0][2] 
//  [1][0] [1][1] [1][2] 
//  [2][0] [2][1] [2][2] 
//---------------------------------------------------------------------------
MkStress NullStress;

MkStress::MkStress()
{
  FMeanStress = 0;
  FI_1 = 0; 
  FI_2 = 0; 
  FI_3 = 0;
  FJ_2 = 0;
  FPhi = 0;

  for (int i=0;i<3;i++) {
    FPrincipal[i] = 0;
    for (int j=0;j<3;j++) {
      FStress[i][j] = 0;
      FDeviatoric[i][j] = 0;
      LMN[i][j] = 0;
    }
  }
}

MkStress::MkStress(double stress[3][3])
{
  FMeanStress = 0;
  FI_1 = 0; 
  FI_2 = 0; 
  FI_3 = 0;
  FJ_2 = 0;
  FPhi = 0;

  for (int i=0;i<3;i++) {
    FPrincipal[i] = 0;
    for (int j=0;j<3;j++) {
      FStress[i][j] = stress[i][j];
      FDeviatoric[i][j] = 0;
      LMN[i][j] = 0;
    }
  }
  CalcInvariant();
  CalcPrincipal();
  CalcLMN();
}

void MkStress::Identity()
{
   int i;
   for (i=0;i<3;i++)
       for (int j=0;j<3;j++)
           FStress[i][j] = 0;

   for (i=0;i<3;i++)
       FStress[i][i] = 1;
}

void MkStress::Clear()
{
  FMeanStress = 0;
  FI_1 = 0; 
  FI_2 = 0; 
  FI_3 = 0;
  FJ_2 = 0;
  FPhi = 0;

  for (int i=0;i<3;i++) {
    FPrincipal[i] = 0;
    for (int j=0;j<3;j++) {
      FStress[i][j] = 0;
      FDeviatoric[i][j] = 0;
      LMN[i][j] = 0;
    }
  }
}

void MkStress::CalcMeanStress()
{
  CalcFI_1();
  FMeanStress = FI_1/3.0;
}

void MkStress::CalcDeviatoric()
{
  CalcMeanStress();
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      FDeviatoric[i][j] = FStress[i][j] - (i==j)?FMeanStress:0;
    }
  }
}

void MkStress::CalcFI_1()
{
  FI_1 = FStress[0][0]+FStress[1][1]+FStress[2][2];
}

void MkStress::CalcFI_2()
{
  FI_2 = FStress[0][0]*FStress[1][1]+FStress[1][1]*FStress[2][2]+FStress[2][2]*FStress[0][0]
        -FStress[0][1]*FStress[0][1]+FStress[1][2]*FStress[1][2]+FStress[2][0]*FStress[2][0];
}

void MkStress::CalcFI_3()
{
  FI_3 = FStress[0][0]*FStress[1][1]*FStress[2][2]
        -FStress[0][0]*FStress[1][2]+FStress[1][2]
        -FStress[1][1]*FStress[2][0]+FStress[2][0]
        -FStress[2][2]*FStress[0][1]+FStress[0][1]
        +FStress[0][1]*FStress[1][2]+FStress[2][0]
        +FStress[0][1]*FStress[1][2]+FStress[2][0];
}

void MkStress::CalcFJ_2()
{
  FJ_2 = FI_1*FI_1/3.0-FI_2;
}

void MkStress::CalcInvariant()
{
  FI_1 = FStress[0][0]+FStress[1][1]+FStress[2][2];
  FI_2 = FStress[0][0]*FStress[1][1]+
         FStress[1][1]*FStress[2][2]+
         FStress[2][2]*FStress[0][0]-
         FStress[0][1]*FStress[0][1]+
         FStress[1][2]*FStress[1][2]+
         FStress[2][0]*FStress[2][0];
  FI_3 = FStress[0][0]*FStress[1][1]*FStress[2][2]
        -FStress[0][0]*FStress[1][2]+FStress[1][2]
        -FStress[1][1]*FStress[2][0]+FStress[2][0]
        -FStress[2][2]*FStress[0][1]+FStress[0][1]
        +FStress[0][1]*FStress[1][2]+FStress[2][0]
        +FStress[0][1]*FStress[1][2]+FStress[2][0];
  FJ_2 = FI_1*FI_1/3.0-FI_2;
  FMeanStress = FI_1/3.0;
//MkDebug("MkStress::CalcInvariant FI_1 %12.6f, FI_2 %12.6f, FI_3 %12.6f,FJ_2 %12.6f\n",FI_1, FI_2, FI_3,FJ_2);
}


void MkStress::CalcPhi()
{
  double a,b,q;

  CalcFI_1();
  CalcFI_2();
  CalcFI_3();

  a = (2*FI_1*FI_1*FI_1 - 9*FI_1*FI_2+27*FI_3)/54;
  q = fabs(FI_1*FI_1-3*FI_2)/9.0;
  b = sqrt(q*q*q);
  if (fabs(a) < 1.0e-3) {FPhi = -3.14159/2.0;}
  else if(fabs(a) > 1.0e-3 && fabs(b) < 1.0e-3) {
    FPhi = 0;
  }
  else {
    FPhi = acos(a/b);
  }
//MkDebug("MkStress::CalcPhi a %f b %f FPhi %f\n",a,b,FPhi);
}

void MkStress::CalcPrincipal()
{
  double q = fabs(FI_1*FI_1-3*FI_2)/9.0;
  CalcPhi();
  FPrincipal[0] = FI_1/3.0 + 2*sqrt(q)*cos(FPhi/3);
  FPrincipal[1] = FI_1/3.0 + 2*sqrt(q)*cos(FPhi/3+2*3.141592654/3.0);
  FPrincipal[2] = FI_1/3.0 + 2*sqrt(q)*cos(FPhi/3+4*3.141592654/3.0);
//MkDebug("MkStress::CalcPrincipal FPhi %f FPrincipal x %12.6f, FPrincipal y %12.6f, FPrincipal z %12.6f\n",FPhi, FPrincipal[0], FPrincipal[1], FPrincipal[2]);
}

void MkStress::CalcLMN()
{
  double theta,phi;
  double sum = fabs(FPrincipal[0])+fabs(FPrincipal[1])+fabs(FPrincipal[2]);
  double prc[3],str[3][3];

  if(sum < 1.e-3) {
    theta =0;
    phi = 0;
  }
  else {
    prc[0] = FPrincipal[0]/sum;
    prc[1] = FPrincipal[1]/sum;
    prc[2] = FPrincipal[2]/sum;

    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
	str[i][j] = FStress[i][j]/sum;
      }
    }

    if(fabs(prc[0]-prc[1])<1.e-3) {
      phi = 0;
    }
    else if(fabs(prc[0]-prc[1])>1.e-3) {
      phi = asin(sqrt(fabs((str[1][1]-prc[1])/(prc[0]-prc[1]))));
    }
    
    if (fabs(phi) < 1.e-3 && fabs(prc[2]-prc[0]) < 1.e-3 ) {
      theta = 0;
    } 
    else if (fabs(phi) < 1.e-3 && fabs(prc[2]-prc[0]) > 1.e-3 ) {
      theta = asin(sqrt(fabs((str[0][0]-prc[0])/(prc[2]-prc[0]))));
    } 
    else if (fabs(phi) >= 1.e-3) {
      theta = acos(str[0][1]/(prc[1]-prc[0])/sin(phi)/cos(phi));
    } 

    /*
    if(fabs(FPrincipal[0]-FPrincipal[1])<1.e-3) {
      phi = 0;
    }
    else if(fabs(FPrincipal[0]-FPrincipal[1])>1.e-3) {
      phi = asin(sqrt(fabs((FStress[1][1]-FPrincipal[1])/(FPrincipal[0]-FPrincipal[1]))));
    }
    
    if (fabs(phi) < 1.e-3 && fabs(FPrincipal[2]-FPrincipal[0]) < 1.e-3 ) {
      theta = 0;
    } 
    else if (fabs(phi) < 1.e-3 && fabs(FPrincipal[2]-FPrincipal[0]) > 1.e-3 ) {
      theta = asin(sqrt(fabs((FStress[0][0]-FPrincipal[0])/(FPrincipal[2]-FPrincipal[0]))));
    } 
    else if (fabs(phi) >= 1.e-3) {
      theta = acos(FStress[0][1]/(FPrincipal[1]-FPrincipal[0])/sin(phi)/cos(phi));
    } 
*/
  }

  LMN[0][0] = cos(theta)*cos(phi);  //l_1
  LMN[0][1] = cos(theta)*sin(phi);  //m_1
  LMN[0][2] = -sin(theta);          //n_1
  LMN[1][0] = -sin(phi);            //l_2
  LMN[1][1] = cos(phi);             //m_2
  LMN[1][2] = 0;                    //n_2
  LMN[2][0] = sin(theta)*cos(phi);  //l_3
  LMN[2][1] = sin(theta)*sin(phi);  //m_3
  LMN[2][2] = cos(theta);           //n_3 
//MkDebug("MkStress::CalcLMN phi %12.6f theta %12.6f \n",phi,theta);
//MkDebug("MkStress::CalcLMN 11 %12.6f 12 %12.6f 13 %12.6f\n",LMN[0][0],LMN[0][1],LMN[0][2]);
//MkDebug("MkStress::CalcLMN 21 %12.6f 22 %12.6f 23 %12.6f\n",LMN[1][0],LMN[1][1],LMN[1][2]);
//MkDebug("MkStress::CalcLMN 31 %12.6f 32 %12.6f 33 %12.6f\n",LMN[2][0],LMN[2][1],LMN[2][2]);
}

void MkStress::Eigen() 
{
  int i,j,m,n,k,l;
  MkMatrix a(3,3), b(3,3), d(3,3), q(3,3), r(3,3);
  double sum, theta,c,s,MP = 3.141592654;

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      a(i,j) = FStress[i][j];
    }
  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      r(i,j) = 0.0;
    }
    r(i,i) = 1.0;
  }

  i = 0; j = 1;

  if (fabs(a(0,2)) < fabs(a(0,1))) j = 2;
  sum = 1;

  while (sum > 1.0E-10) {
    m = 3 - j;
    n = 1 - i;
    if (fabs(a(0,m)) > fabs(a(n,2))) {
      i = 0;
      j = m;
    }
    else {
      i = n;
      j = 2;
    }
    theta = 0;
    if (a(i,j) != 0) {
      theta = PI / 4;
    }

    if (a(i,i) != a(j,j)) {
      theta = 0.5 * atan(2*a(i,j)/(a(i,i)-a(j,j)));
    }
    c = cos(theta);
    s = sin(theta);

    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
	q(k,l) = 0.0;
      }
      q(k,k) = 1.0;
    }

    q(i,i) = c;
    q(i,j) = s;
    q(j,i) =-s;
    q(j,j) = c;

    d = q*a;
    b=q;b.Transpose();
    a = d*b;
    b = q*r;

    r = b;
    sum = fabs(a(0,1)) + fabs(a(0,2)) + fabs(a(1,2));
  }

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      LMN[i][j] = r(i,j);
    }
    FPrincipal[i] = a(i,i);
  }

//MkDebug("MkStress::Stress 11 %12.6f 12 %12.6f 13 %12.6f\n",FStress[0][0],FStress[0][1],FStress[0][2]);
//MkDebug("MkStress::Stress 21 %12.6f 22 %12.6f 23 %12.6f\n",FStress[1][0],FStress[1][1],FStress[1][2]);
//MkDebug("MkStress::Stress 31 %12.6f 32 %12.6f 33 %12.6f\n",FStress[2][0],FStress[2][1],FStress[2][2]);

//MkDebug("MkStress::CalcPrincipal FPrincipal x %12.6f, FPrincipal y %12.6f, FPrincipal z %12.6f\n", FPrincipal[0], FPrincipal[1], FPrincipal[2]);
//MkDebug("MkStress::CalcLMN 11 %12.6f 12 %12.6f 13 %12.6f\n",LMN[0][0],LMN[0][1],LMN[0][2]);
//MkDebug("MkStress::CalcLMN 21 %12.6f 22 %12.6f 23 %12.6f\n",LMN[1][0],LMN[1][1],LMN[1][2]);
//MkDebug("MkStress::CalcLMN 31 %12.6f 32 %12.6f 33 %12.6f\n",LMN[2][0],LMN[2][1],LMN[2][2]);

}

MkStress & MkStress::operator=(const MkStress &str)
{
  for (int i= 0;i<3;i++) {
    FPrincipal[i] = str.FPrincipal[i];
    for (int j=0;j<3;j++) {
      FStress[i][j] = str.FStress[i][j];
      FDeviatoric[i][j] = str.FDeviatoric[i][j];
      LMN[i][j] = str.LMN[i][j];
    }
  }
  FMeanStress = str.FMeanStress;

  FI_1 = str.FI_1;
  FI_2 = str.FI_2;
  FI_3 = str.FI_3;
  FPhi = str.FPhi;
  return *this;
}

void MkStress::Out(char *fname)
{
  FILE *fp;
  char str[256],s[256];

  fp = fopen(fname,"a");
  fprintf(fp,"Output of Stress ");
  sprintf(str,"A-th row [%10.3f %10.3f %10.3f %10.3f ]",0.0,1.0,2.0,3.0);
  fprintf(fp,str);
  for (int i=0;i<4;i++) {
    sprintf(str,"%d-th row [",i);
    for (int j=0;j<3;j++) {
      sprintf(s,"%10.3f ",FStress[i][j]);
      strcat(str,s);
    }
    sprintf(s,"]");
    strcat(str,s);
    fprintf(fp,str);
  }
  fclose(fp);
}

void MkStress::Out()
{
  char str[256],s[256];

  printf("Output of Stress \n");
  sprintf(str,"A-th row [%10.5f %10.5f %10.5f ]\n",0.0,1.0,2.0);
  printf(str);
  for (int i=0;i<3;i++) {
    sprintf(str,"%d-th row [",i);
    for (int j=0;j<3;j++) {
      sprintf(s,"%10.5f ",FStress[i][j]);
      strcat(str,s);
    }
    sprintf(s,"]\n");
    strcat(str,s);
    printf(str);
  }
}
