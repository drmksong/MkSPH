#include "MkLiuKernel.hpp"

MkLiuKernel::MkLiuKernel()
{
  KernelType = knlNone;
  Dim=0;
  SmoothLen=0;
  AlphaD=0;
}

MkLiuKernel::MkLiuKernel(MkKernelType knl, int dim, double smoothlen)
{
  KernelType = knl;
  Dim=dim;
  SmoothLen=smoothlen;
  SetupAlphaD();
}

void MkLiuKernel::SetupAlphaD(void)
{
  if (fabs(SmoothLen) < _MIN_SMOOTH_LENGTH_) {AlphaD = 0; return;}  
  switch(KernelType) {
  case knlCubicSpline:
    switch (Dim) {
    case 1:
      AlphaD = 1/SmoothLen ;
      break;
    case 2:
      AlphaD = 15.0/(7.0*3.14159*SmoothLen*SmoothLen);
      break;
    case 3:
      AlphaD = 3.0/(2.0*3.14159*SmoothLen*SmoothLen*SmoothLen);
      break;
    default:
      break;
    }
    break;

  case knlGaussian:
    switch (Dim) {
    case 1:
      AlphaD = 1/(sqrt(3.14159)*SmoothLen);
      break;
    case 2:
      AlphaD = 1/(3.14159*pow(SmoothLen,2));
      break;
    case 3:
      AlphaD = 1/(pow(3.14159,1.5)*pow(SmoothLen,3));
      break;
    default:
      break;
    }
    break;

  case knlQuintic:
    switch (Dim) {
    case 1:
      AlphaD = 1.0 /(120.0*SmoothLen) ;
      break;
    case 2:
      AlphaD = 7.0 / (478.0*3.14159*pow(SmoothLen,2)) ;
      break;
    case 3:
      AlphaD = 1.0 / (120.0*3.14159*pow(SmoothLen,3));
      break;
    default:
      break;
    }
    break;

  default:
    AlphaD = 0;
    break;
  }
}

double MkLiuKernel::W(double dist)
{
  static double value; 
  double R = dist/SmoothLen;
  switch (KernelType) {
  case knlCubicSpline:
    if (0<=R && R<= 1.0) 
      value = AlphaD * (2.0/3.0 - R*R + R*R*R / 2.);
    else if (1<R && R <=2.0) 
      value = AlphaD * 1.0/6.0 * (2.0-R)*(2.0-R)*(2.0-R);
    else value = 0;
    break;

  case knlGaussian:
    if (0<=R && R<=3.0) 
      value = AlphaD*exp(-pow(R,2));
    else value = 0;
    break;

  case knlQuintic:
    if (0<=R && R<= 1.0) 
      value = AlphaD * ( pow(3-R,5) - 6*pow(2-R,5) + 15*pow(1-R,5) );
    else if (1<=R && R<= 2.0) 
      value = AlphaD * ( pow(3-R,5) - 6*pow(2-R,5) );
    else if (2<=R && R<= 3.0) 
      value = AlphaD * ( pow(3-R,5));
    else value = 0;
    break;

  default:
    value=0;
    break;
  }
  return value;
}

double MkLiuKernel::dWdR(double dist)
{
  static double value; 
  double R = dist/SmoothLen;
  switch (KernelType) {
  case knlCubicSpline:
    if (0<=R && R<= 1.0) 
      value = AlphaD * ( -2.0* R + 3.0*R*R /2.0)/SmoothLen;
    else if (1<R && R <=2.0) 
      value = AlphaD * (-3.0/6.0) * (2.0-R)*(2.0-R)/SmoothLen;
    else value = 0;
    break;
  case knlGaussian:
    if (0<=R && R <= 3) 
      value = -2*AlphaD*R*exp(-pow(R,2))/SmoothLen;
    else 
      value = 0;
    break;
  case knlQuintic:
    if (0<=R && R<= 1.0) 
      value = AlphaD * ( pow(3-R,5) - 6*pow(2-R,5) + 15*pow(1-R,5) )/SmoothLen;
    else if (1<=R && R<= 2.0) 
      value = AlphaD * ( pow(3-R,5) - 6*pow(2-R,5) )/SmoothLen;
    else if (2<=R && R<= 3.0) 
      value = AlphaD * ( pow(3-R,5))/SmoothLen;
    else value = 0;
    break;
  default:
    value=0;
    break;
  }
  return value;
}

double MkLiuKernel::dWdX(double dist, double dx, double dy, double dz)
{
  static double value;
  
  double dist2 = sqrt(dx*dx+dy*dy+dz*dz);
  if (fabs(fabs(dist)-dist2)>1.0e-3) {MkDebug("dist %f is not tally with length of dx (%f, %f, %f) vector \n",dist,dx,dy,dz );return 0;}
  double R = dist/SmoothLen;
  
  switch(KernelType) {
  case knlCubicSpline:
    if (0<=R && R<= 1.0) 
      value = AlphaD * (-2.+3./2.*R)/SmoothLen/SmoothLen * dx;
    else if (1<R && R <=2.0) 
      value =-AlphaD * 1.0/6.0 * 3.0*(2.0-R)*(2.0-R)/SmoothLen * (dx/dist);
    else value = 0;
    break;

  case knlGaussian:
    if (0<=R && R<=3.0) 
      value = AlphaD*exp(-pow(R,2))* ( -2.* dx/SmoothLen/SmoothLen);
    else value = 0;
    break;

  case knlQuintic:
    if (0<=R && R<= 1.0) 
      value = AlphaD * ( (-120 + 120*R - 50*R*R)/ SmoothLen/ SmoothLen * dx );
    else if (1<=R && R<= 2.0) 
      value = AlphaD * (-5*pow((3-R),4) + 30*pow(2-R,4))/ SmoothLen * (dx/dist);
    else if (2<=R && R<= 3.0) 
      value = AlphaD * (-5*pow((3-R),4))/ SmoothLen * (dx/dist);
    else value = 0;
    break;
  default:
    value = 0;
    break;
  }
  return value;
}


double MkLiuKernel::dWdY(double dist, double dx, double dy, double dz)
{
  static double value;
  
  double dist2 = sqrt(dx*dx+dy*dy+dz*dz);
  if (fabs(fabs(dist)-dist2)>1.0e-3) {MkDebug("dist %f is not tally with length of dx (%f, %f, %f) vector \n",dist,dx,dy,dz );return 0;}
  double R = dist/SmoothLen;
  
  switch(KernelType) {
  case knlCubicSpline:
    if (0<=R && R<= 1.0) 
      value = AlphaD * (-2.+3./2.*R)/SmoothLen/SmoothLen * dy;
    else if (1<R && R <=2.0) 
      value =-AlphaD * 1.0/6.0 * 3.0*(2.0-R)*(2.0-R)/SmoothLen * (dy/dist);
    else value = 0;
    break;

  case knlGaussian:
    if (0<=R && R<=3.0) 
      value = AlphaD*exp(-pow(R,2))* ( -2.* dy/SmoothLen/SmoothLen);
    else value = 0;
    break;

  case knlQuintic:
    if (0<=R && R<= 1.0) 
      value = AlphaD * ( (-120 + 120*R - 50*R*R)/ SmoothLen/ SmoothLen * dy );
    else if (1<=R && R<= 2.0) 
      value = AlphaD * (-5*pow((3-R),4) + 30*pow(2-R,4))/ SmoothLen * (dy/dist);
    else if (2<=R && R<= 3.0) 
      value = AlphaD * (-5*pow((3-R),4))/ SmoothLen * (dy/dist);
    else value = 0;
    break;
  default:
    value = 0;
    break;
  }
  return value;
}


double MkLiuKernel::dWdZ(double dist, double dx, double dy, double dz)
{
  static double value;
  
  double dist2 = sqrt(dx*dx+dy*dy+dz*dz);
  if (fabs(fabs(dist)-dist2)>1.0e-3) {MkDebug("dist %f is not tally with length of dx (%f, %f, %f) vector \n",dist,dx,dy,dz );return 0;}
  double R = dist/SmoothLen;
  
  switch(KernelType) {
  case knlCubicSpline:
    if (0<=R && R<= 1.0) 
      value = AlphaD * (-2.+3./2.*R)/SmoothLen/SmoothLen * dz;
    else if (1<R && R <=2.0) 
      value =-AlphaD * 1.0/6.0 * 3.0*(2.0-R)*(2.0-R)/SmoothLen * (dz/dist);
    else value = 0;
    break;

  case knlGaussian:
    if (0<=R && R<=3.0) 
      value = AlphaD*exp(-pow(R,2))* ( -2.* dz/SmoothLen/SmoothLen);
    else value = 0;
    break;

  case knlQuintic:
    if (0<=R && R<= 1.0) 
      value = AlphaD * ( (-120 + 120*R - 50*R*R)/ SmoothLen/ SmoothLen * dz );
    else if (1<=R && R<= 2.0) 
      value = AlphaD * (-5*pow((3-R),4) + 30*pow(2-R,4))/ SmoothLen * (dz/dist);
    else if (2<=R && R<= 3.0) 
      value = AlphaD * (-5*pow((3-R),4))/ SmoothLen * (dz/dist);
    else value = 0;
    break;
  default:
    value = 0;
    break;
  }
  return value;
}


