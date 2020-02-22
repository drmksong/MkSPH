#include "MkSmoothFunc.h"

MkSmoothFunc::MkSmoothFunc()
{
  SmoothFuncType = smtNone;
  Dimension = 0;
  SmoothLength = 0;
  Kappa = 0;
  AlphaD = 0;
}

MkSmoothFunc::MkSmoothFunc(MkSmoothFuncType smft, int dim, float smoothlen, float kappa)
{
  SmoothFuncType = smft;
  Dimension = dim;
  SmoothLength = smoothlen;
  Kappa = kappa;
  SetupAlphaD();
}

void MkSmoothFunc::SetupAlphaD(void)
{
  if (fabs(SmoothLength) < _MIN_SMOOTH_LENGTH_)
  {
    AlphaD = 0;
    return;
  }
  switch (SmoothFuncType)
  {
  case smtGaussian:
    switch (Dimension)
    {
    case 1:
      AlphaD = 1 / (sqrt(3.14159) * SmoothLength);
      break;
    case 2:
      AlphaD = 1 / (3.14159 * pow(SmoothLength, 2));
      break;
    case 3:
      AlphaD = 1 / (pow(3.14159, 1.5) * pow(SmoothLength, 3));
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

float MkSmoothFunc::W(float dist)
{
  static float value;
  float R = dist / SmoothLength;
  switch (SmoothFuncType)
  {
  case smtGaussian:
    value = AlphaD * exp(-pow(R, 2));
    break;
  default:
    value = 0;
    break;
  }
  return value;
}

float MkSmoothFunc::dWdR(float dist)
{
  static float value;
  float R = dist / SmoothLength;
  switch (SmoothFuncType)
  {
  case smtGaussian:
    value = -2 * AlphaD * R * exp(-pow(R, 2));
    break;
  default:
    value = 0;
    break;
  }
  return value;
}

float MkSmoothFunc::d2WdR2(float dist)
{
  static float value;
  float R = dist / SmoothLength;
  switch (SmoothFuncType)
  {
  case smtGaussian:
    value = -2 * AlphaD * exp(-pow(R, 2)) + 4 * AlphaD * pow(R, 2) * exp(-pow(R, 2));
    break;
  default:
    value = 0;
    break;
  }
  return value;
}

float MkSmoothFunc::dRdX(float dx, float dist)
{
  float R = dist / SmoothLength;
  return dx / R / SmoothLength;
}
