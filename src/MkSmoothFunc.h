#ifndef MkSmoothFuncH
#define MkSmoothFuncH
#include <math.h>

#define _MIN_SMOOTH_LENGTH_ 1e-3

typedef enum
{
  smtNone,
  smtQuartic,
  smtPieceCubicSpline,
  smtPiceQuartic,
  smtPieceQuintic,
  smtQuadratic,
  smtSuperGaussian,
  smtGaussian,
  smtDomeShape,
  smtNewQuartic
} MkSmoothFuncType;
class MkSmoothFunc
{
private:
protected:
  MkSmoothFuncType SmoothFuncType;
  int Dimension;
  float SmoothLength;
  float Kappa;
  float AlphaD;

public:
  MkSmoothFunc();
  MkSmoothFunc(MkSmoothFuncType smft, int dim, float smoothlen, float kappa);
  void SetupSmoothFunc(MkSmoothFuncType smft, int dim, float smoothlen, float kappa)
  {
    SetSmoothFuncType(smft);
    SetDim(dim);
    SetSmoothLen(smoothlen);
    SetKappa(kappa);
    SetupAlphaD();
  }
  void SetSmoothFuncType(MkSmoothFuncType smft) { SmoothFuncType = smft; }
  void SetDim(int dim)
  {
    Dimension = dim;
    SetupAlphaD();
  }
  void SetSmoothLen(float sml)
  {
    SmoothLength = sml;
    SetupAlphaD();
  }
  void SetKappa(float kappa)
  {
    Kappa = kappa;
    SetupAlphaD();
  }
  void SetupAlphaD(void);
  MkSmoothFuncType GetSmoothFuncType(void) { return SmoothFuncType; }
  int GetDim(void) { return Dimension; }
  float GetSmoothLen(void) { return SmoothLength; }
  float GetKappa(void) { return Kappa; }
  float W(float dist);
  float dWdR(float dist);
  float d2WdR2(float dist);
  float dRdX(float dx, float dist);
};

#endif
