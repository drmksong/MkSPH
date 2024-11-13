#ifndef MkLiuKernelH
#define MkLiuKernelH
#include <stdlib.h>
#include <stddef.h>
#include <string>
#include <math.h>
#include "MkMisc.hpp"

#define _MIN_SMOOTH_LENGTH_ 1e-6

typedef enum
{
  knlNone,
  knlCubicSpline,
  knlGaussian,
  knlQuintic
} MkKernelType;
class MkLiuKernel
{
private:
protected:
  MkKernelType KernelType;
  int Dim;
  double SmoothLen;
  double AlphaD;

public:
  MkLiuKernel();
  MkLiuKernel(MkKernelType knl, int dim, double smoothlen);
  void Clear()
  {
    KernelType = knlNone;
    Dim = 0;
    AlphaD = 0;
  };
  void SetupLiuKernel(MkKernelType knl, int dim, double smoothlen)
  {
    SetKernelType(knl);
    SetDim(dim);
    SetSmoothLen(smoothlen);
    SetupAlphaD();
  }
  void SetKernelType(MkKernelType knl) { KernelType = knl; }
  void SetDim(int dim)
  {
    Dim = dim;
    SetupAlphaD();
  }
  void SetSmoothLen(double sml)
  {
    SmoothLen = sml;
    SetupAlphaD();
  }
  void SetupAlphaD(void);

  MkKernelType GetKernelType(void) { return KernelType; }
  int GetDim(void) { return Dim; }
  double GetSmoothLen(void) { return SmoothLen; }
  double GetAlphaD() { return AlphaD; }

  double W(double dist);
  double dWdR(double dist);
  double dWdX(double dist, double dx, double dy, double dz);
  double dWdY(double dist, double dx, double dy, double dz);
  double dWdZ(double dist, double dx, double dy, double dz);

  class Alloc
  {
  public:
    std::string What;
    Alloc(std::string what) : What(what) {}
    std::string what() { return What; }
  };
  class Size
  {
  public:
    std::string What;
    int N;
    Size(std::string what, int n) : What(what), N(n) {}
    std::string what() { return What; }
  };
  class Range
  {
  public:
    std::string What;
    int N;
    Range(std::string what, int n) : What(what), N(n) {}
    std::string what() { return What; }
  };
};
#endif
