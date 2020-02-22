#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "MkPoint.hpp"
#include "MkVect.hpp"
#include "MkColor.hpp"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

extern double Zero;

class MkParticle : public MkPoint
{
private:
  double R, G, B;
  double Radius;
  double Density, Density1; //??
  double Mass;
  double Viscosity;
  double Temp, Temp1; //??
  double Pressure, Pressure1;
  double XForceInt, YForceInt;
  double XForceVis, YForceVis;
  double XForceArtVis, YForceArtVis;
  double XForceExt, YForceExt;
  double XAccel[3], YAccel[3];       // 1 : t-dt, 2 : t, 3 : t+dt
  double XVelocity[3], YVelocity[3]; // 1 : t-dt, 2 : t, 3 : t+dt
  double XPosition[3], YPosition[3]; // 1 : t-dt, 2 : t, 3 : t+dt
  int GridNum;
  bool isFixed;

public:
  double XX, YY, ZZ;
  MkParticle(void) : MkPoint() { Initialize(); }
  MkParticle(double x, double y) : MkPoint(x, y) { Initialize(); }
  MkParticle(double x, double y, double z) : MkPoint(x, y, z) { Initialize(); }
  ~MkParticle(void)
  {
  }
  void Initialize();
  void SetR(double r) { R = r; }
  void SetG(double g) { G = g; }
  void SetB(double b) { B = b; }
  void SetDensity(double d) { Density = d; }
  void SetDensity1(double d) { Density1 = d; }
  void SetTemp(double t) { Temp = t; }
  void SetTemp1(double t) { Temp1 = t; }
  void SetPressure(double p) { Pressure = p; }
  void SetPressure1(double p) { Pressure1 = p; }
  void SetGridNum(int num) { GridNum = num; }
  void SetRadius(double r) { Radius = r; }
  void SetFixity(bool flag) { isFixed = flag; }
  void SetMass(double mass) { Mass = mass; }
  void SetViscosity(double v) { Viscosity = v; }
  void SetXForceInt(double xf) { XForceInt = xf; }
  void SetYForceInt(double yf) { YForceInt = yf; }
  void SetXForceVis(double xf) { XForceVis = xf; }
  void SetYForceVis(double yf) { YForceVis = yf; }
  void SetXForceArtVis(double xf) { XForceArtVis = xf; }
  void SetYForceArtVis(double yf) { YForceArtVis = yf; }
  void SetXForceExt(double xf) { XForceExt = xf; }
  void SetYForceExt(double yf) { YForceExt = yf; }
  void SetXVelocity(int i, double xv)
  {
    if (i >= 0 && i <= 2)
      XVelocity[i] = xv;
  }
  void SetYVelocity(int i, double yv)
  {
    if (i >= 0 && i <= 2)
      YVelocity[i] = yv;
  }
  void SetXAccel(int i, double xa)
  {
    if (i >= 0 && i <= 2)
      XAccel[i] = xa;
  }
  void SetYAccel(int i, double ya)
  {
    if (i >= 0 && i <= 2)
      YAccel[i] = ya;
  }
  void SetXPosition(int i, double xp)
  {
    if (i >= 0 && i <= 2)
      XPosition[i] = xp;
  }
  void SetYPosition(int i, double yp)
  {
    if (i >= 0 && i <= 2)
      YPosition[i] = yp;
  }

  double GetR(void) { return R; }
  double GetG(void) { return G; }
  double GetB(void) { return B; }
  double GetDensity(void) { return Density; }
  double GetDensity1(void) { return Density1; }
  double GetTemp(void) { return Temp; }
  double GetTemp1(void) { return Temp1; }
  double GetPressure(void) { return Pressure; }
  double GetPressure1(void) { return Pressure1; }
  int GetGridNum(void) { return GridNum; }
  double GetRadius(void) { return Radius; }
  bool GetFixity(void) { return isFixed; }
  double GetMass(void) { return Mass; }
  double GetViscosity(void) { return Viscosity; }
  double GetXForceInt(void) { return XForceInt; }
  double GetYForceInt(void) { return YForceInt; }
  double GetXForceVis(void) { return XForceVis; }
  double GetYForceVis(void) { return YForceVis; }
  double GetXForceArtVis(void) { return XForceArtVis; }
  double GetYForceArtVis(void) { return YForceArtVis; }
  double GetXForceExt(void) { return XForceExt; }
  double GetYForceExt(void) { return YForceExt; }
  double &GetXVelocity(int i)
  {
    if (i >= 0 && i <= 2)
      return XVelocity[i];
    else
      return Zero;
  }
  double &GetYVelocity(int i)
  {
    if (i >= 0 && i <= 2)
      return YVelocity[i];
    else
      return Zero;
  }
  double &GetXAccel(int i)
  {
    if (i >= 0 && i <= 2)
      return XAccel[i];
    else
      return Zero;
  }
  double &GetYAccel(int i)
  {
    if (i >= 0 && i <= 2)
      return YAccel[i];
    else
      return Zero;
  }
  double &GetXPosition(int i)
  {
    if (i >= 0 && i <= 2)
      return XPosition[i];
    else
      return Zero;
  }
  double &GetYPosition(int i)
  {
    if (i >= 0 && i <= 2)
      return YPosition[i];
    else
      return Zero;
  }

  double RefreshTemp(void)
  {
    double d = fabs(Temp - Temp1);
    Temp = Temp1;
    return d;
  }
  double RefreshPos(void)
  {
    double dx, dy, dz;
    dx = fabs(XX - X);
    dy = fabs(YY - Y);
    dz = fabs(ZZ - Z);
    X = XX;
    Y = YY;
    Z = ZZ;
    return dx + dy + dz;
  }
  void SetupTempColor();
  friend MkParticle &operator*(MkParticle &rp, MkMatrix4<double> &rm);
  friend MkParticle &operator*(MkParticle &rp, double f);
  int Solve(void);
#if defined(__GL_H__)
  void Draw();
#endif
};

class MkParticles
{
protected:
  MkParticle *FParticle;
  MkParticle FCenter;
  int FSize; //Actual size of points
  int FSizeOfArray;
  void FindCenter();
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  MkColor Color;
  double DotRadius;
#endif
public:
  MkParticles(int size, MkParticle *rps);
  MkParticles(int size);
  MkParticles()
  {
    FSizeOfArray = FSize = 0;
    FParticle = NULL;
  }
  ~MkParticles();
  virtual void Initialize(int size);
  virtual void Initialize(int size, MkParticle *);
  int GetSize() { return FSize; };
  int GetNumber() { return FSize; };
  MkParticle *GetParticles() { return FParticle; }
  MkParticle GetCenter()
  {
    FindCenter();
    return FCenter;
  };
  bool Add(MkParticle point); // change of size of point
  bool Add(int index, MkParticle point);
  bool Add(MkParticles &p)
  {
    for (int i = 0; i < p.GetSize(); i++)
      Add(p[i]);
    return true;
  }
  bool Delete(MkParticle point); // change of size of point
  bool Delete(int index);
  int Grow(int Delta);   // change of size of array
  int Shrink(int Delta); // change of size of array
  bool Swap(int i, int j);
  bool Clear();
#ifdef __BCPLUSPLUS__
  TColor GetColor()
  {
    return Color;
  };
  void SetColor(TColor c) { Color = c; }
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void SetColor(MkColor c)
  {
    Color = c;
  }
  void SetDotRadius(double r) { DotRadius = r; }
  MkColor GetColor() { return Color; }
  double GetDotRadius() { return DotRadius; }
#endif
  bool hasParticle(MkParticle &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FParticle[i] == pnt)
        return true;
    return false;
  }
  int numParticle(MkParticle &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FParticle[i] == pnt)
        return i;
    return -1;
  }
  virtual MkParticle &operator[](int);
  MkParticles &operator*=(MkMatrix4<double> &rm);
  friend MkParticles &operator*(MkParticles &rps, MkMatrix4<double> &rm);

  MkParticles &Translate(MkParticle rp);
  MkParticles &Translate(double x, double y, double z);
  MkParticles &Rotate(double alpha, double beta, double gamma);
  MkParticles &RotateInX(double ang);
  MkParticles &RotateInY(double ang);
  MkParticles &RotateInZ(double ang);
  MkParticles &RotateInA(double ang, double l, double m, double n);
  MkParticles &Scale(double sx, double sy, double sz);

  MkParticles &operator=(MkParticles &points);
  bool operator==(MkParticles &points);
  bool operator!=(MkParticles &points) { return !operator==(points); }

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif

  //#if defined(__GL_H__) && defined(_WINDOWS_)
#if defined(__GL_H__)
  void Draw();
#endif
};

extern MkParticle NullParticle;
extern MkParticles NullParticles;

#endif
