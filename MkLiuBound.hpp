#ifndef LiuBOUNDARY_HPP
#define LiuBOUNDARY_HPP

#include <boost/shared_array.hpp>
#include "MkPoint.hpp"
#include "MkPlane.hpp"
#include "MkArray.hpp"
#include "MkColor.hpp"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

extern double Zero;

enum BType
{
  btIdealGas = 1,
  btWater
};

class MkLiuBoundary : public MkPlane
{
public:
  //#ifdef __GL_H__
  double R, G, B;
  double Transparancy; // 0 to 1
                       //#endif
  double Radius;
  CType ColorType;

  //double X,Y,Z from MkPlane checked
  double XVel, YVel, ZVel;    // VX
  double XAVel, YAVel, ZAVel; // AveVel
  //  BType BoundaryType; // IType Boundary Type
  int BoundaryType; // IType Boundary Type

  double Mass;       // real mass checked
  double Rho;        // density checked,
  double Volume;     // volume to calculate the selfdens (density when no neighbor boundarys)  mk
  double Eta;        // dynamic viscosity
  double Press;      // P  checked
  double Temp;       // T  checked
  double Energy;     //U  checked
  double SoundSpeed; //C
  double SmoothLen;  //Hsml checked

  int CountIac;

  double DUDt; //**
  double INDUDt;
  double AVDUDt;
  double AHDUDt;
  double ADDUDt;
  double APDUDt;
  double DRhoDt;
  double DVXDt; //**
  double INDVXDt;
  double ARDVXDt;
  double ADDVXDt; // artificial drag
  double APDVXDt;
  double EXDVXDt;
  double DVYDt; //**
  double INDVYDt;
  double ARDVYDt;
  double ADDVYDt; // artificial drag
  double APDVYDt;
  double EXDVYDt;
  double DVZDt; //**
  double INDVZDt;
  double ARDVZDt;
  double ADDVZDt; // artificial drag
  double APDVZDt;
  double EXDVZDt;
  double TDSDt;

  double W;    // kernal of all interaction pairs
  double DWDX; //derivative of kernal with respect to x
  double DWDY; //derivative of kernal with respect to y
  double DWDZ; //derivative of kernal with respect to z

public:
  MkLiuBoundary(void) : MkPlane() { Initialize(); }
  MkLiuBoundary(int i) : MkPlane(i) { Initialize(); }
  MkLiuBoundary(MkPoint rp1, MkPoint rp2, MkPoint rp3) : MkPlane(rp1, rp2, rp3)
  {
    Initialize();
  }
  ~MkLiuBoundary(void) {}

  void Initialize();

  void SetupColor();
  void SetColorType(CType ct) { ColorType = ct; }
  void SetR(double r) { R = r; }
  void SetG(double g) { G = g; }
  void SetB(double b) { B = b; }
  void SetTransparancy(double t) { Transparancy = t; }
  void SetRadius(double r) { Radius = r; }
  void SetXVel(double xv) { XVel = xv; }
  void SetYVel(double yv) { YVel = yv; }
  void SetZVel(double yv) { YVel = yv; }
  void SetXAVel(double xv) { XAVel = xv; }
  void SetYAVel(double yv) { YAVel = yv; }
  void SetZAVel(double yv) { YAVel = yv; }
  //void SetBoundaryType(BType pt){BoundaryType = (int)pt;}
  void SetBoundaryType(int pt) { BoundaryType = pt; }
  void SetMass(double mass) { Mass = mass; }
  void SetRho(double rho) { Rho = rho; }
  void SetEta(double eta) { Eta = eta; }
  void SetPress(double press) { Press = press; }
  void SetTemp(double temp) { Temp = temp; }
  void SetEnergy(double u) { Energy = u; }
  void SetSoundSpeed(double c) { SoundSpeed = c; }
  void SetSmoothLen(double h) { SmoothLen = h; }
  void SetCountIac(double iac) { CountIac = iac; }
  void SetDRhoDt(double d) { DRhoDt = d; }
  void SetDVXDt(double d) { DVXDt = d; }
  void SetTDSDt(double d) { TDSDt = d; }
  void SetDUDt(double d) { DUDt = d; }
  void SetW(double w) { W = w; }
  void SetDWDX(double d) { DWDX = d; }

  double GetR() { return R; }
  double GetG() { return G; }
  double GetB() { return B; }
  double GetTransparancy() { return Transparancy; }
  double GetRadius() { return Radius; }
  double &GetXVel(int i) { return XVel; }
  double &GetYVel(int i) { return YVel; }
  double &GetZVel(int i) { return YVel; }
  double &GetXAVel(int i) { return XAVel; }
  double &GetYAVel(int i) { return YAVel; }
  double &GetZAVel(int i) { return YAVel; }
  int GetBoundaryType() { return BoundaryType; }
  double GetMass() { return Mass; }
  double GetRho() { return Rho; }
  double GetEta() { return Eta; }
  double GetPress() { return Press; }
  double GetTemp() { return Temp; }
  double GetEnergy() { return Energy; }
  double GetSoundSpeed() { return SoundSpeed; }
  double GetSmoothLen() { return SmoothLen; }
  double GetDRhoDt() { return DRhoDt; }
  double GetDVXDt() { return DVXDt; }
  double GetTDSDt() { return TDSDt; }
  double GetDUDt() { return DUDt; }
  double GetW() { return W; }
  double GetDWDX() { return DWDX; }

  MkLiuBoundary &operator=(const MkLiuBoundary &lb);

#if defined(__GL_H__)
  void Draw();
#endif
};

class MkLiuBoundarys
{
protected:
  boost::shared_array<MkLiuBoundary> FBoundary;
  MkPoint FCenter;
  int FSize; //Actual size of planes
  void FindCenter();
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  MkColor Color;
  double DotRadius;
#endif
public:
  MkLiuBoundarys(int size, MkLiuBoundary *rps);
  MkLiuBoundarys(int size);
  MkLiuBoundarys()
  {
    FSize = 0;
    FBoundary = NULL;
  }
  ~MkLiuBoundarys();
  void Initialize(int size);
  void Initialize(int size, MkLiuBoundary *);
  int GetSize() { return FSize; };
  int GetNumber() { return FSize; };
  boost::shared_array<MkLiuBoundary> GetBoundarys() { return FBoundary; }
  MkPoint GetCenter()
  {
    FindCenter();
    return FCenter;
  };
  bool Add(MkLiuBoundary plane); // change of size of plane
  bool Add(int index, MkLiuBoundary plane);
  bool Add(MkLiuBoundarys &p)
  {
    for (int i = 0; i < p.GetSize(); i++)
      Add(p[i]);
    return true;
  }
  bool Delete(MkLiuBoundary plane); // change of size of plane
  bool Delete(int index);
  int Grow(int Delta);   // change of size of array
  int Shrink(int Delta); // change of size of array
  bool Swap(int i, int j);
  bool Clear();
  void SetColorType(CType ct)
  {
    for (int i = 0; i < FSize; i++)
      FBoundary[i].SetColorType(ct);
  }
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
  bool hasBoundary(MkLiuBoundary &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FBoundary[i] == pnt)
        return true;
    return false;
  }
  int numBoundary(MkLiuBoundary &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FBoundary[i] == pnt)
        return i;
    return -1;
  }
  virtual MkLiuBoundary &operator[](int);
  MkLiuBoundarys &operator*=(MkMatrix4<double> &rm);
  friend MkLiuBoundarys &operator*(MkLiuBoundarys &rps, MkMatrix4<double> &rm);

  void Translate(MkPoint rp);
  void Translate(double x, double y, double z);
  void Rotate(double alpha, double beta, double gamma);
  void RotateInX(double ang);
  void RotateInY(double ang);
  void RotateInZ(double ang);
  void RotateInA(double ang, double l, double m, double n);
  void Scale(double sx, double sy, double sz);

  MkLiuBoundarys &operator=(MkLiuBoundarys &planes);
  bool operator==(MkLiuBoundarys &planes);
  bool operator!=(MkLiuBoundarys &planes) { return !operator==(planes); }

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

extern MkLiuBoundary NullBoundary;
extern MkLiuBoundarys NullBoundarys;

#endif
