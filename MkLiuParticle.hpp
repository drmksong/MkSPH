#ifndef LiuPARTICLE_HPP
#define LiuPARTICLE_HPP

#include <boost/shared_array.hpp>
#include "MkPoint.hpp"
#include "MkArray.hpp"
#include "MkColor.hpp"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

extern double Zero;

enum PType
{
  ptIdealGas = 1,
  ptWater
};

class MkLiuParticle : public MkPoint
{
public:
  double R, G, B;
  double Radius;
  CType ColorType;

  //double X,Y,Z from MkPoint checked
  double XVel, YVel, ZVel;    // VX
  double XAVel, YAVel, ZAVel; // AveVel
  //  PType ParticleType; // IType Particle Type
  int ParticleType; // IType Particle Type
  bool Virtual;

  double Mass;       // real mass checked
  double Rho;        // density checked,
  double Rho_Norm;   // norm density for Norm_Density()
  double Volume;     // volume to calculate the selfdens (density when no neighbor particles)  mk
  double Eta;        // dynamic viscosity
  double Press;      // P  checked
  double Temp;       // T  checked
  double Energy;     //U  checked
  double SoundSpeed; //C
  double SmoothLen;  //Hsml checked

  int CountIac;
  int GI, GJ, GK; // Grid Index for MkLiuGrids

  double DUDt; //**
  double INDUDt;
  double AVDUDt;
  double AHDUDt;
  double ADDUDt; // artificial drag internal energy
  double APDUDt; // artificial repel internal energy
  double DRhoDt;
  double DVXDt; //**
  double INDVXDt;
  double ARDVXDt;
  double ADDVXDt; // artificial drag
  double APDVXDt; // artificial repel
  double EXDVXDt;
  double DVYDt; //**
  double INDVYDt;
  double ARDVYDt;
  double ADDVYDt; // artificial drag
  double APDVYDt; // artificail repel
  double EXDVYDt;
  double DVZDt; //**
  double INDVZDt;
  double ARDVZDt;
  double ADDVZDt; // artificial drag
  double APDVZDt; // artificial repel
  double EXDVZDt;
  double TDSDt;

public:
  MkLiuParticle(void) : MkPoint() { Initialize(); }
  MkLiuParticle(double x, double y) : MkPoint(x, y) { Initialize(); }
  MkLiuParticle(double x, double y, double z) : MkPoint(x, y, z) { Initialize(); }
  ~MkLiuParticle(void) {}

  void Initialize();

  void SetupColor();
  void SetColorType(CType ct) { ColorType = ct; }
  void SetR(double r) { R = r; }
  void SetG(double g) { G = g; }
  void SetB(double b) { B = b; }
  void SetRadius(double r) { Radius = r; }
  void SetXVel(double xv) { XVel = xv; }
  void SetYVel(double yv) { YVel = yv; }
  void SetZVel(double yv) { YVel = yv; }
  void SetXAVel(double xv) { XAVel = xv; }
  void SetYAVel(double yv) { YAVel = yv; }
  void SetZAVel(double yv) { YAVel = yv; }
  void SetParticleType(PType pt) { ParticleType = (int)pt; }
  void SetParticleType(int pt) { ParticleType = pt; }
  void SetMass(double mass) { Mass = mass; }
  void SetRho(double rho) { Rho = rho; }
  void SetEta(double eta) { Eta = eta; }
  void SetPress(double press) { Press = press; }
  void SetTemp(double temp) { Temp = temp; }
  void SetEnergy(double u) { Energy = u; }
  void SetSoundSpeed(double c) { SoundSpeed = c; }
  void SetSmoothLen(double h) { SmoothLen = h; }
  void SetCountIac(double iac) { CountIac = iac; }
  void SetGrid(int gi, int gj, int gk)
  {
    GI = gi;
    GJ = gj;
    GK = gk;
  }
  //  void SetDEDt(double d){DEDt = d;}
  void SetDRhoDt(double d) { DRhoDt = d; }
  void SetDVXDt(double d) { DVXDt = d; }
  void SetTDSDt(double d) { TDSDt = d; }
  void SetDUDt(double d) { DUDt = d; }
  //void SetW(double w){W = w;}
  //void SetDWDX(double d){DWDX = d;}

  double GetR() { return R; }
  double GetG() { return G; }
  double GetB() { return B; }
  double GetRadius() { return Radius; }
  double &GetXVel() { return XVel; }
  double &GetYVel() { return YVel; }
  double &GetZVel() { return YVel; }
  double &GetXAVel() { return XAVel; }
  double &GetYAVel() { return YAVel; }
  double &GetZAVel() { return YAVel; }
  int GetParticleType() { return ParticleType; }
  //  PType GetParticleType(){ return ParticleType;}
  double GetMass() { return Mass; }
  double GetRho() { return Rho; }
  double GetEta() { return Eta; }
  double GetPress() { return Press; }
  double GetTemp() { return Temp; }
  double GetEnergy() { return Energy; }
  double GetSoundSpeed() { return SoundSpeed; }
  double GetSmoothLen() { return SmoothLen; }
  int GetGI() { return GI; }
  int GetGJ() { return GJ; }
  int GetGK() { return GK; }

  //  double GetDEDt(){return DEDt;}
  double GetDRhoDt() { return DRhoDt; }
  double GetDVXDt() { return DVXDt; }
  double GetTDSDt() { return TDSDt; }
  double GetDUDt() { return DUDt; }
  //double GetW(){return W;}
  //double GetDWDX(){return DWDX;}

  MkLiuParticle &operator=(const MkLiuParticle &rp);

#if defined(__GL_H__)
  void Draw();
#endif
};

class MkLiuParticles
{
protected:
  boost::shared_array<MkLiuParticle> FParticle;
  MkLiuParticle FCenter;
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
  MkLiuParticles(int size, MkLiuParticle *rps);
  MkLiuParticles(int size);
  MkLiuParticles()
  {
    FSizeOfArray = FSize = 0;
    FParticle = NULL;
  }
  ~MkLiuParticles();
  void Initialize(int size);
  void Initialize(int size, MkLiuParticle *);
  int GetSize() { return FSize; };
  int GetNumber() { return FSize; };
  boost::shared_array<MkLiuParticle> GetParticles() { return FParticle; }
  MkLiuParticle GetCenter()
  {
    FindCenter();
    return FCenter;
  };
  bool Add(MkLiuParticle point); // change of size of point
  bool Add(int index, MkLiuParticle point);
  bool Add(MkLiuParticles &p)
  {
    for (int i = 0; i < p.GetSize(); i++)
      Add(p[i]);
    return true;
  }
  bool Delete(MkLiuParticle point); // change of size of point
  bool Delete(int index);
  int Grow(int Delta);   // change of size of array
  int Shrink(int Delta); // change of size of array
  bool Swap(int i, int j);
  bool Clear();
  void SetColorType(CType ct)
  {
    for (int i = 0; i < FSize; i++)
      FParticle[i].SetColorType(ct);
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
  bool hasParticle(MkLiuParticle &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FParticle[i] == pnt)
        return true;
    return false;
  }
  int numParticle(MkLiuParticle &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FParticle[i] == pnt)
        return i;
    return -1;
  }
  virtual MkLiuParticle &operator[](int);
  MkLiuParticles &operator*=(MkMatrix4<double> &rm);
  friend MkLiuParticles &operator*(MkLiuParticles &rps, MkMatrix4<double> &rm);

  MkLiuParticles &Translate(MkLiuParticle rp);
  MkLiuParticles &Translate(double x, double y, double z);
  MkLiuParticles &Rotate(double alpha, double beta, double gamma);
  MkLiuParticles &RotateInX(double ang);
  MkLiuParticles &RotateInY(double ang);
  MkLiuParticles &RotateInZ(double ang);
  MkLiuParticles &RotateInA(double ang, double l, double m, double n);
  MkLiuParticles &Scale(double sx, double sy, double sz);

  MkLiuParticles &operator=(MkLiuParticles &points);
  bool operator==(MkLiuParticles &points);
  bool operator!=(MkLiuParticles &points) { return !operator==(points); }

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

extern MkLiuParticle NullParticle;
extern MkLiuParticles NullParticles;

#endif
