//---------------------------------------------------------------------------
#ifndef MkPointHPP
#define MkPointHPP
#include <math.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "MkColor.hpp"

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// class MkPaint;
// //typedef MkColor;
// #endif

#include "MkMisc.hpp"
#include "MkMatrix.hpp"

struct MkPoint;
class MkPoints;
extern MkPoint NullPoint;
extern MkPoints NullPoints;

//class MkMatrix4;

struct MkPoint
{
public:
  double X, Y, Z;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif
#if defined(_MSC_VER) && defined(_WINDOWS_)
  MkColor Color;
#endif

#if defined(__GL_H__)
  MkColor Color;
  double DotRadius;
#endif

  MkPoint()
  {
    X = 0;
    Y = 0;
    Z = 0;
  }
  MkPoint(double x, double y)
  {
    X = x;
    Y = y;
    Z = 0;
  }
  MkPoint(double x, double y, double z)
  {
    X = x;
    Y = y;
    Z = z;
  }
  void SetPoint(double x, double y)
  {
    X = x;
    Y = y;
  }
  void SetPoint(double x, double y, double z)
  {
    X = x;
    Y = y;
    Z = z;
  }
  void Set(double x, double y, double z)
  {
    X = x;
    Y = y;
    Z = z;
  }
  void Set(MkPoint rp)
  {
    X = rp.X;
    Y = rp.Y;
    Z = rp.Z;
  }

#ifdef __BCPLUSPLUS__
  TColor GetColor()
  {
    return Color;
  };
  void SetColor(TColor c) { Color = c; }
#endif

#if !defined(_MSC_VER) && !defined(_WINDOWS_) || defined(__BCPLUSPLUS__)
  MkPoint &operator=(const MkPoint &rp);
#endif

#if defined(__GL_H__)
  void SetColor(MkColor c)
  {
    Color = c;
  }
  void SetDotRadius(double r) { DotRadius = r; }
  MkColor GetColor() { return Color; }
  double GetDotRadius() { return DotRadius; }
#endif

  MkPoint &operator=(MkPoint &rp);
  MkPoint &operator+=(MkPoint &rp)
  {
    X += rp.X;
    Y += rp.Y;
    Z += rp.Z;
    return *this;
  }
  friend MkPoint &operator+(MkPoint &a, MkPoint &b)
  {
    static MkPoint c;
    c.X = a.X + b.X;
    c.Y = a.Y + b.Y;
    c.Z = a.Z + b.Z;
    return c;
  }
  MkPoint &operator-=(MkPoint &rp)
  {
    X -= rp.X;
    Y -= rp.Y;
    Z -= rp.Z;
    return *this;
  }
  friend MkPoint &operator-(MkPoint &a, MkPoint &b)
  {
    static MkPoint c;
    c.X = a.X - b.X;
    c.Y = a.Y - b.Y;
    c.Z = a.Z - b.Z;
    return c;
  }
  friend double CalDist(MkPoint sp, MkPoint ep)
  {
    return sqrt(pow(sp.X - ep.X, 2) + pow(sp.Y - ep.Y, 2) + pow(sp.Z - ep.Z, 2));
  }
  double GetAng();
  void GetAng(double &alpha, double &beta, double &gamma);
  bool operator==(MkPoint &rp);
  bool operator!=(MkPoint &rp);

  friend MkPoint &operator*(MkPoint &rp, MkMatrix4<double> &rm)
  {
    static MkPoint rp_t;
    rp_t = rp;

    rp_t.X = rp.X * rm(0, 0) + rp.Y * rm(0, 1) + rp.Z * rm(0, 2) + 1 * rm(0, 3);
    rp_t.Y = rp.X * rm(1, 0) + rp.Y * rm(1, 1) + rp.Z * rm(1, 2) + 1 * rm(1, 3);
    rp_t.Z = rp.X * rm(2, 0) + rp.Y * rm(2, 1) + rp.Z * rm(2, 2) + 1 * rm(2, 3);

    return rp_t;
  }
  friend MkPoint &operator*(MkPoint &rp, double f)
  {
    static MkPoint rp_t;
    rp_t = rp;

    rp_t.X *= f;
    rp_t.Y *= f;
    rp_t.Z *= f;
    return rp_t;
  }
  MkPoint &operator*=(MkMatrix4<double> &rm);
  MkPoint &operator*=(double f);
  MkPoint &operator/(double f)
  {
    static MkPoint p;
    static MkPoint NullPoint;
    p = *this;
    if (fabs(f) < EPS)
    {
      MkDebug("MkPoint::operator/(double) Try to divide with too small double\n");
      return NullPoint;
    }
    p.X = p.X / f;
    p.Y = p.Y / f;
    p.Z = p.Z / f;
    return p;
  };
  MkPoint &operator/(int i)
  {
    static MkPoint p;
    static MkPoint NullPoint;
    p = *this;
    if (i == 0)
    {
      MkDebug("MkPoint::operator/(int) Try to divide with zero\n");
      return NullPoint;
    }
    p.X = p.X / i;
    p.Y = p.Y / i;
    p.Z = p.Z / i;
    return p;
  }

  void Unify();
  MkPoint &Translate(MkPoint &rp);
  MkPoint &Translate(double x, double y, double z);
  MkPoint &Rotate(double alpha, double beta, double gamma);
  MkPoint &RotateInX(double ang);
  MkPoint &RotateInY(double ang);
  MkPoint &RotateInZ(double ang);
  MkPoint &RotateInA(double ang, double l, double m, double n);
  MkPoint &Scale(double sx, double sy, double sz);
  void Normalize();
  bool IsNear(MkPoint &rp) { return CalDist(*this, rp) < 0.001; }
  bool IsNear(double x, double y)
  {
    static MkPoint rp;
    rp.SetPoint(x, y);
    return CalDist(*this, rp) < 0.001;
  }

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif

#if defined(__GL_H__)
  void Draw();
#endif
};

void Swap(MkPoint &p1, MkPoint &p2);

class MkPoints
{
protected:
  boost::shared_array<MkPoint> FPoint;
  MkPoint FCenter;
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

#if defined(__GL_H__)
  MkColor Color;
  double DotRadius;
#endif

public:
  MkPoints(int size, MkPoint *rps);
  MkPoints(int size);
  MkPoints()
  {
    FSizeOfArray = FSize = 0;
    FPoint.reset();
  }
  ~MkPoints();
  virtual void Initialize(int size);
  virtual void Initialize(int size, MkPoint *);
  int GetSize() { return FSize; };
  int GetNumber() { return FSize; };
  MkPoint *GetPoints() { return FPoint.get(); }
  MkPoint GetCenter()
  {
    FindCenter();
    return FCenter;
  };
  bool Add(MkPoint point); // change of size of point
  bool Add(int index, MkPoint point);
  bool Add(MkPoints &p)
  {
    for (int i = 0; i < p.GetSize(); i++)
      Add(p[i]);
    return true;
  }
  bool Delete(MkPoint point); // change of size of point
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

#if defined(__GL_H__)
  void SetColor(MkColor c)
  {
    Color = c;
  }
  void SetDotRadius(double r) { DotRadius = r; }
  MkColor GetColor() { return Color; }
  double GetDotRadius() { return DotRadius; }
#endif

  bool hasPoint(MkPoint &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FPoint[i] == pnt)
        return true;
    return false;
  }
  int numPoint(MkPoint &pnt)
  {
    for (int i = 0; i < FSize; i++)
      if (FPoint[i] == pnt)
        return i;
    return -1;
  }
  virtual MkPoint &operator[](int);
  MkPoints &operator*=(MkMatrix4<double> &rm);
  friend MkPoints &operator*(MkPoints &rps, MkMatrix4<double> &rm)
  {
    static MkPoints rps_t;
    rps_t = rps;
    for (int i = 0; i < rps.FSize; i++)
      rps_t.FPoint[i] = rps.FPoint[i] * rm;
    return rps_t;
  }

  MkPoints &Translate(MkPoint rp);
  MkPoints &Translate(double x, double y, double z);
  MkPoints &Rotate(double alpha, double beta, double gamma);
  MkPoints &RotateInX(double ang);
  MkPoints &RotateInY(double ang);
  MkPoints &RotateInZ(double ang);
  MkPoints &RotateInA(double ang, double l, double m, double n);
  MkPoints &Scale(double sx, double sy, double sz);

  MkPoints &operator=(MkPoints &points);
  bool operator==(MkPoints &points);
  bool operator!=(MkPoints &points) { return !operator==(points); }

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif

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

//---------------------------------------------------------------------------
#endif
