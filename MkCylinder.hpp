//---------------------------------------------------------------------------
#ifndef MkCylinderHPP
#define MkCylinderHPP

#include <math.h>
#include <string>
#include "MkContainer.hpp"
#include "MkShape.hpp"
#include "MkPoint.hpp"
#include "MkTriangle.hpp"
#include "MkPlane.hpp"

//---------------------------------------------------------------------------
class MkPlane;
class MkJointPlane;
class MkPennyJoint;

class MkCylinder : public MkShape
{
private:
  MkTriangles Surf;
  int SurfDivision;
  bool needUpdate;
  bool UpdateSurf();
  std::string className;

protected:
  MkPoint FCP;
  float FRadius;
  float FLength;
  float Fl, Fm, Fn;
  float Psi, Theta; // Psi rotate in y-axis, Beta rotate in z-axis
  MkPoints FPoints;

public:
  MkCylinder();
  MkCylinder(float cx, float cy, float cz, float radius);
  MkCylinder(MkPoint cp, float radius);
#ifdef __BCPLUSPLUS__
  MkCylinder(float cx, float cy, float cz, float radius, TColor C);
  MkCylinder(MkPoint cp, float radius, TColor C);
  AnsiString ClassName() { return AnsiString("MkCylinder"); };
#else
  std::string ClassName()
  {
    return className;
  }
#endif

  void SetCylinder(float cx, float cy, float cz, float radius);
  void SetCylinder(MkPoint cp, float radius);
#ifdef __BCPLUSPLUS__
  void SetCylinder(float cx, float cy, float cz, float radius, TColor C);
  void SetCylinder(MkPoint cp, float radius, TColor C);
#endif
  void SetOrient(MkPoint orient);
  void SetOrient(float l, float m, float n);
  void SetLength(float len) { FLength = len; }

  virtual void SetCenter(float cx, float cy, float cz);
  virtual void SetCenter(MkPoint cp);
  virtual void SetRadius(float radius);
  void SetDivision(int n) { SurfDivision = n; }

  bool IsInSurface(MkPoint &pnt, float thick);
  bool IsInSpace(MkPoint &pnt);
  bool IsIntersect(MkLine &rl);
  bool isCylinder() { return true; }

  MkPoints &CalcIntPnts(MkLine &rl);
  void GetIntParam(MkLine &rl, float &t1, float &t2);

  MkPoint GetCenter();
  float GetRadius();
  float GetDist(MkPoint &pnt);
  MkPoint &operator[](int);

  MkCylinder &operator=(MkCylinder &rs);
  bool operator&&(MkLine &rl);
  MkPoints &operator&(MkLine &rl);
  bool operator==(MkCylinder &rs);
  bool operator!=(MkCylinder &rs);

  void RotateSpace(MkPoint &rp);
  void RotateSpace(MkLine &rl);
  void RotateSpace(MkPlane &rp);
  void RotateSpace(MkJointPlane &jp);
  void RotateSpace(MkPennyJoint &pj);

  void UnRotateSpace(MkPoint &rp);
  void UnRotateSpace(MkLine &rl);
  void UnRotateSpace(MkPlane &rp);
  void UnRotateSpace(MkJointPlane &jp);
  void UnRotateSpace(MkPennyJoint &pj);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

// class MkCylinders : public MkShape
// {
// protected:
//   MkCylinder *FCylinder;
//   int FSize;

// public:
//   MkCylinders(int Size);
//   MkCylinders()
//   {
//     FSize = 0;
//     FCylinder = NULL;
//   }
//   ~MkCylinders()
//   {
//     if (FCylinder)
//     {
//       delete (MkCylinder *)FCylinder;
//       FCylinder = NULL;
//     }
//   }
//   void Initialize(int Size);
//   void Clear();
//   int GetSize() { return FSize; }
//   virtual MkCylinder &operator[](int);
//   MkCylinders &operator=(MkCylinders &cylinders);

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };

typedef MkContainer<MkCylinder> MkCylinders;
template class MkContainer<MkCylinder>;

extern MkCylinder NullCylinder;
extern MkCylinders NullCylinders;

//---------------------------------------------------------------------------
#endif
