//---------------------------------------------------------------------------
#ifndef MkArcHPP
#define MkArcHPP

#include <string>
#include <math.h>
#include "MkContainer.hpp"
#include "MkShape.hpp"
#include "MkPoint.hpp"
#include "MkTriangle.hpp"
#include "MkCircle.hpp"

// It is used for only 2 dimensional geometry operation
// I should upgrade this for 3 dimensional...
//---------------------------------------------------------------------------

class MkArc : public MkCircle
{
private:
  MkPoint FStartPoint, FEndPoint;
  double FStartAng, FEndAng;
  double FCrownArea;
  std::string className;
  //    void  CalArea();
  void CalCrownArea();
  void CalTriArea();
  void CalStartPoint();
  void CalEndPoint();
  void CalStartAngle();
  void CalEndAngle();

public:
  MkArc();
  MkArc(int);
  MkArc(MkPoint p1, MkPoint p2, MkLine line);
  MkArc(double cx, double cy, double radius, double start_ang, double end_ang);
  MkArc(MkPoint cp, double radius, double start_ang, double end_ang);
#ifdef __BCPLUSPLUS__
  MkArc(MkPoint cp, double radius, double start_ang, double end_ang, TColor C);
  MkArc(double cx, double cy, double radius, double start_ang, double end_ang, TColor C);
#endif
  void ReCalcPoint();
  void ReCalcAng();
  void SetCenter(double cx, double cy);
  void SetCenter(MkPoint cp);
  void SetRadius(double radius);
  void SetStartAng(double start_ang);
  void SetEndAng(double end_ang);
  double GetStartAng();
  double GetEndAng();
  double GetArea();
  double GetAngle() { return FEndAng - FStartAng; };
  double GetCrownArea();
  double GetTriArea();
  bool isWithInArc(MkPoint rp);
  bool isWithInAng(MkPoint rp);
  double CrossProduct();
  MkPoint StartPoint();
  MkPoint EndPoint();
  MkPoint &operator[](int);
  MkLine operator()(int);
#ifdef __BCPLUSPLUS__
  AnsiString ClassName()
  {
    return AnsiString("MkArc");
  }
#else
  std::string ClassName()
  {
    return className;
  }
#endif
  MkArc &operator=(MkArc &ra);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

extern MkArc NullArc;

typedef MkContainer<MkArc> MkArcs;

//---------------------------------------------------------------------------

#endif