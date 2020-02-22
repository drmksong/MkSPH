//---------------------------------------------------------------------------
#ifndef MkPolygonHPP
#define MkPolygonHPP

#include "MkArray.hpp"
#include "MkPoint.hpp"
#include "MkMatrix.hpp"
#include "MkTriangle.hpp"
#include "MkCircle.hpp"
#include "MkLine.hpp"

// How can I distingush between NullPoint and Origin Point.
// if two connected point are the same and is NullPoint then
// isFilled returns false.

enum MkBoolType
{
  btUni,
  btInt,
  btSub
}; // union, intersection, subtraction
class MkPolygons;
class MkPolygon : public MkPoints
{
private:
  std::string className;

protected:
  bool Closeness; // is this polygon closed(true) or opened(false);
  bool Convexity;
  bool Crossing; // if lines are crossed(true) each other or not(false)?
  bool Fullness;
  bool isLengthChanged;
  bool isAreaChanged;
  bool isCrossingChanged;
  bool isClosenessChanged;
  bool isConvexityChanged;

  double FLength;
  double FArea;
  MkArray<int> PointState;
  int FCurrent;

  void CheckConvexity();
  void CheckCrossing();
  void CheckFullness();
  double CalLength();

public:
  MkPolygon(int size, MkPoint *);
  MkPolygon(int size); // empty polygon
  MkPolygon();         // even memory is not allocated
  ~MkPolygon();
#ifdef __BCPLUSPLUS__
  AnsiString ClassName()
  {
    return AnsiString("MkPolygon");
  };
#else
  std::string ClassName()
  {
    return className;
  };
#endif

  void Initialize(int size);
  void Initialize(int size, MkPoint *);
  bool Clear()
  {
    isLengthChanged = true;
    isAreaChanged = true;
    return MkPoints::Clear();
  }
  bool Add(MkPoint point)
  {
    isLengthChanged = true;
    isAreaChanged = true;
    return MkPoints::Add(point);
  }
  bool Add(int index, MkPoint point)
  {
    isLengthChanged = true;
    isAreaChanged = true;
    return MkPoints::Add(index, point);
  }
  bool Add(MkPoints &p)
  {
    isLengthChanged = true;
    isAreaChanged = true;
    return MkPoints::Add(p);
  }
  bool AddInBetween(MkPoints &p);
  bool Delete(MkPoint point)
  {
    isLengthChanged = true;
    isAreaChanged = true;
    return MkPoints::Delete(point);
  }
  bool Delete(int index)
  {
    isLengthChanged = true;
    isAreaChanged = true;
    return MkPoints::Delete(index);
  }

  void SetCloseness(bool b) { Closeness = b; };

  bool GetCloseness() { return Closeness; }
  double GetLength();
  double GetArea();
  int GetAlivePoint();

  bool isClosed() { return Closeness; }
  bool isFilled() { return Fullness; }; // all the data is not NullPoint?
  bool isCrossing() { return Crossing; };
  bool isCrossWith(MkLine &l);
  bool isCrossWithX(double x);
  bool isCrossWithY(double y);
  bool IsIn(MkPoint rp);
  bool IsIn(MkCircle &rc);

  void getCrossWith(MkLine &l, MkPoints &pnts);
  void getCrossWithX(double x, MkPoints &pnts);
  void getCrossWithY(double y, MkPoints &pnts);
  int Current() { return FCurrent; };
  void SetCurrent(int i)
  {
    if (i >= 0 && i < GetSize())
      FCurrent = i;
  };
  int Next();
  int Prev();
  int Next(int);
  int Prev(int);
  int AlivedNext();
  int AlivedPrev();
  int AlivedNext(int);
  int AlivedPrev(int);
  MkPoint &CurrentPoint() { return (*this)[FCurrent]; };
  MkPoint &NextPoint();
  MkPoint &PrevPoint();
  MkPoint &NextPoint(int);
  MkPoint &PrevPoint(int);
  MkPoint &NextPoint(MkPoint &pnt);
  MkPoint &PrevPoint(MkPoint &pnt);
  MkPoint &AlivedNextPoint();
  MkPoint &AlivedPrevPoint();
  MkPoint &AlivedNextPoint(int);
  MkPoint &AlivedPrevPoint(int);

  void Offset(MkPolygon &poly, bool dir, double space);
  void Offset(MkPolygon &poly, MkVector<double> &up, bool dir, double space);
  MkPoint Measure(double dis);
  double Measure(MkPoint pnt); // 0 ~ 1
  bool InverseDirection();
  bool Merge(MkPolygon &poly);
  bool FindInter(MkPolygon &poly, MkPoints &pnts);
  MkVector<double> &GetVector();

  int FindPoly(MkPolygon &b, MkBoolType bt);
  bool BoolSub(MkPolygon &b, MkPolygons &c); // C = A - B subtraction
  bool BoolInt(MkPolygon &b, MkPolygons &c); // C = A * B intersection
  bool BoolUni(MkPolygon &b, MkPolygons &c); // C = A + B union
  void Extract(MkPoints &b, MkPolygon &c);

  MkPoint &operator[](int i);
  MkLine &operator()(int i);
  MkPolygon &operator=(MkPolygon &polygon);
  bool operator!=(MkPolygon &poly);
  double CalArea();
  double CalArea2();

  bool Out(char *fname);
#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

class MkPolygons
{
protected:
  MkPolygon *FPolygon;
  int FSize; //Actual size of polys
  int FSizeOfArray;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif

public:
  MkPolygons(int size, MkPolygon *poly);
  MkPolygons(int size);
  MkPolygons()
  {
    FSizeOfArray = FSize = 0;
    FPolygon = NULL;
  }
  ~MkPolygons();
  virtual void Initialize(int size);
  virtual void Initialize(int size, MkPolygon *);
  int GetSize() { return FSize; };
  int GetNumber() { return FSize; };
  bool Add(MkPolygon &poly); // change of size of poly
  bool Add(MkPolygons &poly)
  {
    bool flag = true;
    MkDebug("MkPolygons::Add(Polygons) is called.\n");
    for (int i = 0; i < poly.GetSize(); i++)
      flag = Add(poly[i]) && flag;
    return flag;
  } // change of size of poly
  bool Add(int index, MkPolygon &poly);
  bool Delete(MkPolygon &poly); // change of size of poly
  bool Delete(int index);
  int Grow(int Delta);   // change of size of array
  int Shrink(int Delta); // change of size of array
  bool Clear();
#ifdef __BCPLUSPLUS__
  TColor GetColor()
  {
    return Color;
  };
  void SetColor(TColor c) { Color = c; }

#endif
  virtual MkPolygon &operator[](int);
  MkPolygons &operator=(MkPolygons &polys);
  bool operator==(MkPolygons &polys);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

void GetSubPolygon(double ymin, double ymax, MkPolygon &inpoly, MkPolygon &outpoly);
void GetSubParam(int i, MkPolygon &in, double &aj, double &bj, double &lj1, double &lj);
//---------------------------------------------------------------------------
extern MkPoint NullPoint;
extern MkLine NullLine;
extern MkPolygon NullPolygon;

#endif
