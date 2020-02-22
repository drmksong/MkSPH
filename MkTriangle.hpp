//---------------------------------------------------------------------------
#ifndef MkTriangleHPP
#define MkTriangleHPP

#include <stdio.h>
#include <string>
#include "MkContainer.hpp"
#include "MkPoint.hpp"
#include "MkLine.hpp"
#include "MkShape.hpp"

double frandom(double f);

class MkTriangle : public MkShape
{
private:
  MkPoint StartPoint, MidPoint, EndPoint;
  double FArea;
  double FRadius;
  MkVector<double> FNormal;
  double A[3], B[3], C[3], AR2, GradX, GradY;
  std::string className;

  void CalArea();
  void CalGrad();
  void CalRadius();
  void CalNormal();
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif

public:
  MkTriangle();
  MkTriangle(int);
  MkTriangle(MkPoint sp, MkPoint mp, MkPoint ep);
  ~MkTriangle();

  void Reset(MkPoint sp, MkPoint mp, MkPoint ep);
  void Reset(MkPoint rps[3]);

  void Translate(MkPoint rp);
  void Translate(double x, double y, double z);
  double GetArea();
  double GetXGrad(); // x-direction gradient to z value
  double GetYGrad(); // y-direction gradient to z value
  double GetRadius();
  MkVector<double> &GetNormal();
#ifdef __BCPLUSPLUS__
  void SetColor(TColor c)
  {
    Color = c;
  }
  TColor GetColor() { return Color; }
  virtual AnsiString ClassName() { return AnsiString("MkTriangle"); }
#else
  std::string ClassName()
  {
    return className;
  }
#endif
  bool isIntersect(MkTriangle &rt);
  bool isValid();
  bool isIn(double x, double y);
  bool isIn(MkPoint);
  bool isInside(MkPoint);
  bool isTriangle() { return true; }

  MkLine &FirstLine();
  MkLine &SecondLine();
  MkLine &LastLine();

  MkPoint &operator[](int);
  MkLine &operator()(int);
  double operator()(double, double);
  bool operator&&(MkTriangle &rt);
  bool operator==(MkTriangle &rt);
  bool operator!=(MkTriangle &rt);
  MkTriangle &operator=(MkTriangle &rt);
  double CrossProduct();

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};
//---------------------------------------------------------------------------
// class MkTriangles
// {
// protected:
//   MkTriangle *FTriangle;
//   int FSize;
//   int FSizeOfArray;
// #ifdef __BCPLUSPLUS__
//   TColor Color;
// #endif
// public:
//   MkTriangles(int Size);
//   MkTriangles()
//   {
//     FSize = 0;
//     FSizeOfArray = 0;
//     FTriangle = NULL;
//   }
//   ~MkTriangles()
//   {
//     if (FTriangle)
//     {
//       delete (MkTriangle *)FTriangle;
//       FTriangle = NULL;
//     }
//   }
//   void Initialize(int Size);
//   void Clear();

//   bool Add(MkTriangle tri);    // change of size of tri
//   bool Delete(MkTriangle tri); // change of size of tri
//   int Grow(int Delta);         // change of size of array
//   int Shrink(int Delta);       // change of size of array

//   int GetSize() { return FSize; }
//   int SaveUCD(char *filename);
//   int SaveLID(char *filename); // Auto-CAD Lisp Input Data file

// #ifdef __BCPLUSPLUS__
//   void SetColor(TColor c)
//   {
//     Color = c;
//   }
//   TColor GetColor() { return Color; }
// #endif

//   virtual MkTriangle &operator[](int);
//   MkTriangles &operator=(MkTriangles &Triangles);

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

typedef MkContainer<MkTriangle> MkTriangles;
template class MkContainer<MkTriangle>;

extern MkTriangle NullTriangle;
extern MkTriangles NullTriangles;
#endif
