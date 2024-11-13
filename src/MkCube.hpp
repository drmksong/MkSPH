//---------------------------------------------------------------------------
#ifndef MkCubeH
#define MkCubeH
//---------------------------------------------------------------------------
#include <string>
#include "MkPoint.hpp"
#include "MkShape.hpp"
#include "MkContainer.hpp"

class MkCube : public MkShape
{
protected:
  MkPoint FCenter;
  float Psi, Theta; // degree
  float XLength;
  float YLength;
  float ZLength;
  //  std::string className;

public:
  MkCube();
  MkCube(MkPoint &rp);
#ifdef __BCPLUSPLUS__
  AnsiString ClassName()
  {
    return AnsiString("MkCube");
  };
#else
  std::string ClassName()
  {
    return className;
  };
#endif
  void SetCenter(MkPoint &center)
  {
    FCenter = center;
  };
  void SetRotation(float psi, float theta)
  {
    Psi = psi;
    Theta = theta;
  }
  void SetLength(float xl, float yl, float zl)
  {
    XLength = xl;
    YLength = yl;
    ZLength = zl;
  }
  void SetXLength(float xl) { XLength = xl; }
  void SetYLength(float yl) { YLength = yl; }
  void SetZLength(float zl) { ZLength = zl; }

public:
  float GetNorm(MkPoint &rp);
  MkPoint GetNormPoint(MkPoint &rp);
  float GetXLength() { return XLength; }
  float GetYLength() { return YLength; }
  float GetZLength() { return ZLength; }

public:
  bool IsInSurface(MkPoint &pnt, float thick);
  bool IsInSpace(MkPoint &pnt);
  bool IsIn(MkPoint &pnt);
  bool isCube() { return true; }

  MkPoint operator[](int i);
  MkCube &operator=(MkCube &rc)
  {
    FCenter = rc.FCenter;
    Psi = rc.Psi;
    Theta = rc.Theta;
    XLength = rc.XLength;
    YLength = rc.YLength;
    ZLength = rc.ZLength;
    return (*this);
  }
  bool operator==(MkCube &rc);
  bool operator!=(MkCube &rc);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

// class MkCubes
// {
// private:
//   int FSize;
//   MkCube *FCube;

// public:
//   MkCubes();
//   MkCubes(int);
//   ~MkCubes();
//   bool Initialize(int size);
//   bool Initialize(int size, MkCube *fault);
//   void Clear();

//   MkCube &operator()(int);
//   MkCube &operator[](int);
//   MkCubes &operator=(MkCubes &a);

//   int GetSize() { return FSize; };

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };

//---------------------------------------------------------------------------
typedef MkContainer<MkCube> MkCubes;
//template class MkContainer<MkCube>;

extern MkCube NullCube;
extern MkCubes NullCubes;

#endif
