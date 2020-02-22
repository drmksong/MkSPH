//---------------------------------------------------------------------------
#ifndef MkPlaneHPP
#define MkPlaneHPP

#include <string>
#include <algorithm>
#include <assert.h>
#include <boost/shared_array.hpp>

#include "MkContainer.hpp"
#include "MkShape.hpp"
#include "MkPoint.hpp"
#include "MkLine.hpp"
#include "MkMisc.hpp"
#include "MkMatrix.hpp"
#include "MkArray.hpp"
#include "MkTriangle.hpp"
#include "MkSphere.hpp"
#include "MkCylinder.hpp"
#include "MkPolygon.hpp"

// plane primitive�� (-1,-1),(1,-1),(1,1),(-1,1)�� ������, ȸ��, �̵��� ������
// ���ؼ� ���ϴ� �簢���� ����� ����� �� �� �ִ�.
//---------------------------------------------------------------------------
class MkPlane : public MkShape
{
protected:
  MkPoint LD, RD, LU, RU;
  MkLine Edge[4];
  MkPoint FCenter;
  double ScaleX;
  double ScaleY;
  double Alpha, Beta, Gamma; // �� x,y,z �������� ȸ��, ������ x,y,z ��
  double A, B, C, D;         // Ax + By + Cz = D
  bool isFinite;
  std::string className;

public:
  MkPlane(); // ������� ��ǥ
  MkPlane(int i);
  MkPlane(MkPoint rp1, MkPoint rp2, MkPoint rp3);
  MkPlane(MkPoint *rps);

#ifdef __BCPLUSPLUS__
  AnsiString ClassName()
  {
    return AnsiString("MkPlane");
  }
#else
  std::string ClassName()
  {
    return className;
  }
#endif

  void ResetAll();
  void ResetScale();
  void ResetRotate();
  void ResetTranslate();
  void SetPoints(MkPoint rp1, MkPoint rp2, MkPoint rp3);
  void SetPoints(MkPoint *rps);
  void SetScale(double scalex, double scaley);
  void SetFiniteness(bool flag) { isFinite = flag; }
  void SetRotate(double Alpha, double Beta, double Gamma); // ������ �����ִ� ��...
  void SetTranslate(MkPoint trans);
  void SetTranslate(double x, double y, double z);
  void SetCenter(MkPoint center) { SetTranslate(center); }
  void SetCenter(double x, double y, double z) { SetTranslate(x, y, z); }

  void Rotate(double alpha, double beta, double gamma); // ������ ���� ȸ��
  void RotateInX(double ang);
  void RotateInY(double ang);
  void RotateInZ(double ang);
  void RotateInA(double ang, double l, double m, double n);
  void RotateCen(double alpha, double beta, double gamma); // Trans�� ���� ȸ��
  void RotateInXCen(double ang);
  void RotateInYCen(double ang);
  void RotateInZCen(double ang);
  void RotateInACen(double ang, double l, double m, double n);

  void Translate(double x, double y, double z)
  {
    FCenter.X += x;
    FCenter.Y += y;
    FCenter.Z += z;
  }
  void Translate(MkPoint rp) { FCenter = FCenter + rp; }

  double GetA(void) { return A; }
  double GetB(void) { return B; }
  double GetC(void) { return C; }
  double GetD(void) { return D; }

  bool GetIntParam(MkLine &rl, double &t1, double &t2);
  bool GetFiniteness() { return isFinite; }
  MkPoint GetTrans() { return FCenter; }
  MkPoint GetCenter() { return FCenter; }
  double GetDistance(MkPoint &pnt)
  {
    assert(A * A + B * B + C * C > 0.0001);
    return (A * pnt.X + B * pnt.Y + C * pnt.Z - D) / sqrt(A * A + B * B + C * C);
  }
  MkPoint GetNearestPoint(MkPoint &pnt)
  {
    double dist = GetDistance(pnt);
    return MkPoint(pnt.X + A * dist, pnt.Y + B * dist, pnt.Z + C * dist);
  }

  void ApplyMatrix();

  MkPoint LeftDown() { return LD; }
  MkPoint RightDown() { return RD; }
  MkPoint LeftUp() { return LU; }
  MkPoint RightUp() { return RU; }

  void CalcABCD();
  bool IsIntersect(MkLine); // updated at 2001.8.25.
  bool IsIntersect(MkTriangle &rt);
  bool IsInPlane(MkLine rl);
  bool IsIn(MkPoint rp);
  bool IsIn(MkLine rl);
  bool IsInSurface(MkPoint &pnt, double thick);
  bool IsInnerSpace(MkPoint &pnt);

  MkPoint &CalcIntPnt(MkLine);
  MkLine &CalcIntLine(MkTriangle &rt);
  MkPoint &CalcIntPnt(MkTriangle &rt);
  MkPoint &CalcIntPnt(MkTriangles &rts);

  void operator*(MkMatrix4<double> &rm);

  bool operator==(MkPlane &);
  bool operator!=(MkPlane &);
  MkPlane &operator=(MkPlane &rp);

  MkPoint &operator[](int i);
  MkLine &operator()(int i);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};
//---------------------------------------------------------------------------

// MkPointsPlane �� ������ ��� Ȥ�� ���� �����ϰ� �� ������ z��ǥ��
// ù��°�� z��ǥ�� �������� ������ ����κ��� ���ȴ�. �̰��� drawing�� ����
// ����ü�� ����� ���� �ƴϱ� ������ MkShape���κ��� inherite���� �ʴ´�.
// ó�� �����ÿ��� ��ǥ���� ���� �� ���並 ǥ���ϱ� ���� �Է¿� ���ǵ����Ͽ���.
// MkPoints�� 3���� ��ǥ�� ǥ���ϹǷ� z���� ������ ������ ù��° point��
// ������ ������ point�� z���� �����Ѵ�.
// �ʱ�������� ������ �� �ִ� ����� ������ �����̳��� ������ �͸� ������
// Z ���� �ٲٴ� ����� ������ �� �ִ�. ���� �Է�ġ�� ����� ���� ���� ��쿡��
// �ǵ��� ��� ���� ���� �� �ְ� ���� ���������� ���ľ� �ϰ� ��Ȯ�� �Է�ġ��
// ��Ⱑ ���� �� �� �ִ�.

enum THeightMode
{
  hmAbs,
  hmCut,
  hmPaste
};
// hmAbs : ������ ������ ���ؼ��� ������ ������ ��� �Ѵ�.
// hmCut : ���� ������ ���ؼ��� ���ݼ����� �������� ������ ���߰� �������� ���д�.
// hmPaste : ���䱸���� ���ؼ��� ���ݼ����� ���� ���� ������ ���̰� ���� ���� ���д�.
// hmAbs�� ���� ������ ���ؼ� �ϴ� ���� �ٶ����ϸ�, hmCut�� ���䰡 ���ݼ���
// �´پ� ���� ���, ��谡 �Ǵ� ��ǥ�� ã�� ���� �ʱ� ������ ������ ��Ȯ�� ��ǥ��
// ���ݼ� ���� �������ϰ� ������ ��Ƽ� �ϸ� ���� ���� �ذ��� �� �ִ�.
// ���䰡 ����� �����ϴ� ��쿡�� ?? -> ���� �ذ�å : �ϴ� ���� ���������� ����
// �����ϴ� ����� ������ �� �ִ�.

class MkPointsPlane : public MkPoints
{
protected:
  MkOrient FOrient;
  THeightMode FHeightMode;
  double Alpha, Beta, Gamma; // X, Y, Z ����������� ȸ������(deg)
  double A, B, C, D;         // Ax + By + Cz = D
public:
  MkPointsPlane();
  MkPointsPlane(int size);
  MkPointsPlane(int size, MkPoint *rps);
  MkPointsPlane(int size, MkPoint *rps, MkOrient fo);
  void SetOrient(MkOrient fo);

  void SetHeightMode(THeightMode fhm) { FHeightMode = fhm; }
  THeightMode GetHeightMode() { return FHeightMode; }
  bool operator==(MkPointsPlane &);
  bool operator!=(MkPointsPlane &);

  bool IsIn(MkPoint rp);         // ���� x-y �����ǥ�� ������ �Ǵ��Ѵ�.
  double CalcHeight(MkPoint rp); // ����� ���������κ��� ù��° ���� ���⼺�� ������
                                 // �������� �������� ������� ���̸� ����� �Ѵ�.
                                 // ���� ����� �������̶�� ��� �� ���ΰ�?

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
//---------------------------------------------------------------------------
class MkJointPlane : public MkPlane
{
protected:
  double Dip, DipDir; // y axis = north, x axis = east.
  double Aperture;
  int Sort, Form, Condition;

public:
  MkJointPlane();
  MkJointPlane(double dip, double dipdir);
  void SetDipNDir(double dip, double dipdir);
  void SetDipDir(double dipdir) { DipDir = dipdir; }
  void SetDip(double dip) { Dip = dip; }
  void SetAperture(double a) { Aperture = a; }
  void SetSort(int s) { Sort = s; }
  void SetForm(int f) { Form = f; }
  void SetCondition(int c) { Condition = c; }
  double GetAperture() { return Aperture; };
  double GetDipDir() { return DipDir; }
  double GetDip() { return Dip; }
  MkOrient GetOrient() { return MkOrient(Dip, DipDir); }
  int GetSort() { return Sort; }
  void CalcABCD();
  MkJointPlane &operator=(MkJointPlane &p)
  {
    MkPlane::operator=((MkPlane &)p);
    Dip = p.Dip;
    DipDir = p.DipDir;
    Aperture = p.Aperture;
    Sort = p.Sort;
    Form = p.Form;
    Condition = p.Condition;
    return *this;
  }
  bool operator==(MkJointPlane &);
  bool operator!=(MkJointPlane &);
};
//---------------------------------------------------------------------------
class MkCylinder;
class MkPennyJoint : public MkJointPlane
{
protected:
  double FRadius;
  MkLine FLine; // for internal use
public:
  MkPennyJoint();
  MkPennyJoint(double dip, double dipdir);
  void SetRadius(double r) { FRadius = r; }
  double GetRadius() { return FRadius; }
  bool GetOutline(MkPolygon &poly, int div);
  bool IsIntersect(MkLine); // updated at 2001.8.25.
  bool IsIntersect(MkTriangle &rt);
  bool IsIntersect(MkPennyJoint &pj);
  bool IsIntersect(MkPlane &rp);
  bool IsIntersect(MkCylinder &rc);
  MkLine &CalcIntLine(MkTriangle &rt);
  MkLine &CalcIntLine(MkPennyJoint &pj);
  MkLine &CalcIntLine(MkPlane &rp);
  MkPoint &CalcIntPnt(MkLine &);
  MkPoint &CalcIntPnt(MkCylinder &rc);
  MkPoint &CalcIntPnt(MkTriangle &rt);
  MkPoint &CalcIntPnt(MkTriangles &rts);

  bool operator&&(MkLine &rl) { return IsIntersect(rl); }
  bool operator&&(MkTriangle &rt) { return IsIntersect(rt); }
  bool operator&&(MkPennyJoint &pj) { return IsIntersect(pj); }
  bool operator&&(MkPlane &rp) { return IsIntersect(rp); }
  bool operator&&(MkCylinder &rc) { return IsIntersect(rc); }
  MkPoint &operator&(MkLine &rl) { return CalcIntPnt(rl); }
  MkLine &operator&(MkPennyJoint &pj) { return CalcIntLine(pj); }
  MkLine &operator&(MkPlane &rp) { return CalcIntLine(rp); }
  MkPoint &operator&(MkCylinder &rc) { return CalcIntPnt(rc); }
  MkPoint &operator&(MkTriangle &rt) { return CalcIntPnt(rt); }
  MkPoint &operator&(MkTriangles &rts) { return CalcIntPnt(rts); }
  MkPennyJoint &operator=(MkPennyJoint &p)
  {
    MkJointPlane::operator=((MkJointPlane &)p);
    FRadius = p.FRadius;
    FLine = p.FLine;
    return *this;
  }
  bool operator==(MkPennyJoint &);
  bool operator!=(MkPennyJoint &);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

//---------------------------------------------------------------------------
template <class T>
class MkPlaneContainer
{
protected:
  boost::shared_array<T> FPlane;
  int FSize;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif
public:
  MkPlaneContainer(int size, T *jp)
  {
    if (size < 0)
    {
      MkDebug("::MkPlaneContainer - MkPlaneContainer(int size)");
      throw Size(std::string("MkPlaneContainer<T>::size below zero strange error"), size);
    }
    if (jp == NULL)
    {
      MkDebug("::MkPlaneContainer - MkPlaneContainer(int size)");
      throw Size(std::string("MkPlaneContainer<T>::attempt copy null pointer"), size);
    }

    FSize = size;

    if (FSize == 0)
    {
      FPlane.reset();
      return;
    }

    try
    {
      FPlane.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkPlaneContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::Constructor throw bad_alloc exception"));
    }

    for (int i = 0; i < FSize; i++)
      FPlane[i] = jp[i];
  }

  MkPlaneContainer(int size)
  {
    if (size < 0)
    {
      MkDebug("::MkPlaneContainer - MkPlaneContainer(int size)");
      throw Size(std::string("MkPlaneContainer<T>::size below zero strange error"), size);
    }

    FSize = size;
    if (FSize == 0)
    {
      FPlane.reset();
      return;
    }

    try
    {
      FPlane.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkPlaneContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::Constructor throw bad_alloc exception"));
    }
  }
  MkPlaneContainer()
  {
    FSize = 0;
    FPlane.reset();
  }
  ~MkPlaneContainer()
  {
    FSize = 0;
    FPlane.reset();
  }
  virtual void Initialize(int size)
  {
    if (size < 0)
    {
      MkDebug("::MkPlaneContainer - Initialize(int size)");
      throw Size(std::string("MkPlaneContainer<T>::Initialize size below zero strange error"), size);
    }
    if (FSize == size)
      return;

    FSize = size;
    if (FSize == 0)
    {
      FPlane.reset();
      return;
    }

    try
    {
      FPlane.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkPlaneContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::Constructor throw bad_alloc exception"));
    }
  }
  virtual void Initialize(int size, T *jp)
  {
    if (size < 0)
    {
      MkDebug("::MkPlaneContainer - Initialize(int size)");
      throw Size(std::string("MkPlaneContainer<T>::size below zero strange error"), size);
    }
    if (jp == NULL)
    {
      MkDebug("::MkPlaneContainer - Initialize(int size)");
      throw Size(std::string("MkPlaneContainer<T>::Initialize attempt copy null pointer"), size);
    }

    FSize = size;
    if (FSize == 0)
    {
      FPlane.reset();
      return;
    }

    try
    {
      FPlane.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkPlaneContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::Constructor throw bad_alloc exception"));
    }

    for (int i = 0; i < FSize; i++)
      FPlane[i] = jp[i];
  }
  int GetSize() { return FSize; };
  int GetNumber() { return FSize; };
  T *GetReal() { return FPlane; }

  bool Clear()
  {
    FSize = 0;
    FPlane.reset();

    return true;
  }

  void Translate(double x, double y, double z)
  {
    for (int i = 0; i < GetSize(); i++)
      FPlane[i].Translate(x, y, z);
  }

  void Translate(MkPoint rp)
  {
    for (int i = 0; i < GetSize(); i++)
      FPlane[i].Translate(rp);
  }

  void Rotate(double alpha, double beta, double gamma)
  {
    for (int i = 0; i < GetSize(); i++)
      FPlane[i].Rotate(alpha, beta, gamma);
  }

  void RotateInX(double ang)
  {
    for (int i = 0; i < GetSize(); i++)
      FPlane[i].RotateInX(ang);
  }

  void RotateInY(double ang)
  {
    for (int i = 0; i < GetSize(); i++)
      FPlane[i].RotateInY(ang);
  }

  void RotateInZ(double ang)
  {
    for (int i = 0; i < GetSize(); i++)
      FPlane[i].RotateInZ(ang);
  }

  void RotateInA(double ang, double l, double m, double n)
  {
    for (int i = 0; i < GetSize(); i++)

      FPlane[i].RotateInA(ang, l, m, n);
  }
#ifdef __BCPLUSPLUS__
  TColor GetColor()
  {
    return Color;
  };
  void SetColor(TColor c) { Color = c; }
#endif
  T &operator[](int i)
  {
    if (FSize == 0)
    {
      MkDebug("MKPlaneContainer::operator[] when FSize == 0");
      throw Range(std::string("MkPlaneContainer::operator[] attempt access non-existing memory "), i);
    }
    else if (i >= 0 && i < FSize)
    {
      return FPlane[i];
    }
    else
    {
      MkDebug("MkPlaneContainer::operator[] range error");
      throw Range(std::string("MkPlaneContainer::operator[] range error"), i);
    }
  }

  MkPlaneContainer &operator=(MkPlaneContainer &jps)
  {
    int i;

    Clear();
    FSize = jps.FSize;
    if (FSize == 0)
    {
      FPlane.reset();
      return *this;
    }

    try
    {
      FPlane.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkPlaneContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::Constructor throw bad_alloc exception"));
    }

    for (i = 0; i < FSize; i++)
      this->FPlane[i] = jps.FPlane[i];

    return *this;
  }

  bool operator==(MkPlaneContainer<T> &p)
  {
    int i;

    if (FSize != p.FSize)
    {
      return false;
    }
    for (i = 0; i < FSize; i++)
    {
      if (this->FPlane[i] != p.FPlane[i])
      {
        return false;
      }
    }

    return true;
  }

#ifdef __BCPLUSPLUS__
  void Draw(TObject *Sender)
  {
    for (int i = 0; i < GetSize(); i++)
    {
      FPlane[i].Draw(Sender);
    }
  }
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *pb)
  {
    for (int i = 0; i < GetSize(); i++)
    {
      FPlane[i].Draw(pb);
    }
  }
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

typedef MkPlaneContainer<MkPlane> MkPlanes;
typedef MkPlaneContainer<MkPointsPlane> MkPointsPlanes;
typedef MkPlaneContainer<MkJointPlane> MkJointPlanes;
typedef MkPlaneContainer<MkPennyJoint> MkPennyJoints;

// //---------------------------------------------------------------------------
// class MkPlanes
// {
// protected:
//   boost::share_array<MkPlane> FPlane;
//   int FSize;
// #ifdef __BCPLUSPLUS__
//   TColor Color;
// #endif
// public:
//   MkPlanes(int size, MkPlane *jp);
//   MkPlanes(int FSize);
//   MkPlanes()
//   {
//     FSize = 0;
//     FPlane = NULL;
//   }
//   ~MkPlanes();
//   virtual void Initialize(int size);
//   virtual void Initialize(int size, MkPlane *jp);
//   int GetSize() { return FSize; };
//   int GetNumber() { return FSize; };
//   MkPlane *GetReal() { return FPlane; }
//   bool Clear();
//   void Translate(double x, double y, double z);
//   void Translate(MkPoint rp);
//   void Rotate(double alpha, double beta, double gamma);
//   void RotateInX(double ang);
//   void RotateInY(double ang);
//   void RotateInZ(double ang);
//   void RotateInA(double ang, double l, double m, double n);
// #ifdef __BCPLUSPLUS__
//   TColor GetColor()
//   {
//     return Color;
//   };
//   void SetColor(TColor c) { Color = c; }
// #endif
//   virtual MkPlane &operator[](int);
//   MkPlanes &operator=(MkPlanes &);
//   bool operator==(MkPlanes &);

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };

// //---------------------------------------------------------------------------
// #ifdef __BCPLUSPLUS__
// #endif

// class MkPointsPlanes
// {
// protected:
//   boost::shared_array<MkPointsPlane> FPoints;
//   int FSize;
// #ifdef __BCPLUSPLUS__
//   TColor Color;
// #endif

// public:
//   MkPointsPlanes(int size, MkPointsPlane *jp);
//   MkPointsPlanes(int FSize);
//   MkPointsPlanes()
//   {
//     FSize = 0;
//     FPoints = NULL;
//   }
//   ~MkPointsPlanes();
//   virtual void Initialize(int size);
//   virtual void Initialize(int size, MkPointsPlane *jp);
//   int GetSize() { return FSize; };
//   int GetNumber() { return FSize; };
//   MkPointsPlane *GetPoints() { return FPoints; }
//   bool Clear();

// #ifdef __BCPLUSPLUS__
//   TColor GetColor()
//   {
//     return Color;
//   };
//   void SetColor(TColor c) { Color = c; }
// #endif

//   virtual MkPointsPlane &operator[](int);

//   MkPointsPlanes &operator=(MkPointsPlanes &);
//   bool operator==(MkPointsPlanes &);

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };

// //---------------------------------------------------------------------------

// class MkJointPlanes
// {
// protected:
//   boost::shared_array<MkJointPlan> e FJoint;
//   int FSize;
// #ifdef __BCPLUSPLUS__
//   TColor Color;
// #endif
// public:
//   MkJointPlanes(int size, MkJointPlane *jp);
//   MkJointPlanes(int FSize);
//   MkJointPlanes()
//   {
//     FSize = 0;
//     FJoint = NULL;
//   }
//   ~MkJointPlanes();
//   virtual void Initialize(int size);
//   virtual void Initialize(int size, MkJointPlane *jp);
//   int GetSize() { return FSize; };
//   int GetNumber() { return FSize; };
//   MkJointPlane *GetJoint() { return FJoint; }
//   bool Clear();
// #ifdef __BCPLUSPLUS__
//   TColor GetColor()
//   {
//     return Color;
//   };
//   void SetColor(TColor c) { Color = c; }
// #endif
//   virtual MkJointPlane &operator[](int);
//   MkJointPlanes &operator=(MkJointPlanes &);
//   bool operator==(MkJointPlanes &);

//   void Translate(double x, double y, double z);
//   void Translate(MkPoint rp);

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };

// //---------------------------------------------------------------------------

// class MkPennyJoints
// {
// protected:
//   boost::shared_array<MkPennyJoint> FPenny;
//   int FSize;
// #ifdef __BCPLUSPLUS__
//   TColor Color;
// #endif
// public:
//   MkPennyJoints(int size, MkPennyJoint *jp);
//   MkPennyJoints(int FSize);
//   MkPennyJoints()
//   {
//     FSize = 0;
//     FPenny = NULL;
//   }
//   ~MkPennyJoints();
//   virtual void Initialize(int size);
//   virtual void Initialize(int size, MkPennyJoint *jp);
//   int GetSize() { return FSize; };
//   int GetNumber() { return FSize; };
//   MkPennyJoint *GetJoint() { return FPenny; }
//   bool Clear();
// #ifdef __BCPLUSPLUS__
//   TColor GetColor()
//   {
//     return Color;
//   };
//   void SetColor(TColor c) { Color = c; }
// #endif
//   virtual MkPennyJoint &operator[](int);

//   MkPennyJoints &operator=(MkPennyJoints &);
//   bool operator==(MkPennyJoints &);

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };
//---------------------------------------------------------------------------
extern MkJointPlane NullJoint;
extern MkJointPlanes NullJoints;
extern MkPointsPlane NullPointPlane;
extern MkPlane NullPlane;
extern MkPennyJoint NullPenny;

#endif

// Scale->Rotate->Translate
