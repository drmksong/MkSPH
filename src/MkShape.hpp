//---------------------------------------------------------------------------
#ifndef MkShapeH
#define MkShapeH
#include <string>
#include "MkPoint.hpp"

class MkShape
{
private:
     //     virtual void CalArea(){};
protected:
     std::string className;

public:
#ifdef __GL_H__
     double R, G, B;
     double Transparancy;
#endif

#ifdef __BCPLUSPLUS__
     TColor Color;
     TPenStyle PenStyle;
#endif

public:
     MkShape()
     {
          className = "MkShape";
          NextShape = NULL;
          PrevShape = NULL;
     };
     //     ~MkShape();
     MkShape *Next() { return NextShape; };
     MkShape *Prev() { return PrevShape; };
#ifdef __BCPLUSPLUS__
     void SetColor(TColor C)
     {
          Color = C;
     }
     TColor GetColor(void) { return Color; }
     void SetStyle(TPenStyle ps) { PenStyle = ps; }
     TPenStyle GetStyle() { return PenStyle; }
#endif

#ifdef __GL_H__
     void SetR(double r)
     {
          R = r;
     }
     void SetG(double g) { G = g; }
     void SetB(double b) { B = b; }
     void SetTransparancy(double t) { Transparancy = t; }
     double GetR() { return R; }
     double GetG() { return G; }
     double GetB() { return B; }
     double GetTransparancy() { return Transparancy; }
#endif

     virtual double GetArea()
     {
          return 0;
     }
#ifdef __BCPLUSPLUS__
     virtual AnsiString ClassName()
     {
          return AnsiString("MkShape");
     }
#else
     std::string ClassName()
     {
          return className;
     }
#endif

     virtual bool isTriangle()
     {
          return false;
     }
     virtual bool isRect() { return false; }
     virtual bool isCube() { return false; }
     virtual bool isCircle() { return false; }
     virtual bool isSphere() { return false; }
     virtual bool isCylinder() { return false; }

     virtual void operator=(MkShape &ms)
     {
#ifdef __BCPLUSPLUS__
          Color = ms.Color;
          PenStyle = ms.PenStyle;
#endif
          NextShape = ms.NextShape;
          PrevShape = ms.PrevShape;
          return;
     }
     virtual bool operator==(MkShape &ms);
     virtual bool operator!=(MkShape &ms);
     virtual bool IsInSurface(MkPoint &pnt, double thick) { return false; }
     virtual bool IsInSpace(MkPoint &pnt) { return false; }
     //     virtual MkShape operator=(MkShape &ms);
public:
     MkShape *NextShape;
     MkShape *PrevShape;

#ifdef __BCPLUSPLUS__
     virtual void Draw(TObject *)
     {
          MkDebug("Pure virtual function MkShape::Draw() is called\n");
     }
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
     virtual void Draw(MkPaint *)
     {
          MkDebug("Pure virtual function MkShape::Draw() is called\n");
     }
#endif
};

//---------------------------------------------------------------------------
// �� class������ � ������ �޸� �Ҵ絵 �Ͼ�� �ʴ´�.
class MkShapeList
{ // �ͳγ��� �� �κ� : ���̴�, ������ ��
private:
     MkShape *FirstMyShapes;
     MkShape *CurrentMyShapes;
     MkShape *LastMyShapes;
     int NumberOfShape;

public:
     MkShapeList(MkShape *);
     MkShapeList();
     ~MkShapeList();
     bool Add(MkShape *);
     bool Insert(MkShape *);
     bool Delete(MkShape *); // �� ������ ����Ʈ&�޸� ���Ŵ� �ƴ�.
     bool Clear();           // ��� ����Ʈ �Ӹ� �ƴ϶� �޸𸮵� ����
     double GetArea();
     int GetNumberOfShape() { return NumberOfShape; }
     MkShape &operator[](int i);
     MkShapeList &operator=(MkShapeList &sl);

#ifdef __BCPLUSPLUS__
     void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
     void Draw(MkPaint *);
#endif
};
//---------------------------------------------------------------------------
#endif
