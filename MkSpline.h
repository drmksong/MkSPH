//---------------------------------------------------------------------------
#ifndef MkSplineH
#define MkSplineH

#include "MkFloat.h"

class MkSpline {
private:
   MkFloat X,Y,Y2;
   float YP1,YPn;
   int   N;
public:
   MkSpline();
   MkSpline(int);
   MkSpline(MkFloat &x,MkFloat &y);
   MkSpline(MkFloat &x,MkFloat &y,float yp1,float ypn);
   void Set(MkFloat &x,MkFloat &y);
   void Set(float yp1,float ypn);
   void Set(MkFloat &x,MkFloat &y,float yp1,float ypn);

   void CalcY2();
   MkSpline & operator=(MkSpline &);
   bool operator==(MkSpline &);
   bool operator!=(MkSpline &);
   float operator [](float);
   float operator ()(float);
};

class MkSplines {
protected:
    MkSpline *FSpline;
    int FSize;
#ifdef __BCPLUSPLUS__
    TColor Color;
#endif
public:
    MkSplines(int FSize);
    MkSplines(){FSize = 0;FSpline = NULL;}
     ~MkSplines();
    virtual void Initialize(int size);
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    bool Clear();
#ifdef __BCPLUSPLUS__
    TColor GetColor(){return Color;};
    void SetColor(TColor c){Color = c;}
#endif 
    MkSpline & operator[](int);
    MkSplines & operator=(MkSplines &splines);
    bool operator==(MkSplines &splines);
#ifdef __BCPLUSPLUS__
    virtual void Draw(TObject *);
#endif
};
//---------------------------------------------------------------------------
#endif
