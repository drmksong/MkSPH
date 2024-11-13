//---------------------------------------------------------------------------
#ifndef MkCircleHPP
#define MkCircleHPP

#include <string>
#include <math.h>
#include "MkContainer.hpp"
#include "MkShape.hpp"
#include "MkPoint.hpp"
#include "MkTriangle.hpp"

// It is used for only 2 dimensional geometry operation
// I should upgrade this for 3 dimensional...
//---------------------------------------------------------------------------

class MkCircle : public MkShape
{
protected:
    MkPoint FCP;
    double FRadius;
    double FCircleArea;
    MkPoints FRealPoints;
    virtual void CalArea();
    std::string className;

public:
    MkCircle();
    MkCircle(int);
    MkCircle(double cx, double cy, double radius);
    MkCircle(MkPoint cp, double radius);
#ifdef __BCPLUSPLUS__
    MkCircle(double cx, double cy, double radius, TColor C);
    MkCircle(MkPoint cp, double radius, TColor C);
#endif
    void SetCircle(double cx, double cy, double radius);
    void SetCircle(MkPoint cp, double radius);
#ifdef __BCPLUSPLUS__
    void SetCircle(double cx, double cy, double radius, TColor C);
    void SetCircle(MkPoint cp, double radius, TColor C);
#endif
    virtual void SetCenter(double cx, double cy);
    virtual void SetCenter(MkPoint cp);
    virtual void SetRadius(double radius);

    MkPoint &GetCenter() { return FCP; }
    double GetRadius() { return FRadius; }
    double GetArea()
    {
        CalArea();
        return FCircleArea;
    }
    MkPoint &operator[](int);
#ifdef __BCPLUSPLUS__
    AnsiString ClassName()
    {
        return AnsiString("MkCircle");
    }
#else
    std::string ClassName()
    {
        return className;
    }
#endif

    bool isCircle()
    {
        return true;
    }
    bool IsInSurface(MkPoint &pnt, double thick);
    bool IsInSpace(MkPoint &pnt);

    MkCircle &operator=(MkCircle &rc);
    bool operator&&(MkLine &rl);
    MkPoints operator&(MkLine &rl);
    MkPoints &operator&(MkPoint &rp);

    bool operator==(MkCircle &c);
    bool operator!=(MkCircle &c);

#ifdef __BCPLUSPLUS__
    void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
    void Draw(MkPaint *);
#endif
};

// typedef MkContainer<MkCircle> MkCircles;

// extern MkCircle NullCircle;

// class MkCircles : public MkShape
// {
// protected:
//     MkCircle *FCircle;
//     int FSize;

// public:
//     MkCircles(int Size);
//     MkCircles()
//     {
//         FSize = 0;
//         FCircle = NULL;
//     }
//     ~MkCircles()
//     {
//         if (FCircle)
//         {
//             delete (MkCircle *)FCircle;
//             FCircle = NULL;
//         }
//     }
//     void Initialize(int Size);
//     void Clear();
//     int GetSize() { return FSize; }
//     virtual MkCircle &operator[](int);
//     MkCircles &operator=(MkCircles &circles);

// #ifdef __BCPLUSPLUS__
//     void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//     void Draw(MkPaint *);
// #endif
// };

typedef MkContainer<MkCircle> MkCircles;
//template class MkContainer<MkCircle>;

extern MkCircle NullCircle;
extern MkCircles NullCircles;

//---------------------------------------------------------------------------

#endif
