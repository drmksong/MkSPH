//---------------------------------------------------------------------------
#ifndef MkLineHPP
#define MkLineHPP

#ifdef __BCPLUSPLUS__
#ifdef _MSC_VER
#undef _MSC_VER
#endif
#endif

#include <string>
#include "math.h"
#include "MkContainer.hpp"
#include "MkShape.hpp"
#include "MkPoint.hpp"

#if defined(_MSC_VER) && defined(_WINDOWS_)
class MkPaint;
#endif

class MkArc;
//---------------------------------------------------------------------------
class MkLine : public MkShape
{
private:
    MkPoint StartPoint, EndPoint;
    double Theta, Length;
    double A, B, C; // A*X + B*Y = C  (C = 0 or 1)
    double L, M, N; // direction cosine
    bool isSelected;
    bool isFinite;

    void CalLength();
    void CalTheta();
    void CalCoeff();
    double Tol;

public:
    MkLine();
    MkLine(MkPoint sp, MkPoint ep);
    MkLine(double sx, double sy, double ex, double ey);

#ifdef __BCPLUSPLUS__
    MkLine(MkPoint sp, MkPoint ep, TColor C);
    MkLine(double sx, double sy, double ex, double ey, TColor C);
    void SetLine(MkPoint sp, MkPoint ep, TColor C);
    void SetLine(double sx, double sy, double ex, double ey, TColor C);
#endif

#ifdef __BCPLUSPLUS__
    AnsiString ClassName()
    {
        return AnsiString("MkLine");
    }
#else
    std::string ClassName()
    {
        return className;
    }
#endif

    void SetLine(MkPoint sp, MkPoint ep);
    void SetLine(double sx, double sy, double ex, double ey);
    void SetTol(double tol) { Tol = tol; }
    double GetLength();
    double GetTheta();
    bool GetFiniteness() { return isFinite; }
    void SetFiniteness(bool fin) { isFinite = fin; }
    void AdjustTheta()
    {
        if (Theta > 360)
            Theta -= 360;
        else if (Theta < 0)
            Theta += 360;
    }
    MkPoint &GetMiddlePoint()
    {
        static MkPoint mp;
        mp.SetPoint(
            ((*this)[0].X + (*this)[1].X) / 2,
            ((*this)[0].Y + (*this)[1].Y) / 2,
            ((*this)[0].Z + (*this)[1].Z) / 2);
        return mp;
    };
    MkPoint &GetDivision(double f);

    void Extend(double f);
    double DeltaX() { return EndPoint.X - StartPoint.X; };
    double DeltaY() { return EndPoint.Y - StartPoint.Y; };
    double DeltaZ() { return EndPoint.Z - StartPoint.Z; };
    double NormDeltaX() { return fabs(StartPoint.X) < 1.0e-6 ? 0 : ((EndPoint.X - StartPoint.X) / StartPoint.X); };
    double NormDeltaY() { return fabs(StartPoint.Y) < 1.0e-6 ? 0 : ((EndPoint.Y - StartPoint.Y) / StartPoint.Y); };
    bool IsIntersect(MkLine &realline);
    bool IsIn(MkPoint rp);
    bool IsIn(double x, double y);
    bool IsInLine(MkPoint rp);
    bool IsInLine(double x, double y);
    bool IsInSamePlane(MkLine rl);
    bool GetIntParam(MkLine &rl, double &t);
    bool GetIntParam(MkPoint p, double &t);
    MkPoint &GetIntPoint(MkLine &rl);
    double GetArea() { return 0; };
    double GetA() { return A; }
    double GetB() { return B; }
    double GetC() { return C; }
    double GetL() { return L; }
    double GetM() { return M; }
    double GetN() { return N; }

    MkVector<double> &GetVector()
    {
        static MkVector<double> vec(3);
        vec[0] = EndPoint.X - StartPoint.X;
        vec[1] = EndPoint.Y - StartPoint.Y;
        vec[2] = EndPoint.Z - StartPoint.Z;
        vec.Normalize();
        return vec;
    }

    void Select() { isSelected = true; }
    void Unselect() { isSelected = false; }
    bool GetSelected() { return isSelected; }

    double CalDist(MkPoint rp) { return fabs(A * rp.X + B * rp.Y - C) / sqrt(A * A + B * B); }
    double CalDist(double x, double y) { return fabs(A * x + B * y - C) / sqrt(A * A + B * B); }
    double CalDist3D(MkPoint rp);
    double CalDist3D(double x, double y, double z);
    MkPoint &GetNearestPnt(MkPoint rp);
    MkPoint &GetNearestPnt(double x, double y, double z);

    MkLine Translate(MkPoint rp);
    MkLine Translate(double x, double y, double z);
    MkLine Rotate(double alpha, double beta, double gamma);
    MkLine RotateInX(double ang);
    MkLine RotateInY(double ang);
    MkLine RotateInZ(double ang);
    MkLine RotateInA(double ang, double l, double m, double n);

    MkLine Scale(double sx, double sy, double sz);

    void Out();
    void Clear();
    bool operator&&(MkLine rline1);
    bool operator==(MkLine rline1);
    bool operator!=(MkLine rline1);
    MkPoint &operator&(MkLine &);
    MkPoint &operator[](int);
    double operator+(MkLine); //dot product
    double operator*(MkLine);
    double operator*(MkPoint rp); //cross product
    MkLine &operator*(MkMatrix4<double> &rm);
    double operator+=(MkPoint);
    double operator-=(MkPoint);

#if !defined(_MSC_VER) && !defined(_WINDOWS_) || defined(__BCPLUSPLUS__)
    MkLine &operator=(const MkLine &);
#endif
    MkLine &operator=(MkLine &);
    MkLine operator!();

    void Draw(void);
#ifdef __BCPLUSPLUS__
    void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
    void Draw(MkPaint *);
#endif
};

// class MkLines
// {
// protected:
//     MkLine *FRealLine;
//     int FSize;
// #ifdef __BCPLUSPLUS__
//     TColor Color;
// #endif
// #if defined(_MSC_VER) && defined(_WINDOWS_)
//     MkColor PColor, BColor;
// #endif
// public:
//     MkLines(int size, MkLine *rl);
//     MkLines(int FSize);
//     MkLines()
//     {
//         FSize = 0;
//         FRealLine = NULL;
//     }
//     ~MkLines();
//     void Initialize(int size);
//     void Initialize(int size, MkLine *rl);
//     void Grow(int sz);
//     void Add(MkLine &l)
//     {
//         Grow(1);
//         FRealLine[FSize - 1] = l;
//     }
//     void DeleteSelected();
//     int GetSize() { return FSize; };
//     int GetNumber() { return FSize; };
//     MkLine *GetLine() { return FRealLine; }
//     bool Clear();
// #ifdef __BCPLUSPLUS__
//     TColor GetColor()
//     {
//         return Color;
//     };
//     void SetColor(TColor c) { Color = c; }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//     MkColor GetColor()
//     {
//         return PColor;
//     };
//     void SetColor(MkColor c) { PColor = c; }
// #endif

//     virtual MkLine &operator[](int);

//     MkLines &RotateInX(double ang);
//     MkLines &RotateInY(double ang);
//     MkLines &RotateInZ(double ang);
//     MkLines &RotateInA(double ang, double l, double m, double n);

//     MkLines &operator*(MkMatrix4<double> &rm);
//     MkLines &operator=(MkLines &);
//     bool operator==(MkLines &);
//     void Draw(void);

// #ifdef __BCPLUSPLUS__
//     void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//     void Draw(MkPaint *);
// #endif
// };

typedef MkContainer<MkLine> MkLines;

extern MkLine NullLine;
extern MkLines NullLines;
//---------------------------------------------------------------------------
#endif
