//---------------------------------------------------------------------------
#ifndef MkPointDataH
#define MkPointDataH
#endif

#include <stdio.h>
#include <math.h>
#include "MkMisc.hpp"
#include "MkPoint.hpp"
#include "MkArray.hpp"
#include "MkCube.hpp"
#include "MkFault.hpp"
#include "MkPlane.hpp"
#include "MkMatrix.hpp"
#include "MkColor.hpp"

struct MkDataPoint : MkPoint
{
public:
    double Data;

    MkDataPoint();
    MkDataPoint(double x, double y);
    MkDataPoint(double x, double y, double z);

    MkDataPoint &operator=(MkDataPoint rp);
    MkDataPoint &operator=(MkPoint rp);
    MkDataPoint &operator=(double data)
    {
        Data = data;
        return *this;
    };
    void operator+=(MkDataPoint &rp)
    {
        X += rp.X;
        Y += rp.Y;
        Z += rp.Z;
    }
    void operator+=(double data) { Data += data; }
    MkDataPoint operator*(MkMatrix4<double> &rm);
    friend MkDataPoint operator+(MkDataPoint a, MkDataPoint b)
    {
        MkDataPoint c;
        c.X = a.X + b.X;
        c.Y = a.Y + b.Y;
        return c;
    }
    bool operator==(MkDataPoint);
    bool operator!=(MkDataPoint);
};

class MkDataPoints
{
protected:
    boost::shared_array<MkDataPoint> FPoint;
    int FSize;
    MkColor Color;

public:
    MkDataPoints(int size, MkDataPoint *rps);
    MkDataPoints(int FSize);
    MkDataPoints()
    {
        FSize = 0;
    }
    ~MkDataPoints();
    virtual void Initialize(int size);
    virtual void Initialize(int size, MkDataPoint *);
    int GetSize() { return FSize; };
    int GetNumber() { return FSize; };
    bool Add(MkDataPoint point);
    //    void Delete(MkDataPoint point);
    void Grow(int Delta);
    //    Shrink(int Delta);
    bool Clear();
    MkColor GetColor() { return Color; };
    void SetColor(MkColor c) { Color = c; }
    virtual MkDataPoint &operator[](int);

    MkDataPoints &operator*(MkMatrix4<double> &rm);

    MkDataPoints &Translate(MkPoint rp);
    MkDataPoints &Translate(MkDataPoint rp);
    MkDataPoints &Translate(double x, double y, double z);
    MkDataPoints &Rotate(double alpha, double beta, double gamma);
    MkDataPoints &Scale(double sx, double sy, double sz);

    MkDataPoints &operator=(MkDataPoints &points);
    bool operator==(MkDataPoints &points);
    //virtual void Draw(TObject *);

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

extern MkDataPoint NullDataPoint;
extern MkDataPoints NullDataPoints;
//---------------------------------------------------------------------------
