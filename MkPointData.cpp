//---------------------------------------------------------------------------
// This module is general purposed simple graphic class to store, draw,
// manipulate object. It is well suited to VCL component, but not restricted.
// It forms the base for the higher level class, such as tunnel component.
//
// Copyright (c) 1999-2000 Myung Kyu Song, ESCO Consultant Co., Ltd.

#include "MkPointData.hpp"
#define ShowMessage printf

MkDataPoint NullDataPoint(0, 0, 0);
MkDataPoints NullDataPoints(0);

MkDataPoint::MkDataPoint()
    : MkPoint()
{
    Data = 0;
};

MkDataPoint::MkDataPoint(double x, double y)
    : MkPoint(x, y)
{
    Data = 0;
};

MkDataPoint::MkDataPoint(double x, double y, double z)
    : MkPoint(x, y, z)
{
    Data = 0;
};

MkDataPoint &MkDataPoint::operator=(MkPoint rp)
{
    X = rp.X;
    Y = rp.Y;
    Z = rp.Z;
    Data = 0;
    return (*this);
}

MkDataPoint &MkDataPoint::operator=(MkDataPoint rp)
{
    X = rp.X;
    Y = rp.Y;
    Z = rp.Z;
    Data = rp.Data;
    return (*this);
}

bool MkDataPoint::operator==(MkDataPoint rp)
{
    return (fabs(X - rp.X) < FTOL) &&
           (fabs(Y - rp.Y) < FTOL) &&
           (fabs(Z - rp.Z) < FTOL) &&
           (fabs(Data - rp.Data) < FTOL);
}

bool MkDataPoint::operator!=(MkDataPoint rp)
{
    return (fabs(X - rp.X) > FTOL) ||
           (fabs(Y - rp.Y) > FTOL) ||
           (fabs(Z - rp.Z) > FTOL) ||
           (fabs(Data - rp.Data) > FTOL);
}

MkDataPoint MkDataPoint::operator*(MkMatrix4<double> &rm)
{
    MkDataPoint rp;

    rp.X = X * rm(0, 0) + Y * rm(0, 1) + Z * rm(0, 2);
    rp.Y = X * rm(1, 0) + Y * rm(1, 1) + Z * rm(1, 2);
    rp.Z = X * rm(2, 0) + Y * rm(2, 1) + Z * rm(2, 2);
    return rp;
}

MkDataPoints::MkDataPoints(int size, MkDataPoint *rps)
{

    if (size < 0)
    {
        MkDebug("::MkDataPoints - MkDataPoints(int size)");
        throw Size(std::string("MkDataPoints::Constructor throw negative size exception"), size);
    }
    if (rps == NULL)
    {
        MkDebug("::MkDataPoints - MkDataPoints(int size)");
        throw Size(std::string("MkDataPoints::Constructor throw null pointer argument"), size);
    }

    FSize = size;
    if (FSize == 0)
    {
        FPoint.reset();
        return;
    }
    try
    {
        FPoint.reset(new MkDataPoint[FSize]);
    }
    catch (std::bad_alloc &a)
    {
        MkDebug("MkDataPoints::MkDataPoints bad_alloc thrown!!!\n");
        throw Alloc(a.what());
    }

    for (int i = 0; i < FSize; i++)
        FPoint[i] = rps[i];
}

MkDataPoints::MkDataPoints(int size)
{
    if (size < 0)
    {
        MkDebug("::MkDataPoints - MkDataPoints(int size)");
        throw Size(std::string("MkDataPoints::Constructor throw negative size exception"), size);
    }

    FSize = size;
    if (FSize == 0)
    {
        FPoint.reset();
        return;
    }

    try
    {
        FPoint.reset(new MkDataPoint[FSize]);
    }
    catch (std::bad_alloc &a)
    {
        MkDebug("MkDataPoints::MkDataPoints bad_alloc thrown!!!\n");
        throw Alloc(a.what());
    }
}

MkDataPoints::~MkDataPoints()
{
    FSize = 0;
    FPoint.reset();
}

void MkDataPoints::Initialize(int size)
{

    if (size < 0)
    {
        MkDebug("::MkDataPoints - Initialize(int size)");
        throw Size(std::string("MkDataPoints::Initialize throw negative size exception"), size);
    }
    if (FSize == size)
        return;

    FSize = size;
    if (FSize == 0)
    {
        FPoint.reset();
        return;
    }
    try
    {
        FPoint.reset(new MkDataPoint[FSize]);
    }
    catch (std::bad_alloc &a)
    {
        MkDebug("MkDataPoints::Initialize bad_alloc thrown!!!\n");
        throw Alloc(a.what());
    }
}

void MkDataPoints::Initialize(int size, MkDataPoint *rps)
{

    if (size < 0)
    {
        MkDebug("::MkDataPoints - Initialize(int size)");
        throw Size(std::string("MkDataPoints::Initialize throw negative size exception"), size);
    }
    if (rps == NULL)
    {
        MkDebug("::MkDataPoints - Initialize(int size)");
        throw Size(std::string("MkDataPoints::Initialize throw null pointer argument"), size);
    }
    if (FSize == size)
        return;
    FSize = size;
    if (FSize == 0)
    {
        FPoint.reset();
        return;
    }
    try
    {
        FPoint.reset(new MkDataPoint[FSize]);
    }

    catch (std::bad_alloc &a)
    {
        MkDebug("MkDataPoints::Initialize bad_alloc thrown!!!\n");
        throw Alloc(a.what());
    }

    for (int i = 0; i < FSize; i++)
        (*this)[i] = rps[i];
}

void MkDataPoints::Grow(int Delta)
{
    int i;
    boost::shared_array<MkDataPoint> rp;
    try
    {
        rp.reset(new MkDataPoint[FSize + Delta]);
    }
    catch (std::bad_alloc &a)
    {
        MkDebug("MkDataPoints::Grow bad_alloc thrown!!!\n");
        throw Alloc(a.what());
    }

    for (i = 0; i < FSize; i++)
        rp[i] = FPoint[i];
    for (i = FSize; i < FSize + Delta; i++)
        rp[i] = NullDataPoint;
    FPoint = rp;
    FSize = FSize + Delta;
}

bool MkDataPoints::Add(MkDataPoint point)
{
    Grow(1);
    FPoint[FSize - 1] = point;
    return true;
}

bool MkDataPoints::Clear()
{

    FSize = 0;
    FPoint.reset();

    return true;
}

MkDataPoint &MkDataPoints::operator[](int i)
{
    if (FSize == 0)
        return NullDataPoint;
    else if (i >= 0 && i < FSize)
        return FPoint[i];
    else
        return NullDataPoint;
}

MkDataPoints &MkDataPoints::operator=(MkDataPoints &points)
{
    int i;

    Clear();
    FSize = points.FSize;
    if (FSize == 0)
    {
        FPoint.reset();
        return *this;
    }

    try
    {
        FPoint.reset(new MkDataPoint[FSize]);
    }
    catch (std::bad_alloc &a)
    {
        MkDebug("MkDataPoints::operator= bad_alloc thrown!!!\n");
        throw Alloc(a.what());
    }

    for (i = 0; i < FSize; i++)
        this->FPoint[i] = points.FPoint[i];

    return *this;
}

bool MkDataPoints::operator==(MkDataPoints &points)
{
    int i;

    if (FSize != points.FSize)
        return false;
    for (i = 0; i < FSize; i++)
        if (this->FPoint[i] != points.FPoint[i])
            return false;

    return true;
}

MkDataPoints &MkDataPoints::operator*(MkMatrix4<double> &rm)
{
    for (int i = 0; i < FSize; i++)
        this->FPoint[i] = this->FPoint[i] * rm;
    return *this;
}

MkDataPoints &MkDataPoints::Translate(MkPoint rp)
{
    for (int i = 0; i < FSize; i++)
        this->FPoint[i].Translate(rp);
    return *this;
}

MkDataPoints &MkDataPoints::Translate(MkDataPoint rp)
{
    for (int i = 0; i < FSize; i++)
        this->FPoint[i].Translate(rp);
    return *this;
}

MkDataPoints &MkDataPoints::Translate(double x, double y, double z)
{
    for (int i = 0; i < FSize; i++)
        this->FPoint[i].Translate(x, y, z);
    return *this;
}

MkDataPoints &MkDataPoints::Rotate(double alpha, double beta, double gamma)
{
    for (int i = 0; i < FSize; i++)
        this->FPoint[i].Rotate(alpha, beta, gamma);
    return *this;
}

MkDataPoints &MkDataPoints::Scale(double sx, double sy, double sz)
{
    for (int i = 0; i < FSize; i++)
        this->FPoint[i].Scale(sx, sy, sz);
    return *this;
}

/*
void MkDataPoints::Draw(TObject *Sender)
{
    MkColor C;
    double Offset;
    if (FSize == 0) return;
    if (String(Sender->ClassName()) == String("TWorldPaintBox")) {
       TWorldPaintBox *pb=(TWorldPaintBox*)Sender;
       C = pb->Canvas->Pen->Color;
       pb->Canvas->Pen->Color = Color;
       for (int i = 0; i < FSize;i++) {
           Offset = pb->Offset(3);
//           pb->Cube( (*this)[i].X-Offset,(*this)[i].Y-Offset,(*this)[i].Z-Offset,
//                            (*this)[i].X+Offset,(*this)[i].Y+Offset,(*this)[i].Z+Offset);
//           pb->Cube( (*this)[i].X-Offset,(*this)[i].Y-Offset,0,
//                            (*this)[i].X+Offset,(*this)[i].Y+Offset,0);
       }
       pb->Canvas->Pen->Color = C;
    }
}*/
//---------------------------------------------------------------------------
