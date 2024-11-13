//---------------------------------------------------------------------------
#pragma hdrstop
//#include "stdafx.h"
#include "MkSphere.hpp"

#ifdef __BCPLUSPLUS__
//#include <vcl.h>
#endif

MkSphere NullSphere(0, 0, 0, 0);
MkSpheres NullSpheres;
//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif
//---------------------------------------------------------------------------
MkSphere::MkSphere(double cx, double cy, double cz, double radius)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    CalcVolume();
    needUpdate = true;
    SurfDivision = 4;
    className = "MkSphere";
}

MkSphere::MkSphere(MkPoint cp, double radius)
{
    FCP = cp;
    FRadius = radius;
    CalcVolume();
    needUpdate = true;
    SurfDivision = 4;
    className = "MkSphere";
}

#ifdef __BCPLUSPLUS__
MkSphere::MkSphere(double cx, double cy, double cz, double radius, TColor C)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    Color = C;
    CalcVolume();
    needUpdate = true;
    SurfDivision = 4;
}

MkSphere::MkSphere(MkPoint cp, double radius, TColor C)
{
    FCP = cp;
    FRadius = radius;
    Color = C;
    CalcVolume();
    needUpdate = true;
    SurfDivision = 4;
}
#endif

MkSphere::MkSphere()
{
    FCP.X = 0;
    FCP.Y = 0;
    FCP.Z = 0;
#ifdef __BCPLUSPLUS__
    Color = clBlack;
#endif
    FRadius = 0;
    needUpdate = true;
    SurfDivision = 4;
    className = "MkSphere";
}

void MkSphere::SetSphere(double cx, double cy, double cz, double radius)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    CalcVolume();
    needUpdate = true;
}

void MkSphere::SetSphere(MkPoint cp, double radius)
{
    FCP = cp;
    FRadius = radius;
    CalcVolume();
    needUpdate = true;
}

#ifdef __BCPLUSPLUS__
void MkSphere::SetSphere(double cx, double cy, double cz, double radius, TColor C)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    Color = C;
    CalcVolume();
    needUpdate = true;
}

void MkSphere::SetSphere(MkPoint cp, double radius, TColor C)
{
    FCP = cp;
    FRadius = radius;
    Color = C;
    CalcVolume();
    needUpdate = true;
}
#endif

void MkSphere::CalcVolume()
{
    FSphereVolume = 4 * M_PI * FRadius * FRadius * FRadius / 3;
}

void MkSphere::SetCenter(double cx, double cy, double cz)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;
    needUpdate = true;
}

void MkSphere::SetCenter(MkPoint cp)
{
    FCP = cp;
    needUpdate = true;
}

void MkSphere::SetRadius(double radius)
{
    FRadius = radius;
    CalcVolume();
    needUpdate = true;
}

bool MkSphere::IsInSurface(MkPoint &pnt, double thick)
{
    return FRadius - thick < CalDist(FCP, pnt) && CalDist(FCP, pnt) < FRadius + thick;
}

bool MkSphere::IsInSpace(MkPoint &pnt)
{
    return FRadius > CalDist(FCP, pnt);
}

MkPoint MkSphere::GetCenter()
{
    return FCP;
}

double MkSphere::GetRadius()
{
    return FRadius;
}

double MkSphere::GetVolume()
{
    return FSphereVolume;
}

MkPoint &MkSphere::operator[](int i)
{
    if (i == 0)
        return FCP;
    else
        return NullPoint;
}

MkSphere &MkSphere::operator=(MkSphere &rc)
{
    MkShape::operator=((MkShape &)rc);
    FCP = rc.FCP;
    FRadius = rc.FRadius;
#ifdef __BCPLUSPLUS__
    Color = rc.Color;
#endif
    CalcVolume();
    return (*this);
}

bool MkSphere::operator==(MkSphere &rs)
{
    if (FCP != rs.FCP || FRadius != rs.FRadius || FSphereVolume != rs.FSphereVolume && FRealPoints != rs.FRealPoints)
        return false;
#ifdef __BCPLUSPLUS__
    if (Color != rs.Color)
        return false;
#endif
    return true;
}

bool MkSphere::operator!=(MkSphere &rs)
{
    return !operator==(rs);
}

bool MkSphere::operator&&(MkLine &rl)
{
    double A, B;
    double l, m, n;
    double len;
    MkPoint sp, ep;
    sp = rl[0];
    ep = rl[1];

    l = ep.X - sp.X;
    m = ep.Y - sp.Y;
    n = ep.Z - sp.Z;

    len = sqrt(l * l + m * m + n * n);

    l = l / len;
    m = m / len;
    n = n / len;

    A = l * (sp.X - FCP.X) + m * (sp.Y - FCP.Y) + n * (sp.Z - FCP.Z);
    B = (sp.X - FCP.X) * (sp.X - FCP.X) + (sp.Y - FCP.Y) * (sp.Y - FCP.Y) + (sp.Z - FCP.Z) * (sp.Z - FCP.Z) - FRadius * FRadius;

    if ((A * A - B) <= 0)
        return false;

    double t[2];
    t[0] = -A - sqrt(A * A - B);
    t[1] = -A + sqrt(A * A - B);

    if (!rl.GetFiniteness() && (A * A - B) > 0)
        return true;
    else
    {
        if (t[0] < 0 && t[1] < 0)
            return false;
        else if (t[0] > 1 && t[1] > 1)
            return false;
        else
            return true;
    }
}

MkPoints MkSphere::operator&(MkLine &rl)
{
    double A, B;
    double l, m, n;
    double len;
    MkPoint sp, ep;
    MkPoints rps;
    sp = rl[0];
    ep = rl[1];

    l = ep.X - sp.X;
    m = ep.Y - sp.Y;
    n = ep.Z - sp.Z;

    len = sqrt(l * l + m * m + n * n);

    l = l / len;
    m = m / len;
    n = n / len;

    A = l * (sp.X - FCP.X) + m * (sp.Y - FCP.Y) + n * (sp.Z - FCP.Z);
    B = (sp.X - FCP.X) * (sp.X - FCP.X) + (sp.Y - FCP.Y) * (sp.Y - FCP.Y) + (sp.Z - FCP.Z) * (sp.Z - FCP.Z) - FRadius * FRadius;

    if ((A * A - B) <= 0)
        return false;

    double t[2];
    t[0] = -A - sqrt(A * A - B);
    t[1] = -A + sqrt(A * A - B);

    rps.Initialize(2);
    rps[0].SetPoint(sp.X + l * t[0], sp.Y + m * t[0], sp.Z + n * t[0]);
    rps[1].SetPoint(sp.X + l * t[1], sp.Y + m * t[1], sp.Z + n * t[1]);

    if (!rl.GetFiniteness() && (A * A - B) > 0)
        return rps;
    else
    {
        if (t[0] < 0 && t[1] < 0)
            return NullPoints;
        else if (t[0] > 1 && t[1] > 1)
            return NullPoints;
        else
            return rps;
    }
}

bool MkSphere::IsIntersect(MkLine &rl)
{
    return *this && rl;
}

MkPoints MkSphere::CalcIntPnts(MkLine &rl)
{
    return *this & rl;
}

void MkSphere::GetIntParam(MkLine &rl, double &t1, double &t2)
{
    double A, B;
    double l, m, n;
    double len;
    MkPoint sp, ep;
    sp = rl[0];
    ep = rl[1];

    l = ep.X - sp.X;
    m = ep.Y - sp.Y;
    n = ep.Z - sp.Z;

    len = sqrt(l * l + m * m + n * n);

    l = l / len;
    m = m / len;
    n = n / len;

    A = l * (sp.X - FCP.X) + m * (sp.Y - FCP.Y) + n * (sp.Z - FCP.Z);
    B = (sp.X - FCP.X) * (sp.X - FCP.X) + (sp.Y - FCP.Y) * (sp.Y - FCP.Y) + (sp.Z - FCP.Z) * (sp.Z - FCP.Z) - FRadius * FRadius;

    if ((A * A - B) <= 0)
    {
        t1 = 0;
        t2 = 0;
        return;
    }
    t1 = -A - sqrt(A * A - B);
    t2 = -A + sqrt(A * A - B);

    double x1, y1, z1;
    x1 = sp.X + l * t1;
    y1 = sp.Y + m * t1;
    z1 = sp.Z + n * t1;
}

bool MkSphere::UpdateSurf()
{
    MkPoint pnt[6];
    Surf.Initialize(8);
    pnt[0].SetPoint(0, 0, FRadius);
    pnt[1].SetPoint(0, 0, -FRadius);
    pnt[2].SetPoint(FRadius, 0, 0);
    pnt[3].SetPoint(-FRadius, 0, 0);
    pnt[4].SetPoint(0, FRadius, 0);
    pnt[5].SetPoint(0, -FRadius, 0);

    for (int i = 0; i < 6; i++)
        pnt[i] += FCP;

    Surf[0].Reset(pnt[0], pnt[3], pnt[5]);
    Surf[1].Reset(pnt[0], pnt[5], pnt[2]);
    Surf[2].Reset(pnt[0], pnt[2], pnt[4]);
    Surf[3].Reset(pnt[0], pnt[4], pnt[3]);
    Surf[4].Reset(pnt[1], pnt[5], pnt[3]);
    Surf[5].Reset(pnt[0], pnt[2], pnt[5]);
    Surf[6].Reset(pnt[0], pnt[4], pnt[2]);
    Surf[7].Reset(pnt[0], pnt[3], pnt[4]);

    needUpdate = false;
    return true;
}

#ifdef __BCPLUSPLUS__
void MkSphere::Draw(TObject *Sender)
{
    /*  TColor C;
  TPenStyle PS;

  if (String(Sender->ClassName()) == String("MkPaintBox")) {
    MkPaintBox *pb=(MkPaintBox*)Sender;
    C = pb->Canvas->Pen->Color;
    PS = pb->Canvas->Pen->Style;
    pb->Canvas->Pen->Color = Color;

    if(needUpdate) UpdateSurf();
    Surf.Draw(pb);

    pb->Canvas->Pen->Color = C;
    pb->Canvas->Pen->Style = PS;
  }*/
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkSphere::Draw(MkPaint *pb)
{
    if (needUpdate)
        UpdateSurf();
    Surf.Draw(pb);
}
#endif
//---------------------------------------------------------------------------
// MkSpheres::MkSpheres(int size)
// {
//     //  try {
//     if (size <= 0)
//     {
//         MkDebug("::MkSpheres - MkSpheres(int size)");
//         return;
//     }

//     FSize = size;
//     FSphere = new MkSphere[FSize];
//     //  }
//     //  catch () {
//     //    MkDebug(AnsiString(E.ClassName())+ E.Message);
//     //  }
// }

// void MkSpheres::Initialize(int size)
// {
//     Clear();
//     FSize = size;
//     if (FSize == 0)
//     {
//         FSphere = NULL;
//         return;
//     }
//     FSphere = new MkSphere[FSize];
// }

// void MkSpheres::Clear()
// {
//     FSize = 0;
//     if (FSphere)
//         delete[] FSphere;
//     FSphere = NULL;
// }

// MkSphere &MkSpheres::operator[](int i)
// {
//     if (i >= 0 && i < FSize)
//         return FSphere[i];
//     else
//         return NullSphere;
// }

// MkSpheres &MkSpheres::operator=(MkSpheres &spheres)
// {
//     int i;
//     FSize = spheres.FSize;
//     this->FSphere = new MkSphere[FSize];

//     for (i = 0; i < FSize; i++)
//         this->FSphere[i] = spheres.FSphere[i];

//     return *this;
// }

// #ifdef __BCPLUSPLUS__
// void MkSpheres::Draw(TObject *Sender)
// {
//     for (int i = 0; i < FSize; i++)
//         FSphere[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkSpheres::Draw(MkPaint *pb)
// {
//     for (int i = 0; i < FSize; i++)
//         FSphere[i].Draw(pb);
// }
// #endif

//---------------------------------------------------------------------------
