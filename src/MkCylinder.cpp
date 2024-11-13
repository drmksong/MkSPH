//---------------------------------------------------------------------------
#include "MkCylinder.hpp"

MkCylinder NullCylinder(0, 0, 0, 0);
MkCylinders NullCylinders;
//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif
//---------------------------------------------------------------------------
MkCylinder::MkCylinder()
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
    className = "MkCylinder";
}

MkCylinder::MkCylinder(float cx, float cy, float cz, float radius)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    needUpdate = true;
    SurfDivision = 4;
    className = "MkCylinder";
}

MkCylinder::MkCylinder(MkPoint cp, float radius)
{
    FCP = cp;
    FRadius = radius;
    needUpdate = true;
    SurfDivision = 4;
    className = "MkCylinder";
}
#ifdef __BCPLUSPLUS__
MkCylinder::MkCylinder(float cx, float cy, float cz, float radius, TColor C)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    Color = C;
    needUpdate = true;
    SurfDivision = 4;
}

MkCylinder::MkCylinder(MkPoint cp, float radius, TColor C)
{
    FCP = cp;
    FRadius = radius;
    Color = C;
    needUpdate = true;
    SurfDivision = 4;
}
#endif

void MkCylinder::SetCylinder(float cx, float cy, float cz, float radius)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    needUpdate = true;
    className = "MkCylinder";
}

void MkCylinder::SetCylinder(MkPoint cp, float radius)
{
    FCP = cp;
    FRadius = radius;
    needUpdate = true;
    className = "MkCylinder";
}
#ifdef __BCPLUSPLUS__
void MkCylinder::SetCylinder(float cx, float cy, float cz, float radius, TColor C)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;

    FRadius = radius;
    Color = C;
    needUpdate = true;
}

void MkCylinder::SetCylinder(MkPoint cp, float radius, TColor C)
{
    FCP = cp;
    FRadius = radius;
    Color = C;
    needUpdate = true;
}
#endif

void MkCylinder::SetCenter(float cx, float cy, float cz)
{
    FCP.X = cx;
    FCP.Y = cy;
    FCP.Z = cz;
    needUpdate = true;
}

void MkCylinder::SetCenter(MkPoint cp)
{
    FCP = cp;
    needUpdate = true;
}

void MkCylinder::SetRadius(float radius)
{
    FRadius = radius;
    needUpdate = true;
}

void MkCylinder::SetOrient(MkPoint orient)
{
    float len;

    Fl = orient.X;
    Fm = orient.Y;
    Fn = orient.Z;

    len = sqrt(Fl * Fl + Fm * Fm + Fn * Fn);

    Fl /= len;
    Fm /= len;
    Fn /= len;
    needUpdate = true;
}

void MkCylinder::SetOrient(float l, float m, float n)
{
    float len;

    Fl = l;
    Fm = m;
    Fn = n;

    len = sqrt(Fl * Fl + Fm * Fm + Fn * Fn);

    Fl /= len;
    Fm /= len;
    Fn /= len;
    needUpdate = true;
}

bool MkCylinder::IsInSurface(MkPoint &pnt, float thick)
{
    float dist = GetDist(pnt);
    return (FRadius - thick < dist && dist < FRadius + thick);
}

bool MkCylinder::IsInSpace(MkPoint &pnt)
{
    return GetDist(pnt) < FRadius;
}

MkPoint MkCylinder::GetCenter()
{
    return FCP;
}

float MkCylinder::GetRadius()
{
    return FRadius;
}

float MkCylinder::GetDist(MkPoint &pnt)
{
    float len;
    float l, m, n; //cross vector
    float x, y, z; //point
    len = sqrt(Fl * Fl + Fm * Fm + Fn * Fn);
    assert(len - 1.0 < 0.001);
    x = pnt.X - FCP.X;
    y = pnt.Y - FCP.Y;
    z = pnt.Z - FCP.Z;

    l = Fm * z - y * Fn;
    m = Fn * x - z * Fl;
    n = Fl * y - x * Fm;

    return sqrt(l * l + m * m + n * n);
}

MkPoint &MkCylinder::operator[](int i)
{
    if (i == 0)
        return FCP;
    else
        return NullPoint;
}

MkCylinder &MkCylinder::operator=(MkCylinder &rc)
{
    MkShape::operator=((MkShape &)rc);
    FCP = rc.FCP;
    FRadius = rc.FRadius;
#ifdef __BCPLUSPLUS__
    Color = rc.Color;
#endif
    return (*this);
}

bool MkCylinder::operator&&(MkLine &rl)
{
    return false;
}

MkPoints &MkCylinder::operator&(MkLine &rl)
{
    //calc FPoints;
    return FPoints;
}

bool MkCylinder::operator==(MkCylinder &rs)
{
    return FCP == rs.FCP && fabs(FRadius - rs.FRadius) < EPS && fabs(FLength - rs.FLength) < EPS && fabs(Fl - rs.Fl) < EPS &&
           fabs(Fm - rs.Fm) < EPS && fabs(Fn - rs.Fn) < EPS && fabs(Psi - rs.Psi) < EPS && fabs(Theta - rs.Theta) && FPoints == rs.FPoints;
}

bool MkCylinder::operator!=(MkCylinder &rs)
{
    return !operator==(rs);
}

bool MkCylinder::IsIntersect(MkLine &rl)
{
    return *this && rl;
}

MkPoints &MkCylinder::CalcIntPnts(MkLine &rl)
{
    return *this & rl;
}

void MkCylinder::GetIntParam(MkLine &rl, float &t1, float &t2)
{
}

void MkCylinder::RotateSpace(MkPoint &rp)
{
    rp -= FCP;
    rp.Rotate(0, Psi, 0);   //rotate in y-axis
    rp.Rotate(Theta, 0, 0); //rotate in z-axis --- check it out
}

void MkCylinder::RotateSpace(MkLine &rl)
{
    rl -= FCP;
    rl.Rotate(0, Psi, 0);
    rl.Rotate(Theta, 0, 0);
}

void MkCylinder::RotateSpace(MkPlane &rp)
{
}

void MkCylinder::RotateSpace(MkJointPlane &jp)
{
}

void MkCylinder::RotateSpace(MkPennyJoint &pj)
{
}

void MkCylinder::UnRotateSpace(MkPoint &rp)
{
}

void MkCylinder::UnRotateSpace(MkLine &rl)
{
}

void MkCylinder::UnRotateSpace(MkPlane &rp)
{
}

void MkCylinder::UnRotateSpace(MkJointPlane &jp)
{
}

void MkCylinder::UnRotateSpace(MkPennyJoint &pj)
{
}

bool MkCylinder::UpdateSurf()
{
    int ldiv, cdiv;
    float len;
    cdiv = int(4 * pow(2, SurfDivision));
    len = 2 * 3.141592 * FRadius / cdiv;
    ldiv = int(FLength / len);
    if (ldiv < 2)
        ldiv = 2;

    needUpdate = false;
    return true;
}

#ifdef __BCPLUSPLUS__
void MkCylinder::Draw(TObject *Sender)
{
    if (needUpdate)
        UpdateSurf();
    Surf.Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkCylinder::Draw(MkPaint *pb)
{
    if (needUpdate)
        UpdateSurf();
    Surf.Draw(pb);
}
#endif

//---------------------------------------------------------------------------
// MkCylinders::MkCylinders(int size)
// {
//     //  try {
//     if (size <= 0)
//     {
//         MkDebug("::MkCylinders - MkCylinders(int size)");
//         return;
//     }

//     FSize = size;
//     FCylinder = new MkCylinder[FSize];
//     //  }
//     //  catch () {
//     //    MkDebug(AnsiString(E.ClassName())+ E.Message);
//     //  }
// }

// void MkCylinders::Initialize(int size)
// {
//     Clear();
//     FSize = size;
//     if (FSize == 0)
//     {
//         FCylinder = NULL;
//         return;
//     }
//     FCylinder = new MkCylinder[FSize];
// }

// void MkCylinders::Clear()
// {
//     FSize = 0;
//     if (FCylinder)
//         delete[] FCylinder;
//     FCylinder = NULL;
// }

// MkCylinder &MkCylinders::operator[](int i)
// {
//     if (i >= 0 && i < FSize)
//         return FCylinder[i];
//     else
//         return NullCylinder;
// }

// MkCylinders &MkCylinders::operator=(MkCylinders &cylinders)
// {
//     int i;
//     FSize = cylinders.FSize;
//     this->FCylinder = new MkCylinder[FSize];

//     for (i = 0; i < FSize; i++)
//         this->FCylinder[i] = cylinders.FCylinder[i];

//     return *this;
// }

// #ifdef __BCPLUSPLUS__
// void MkCylinders::Draw(TObject *Sender)
// {
//     for (int i = 0; i < FSize; i++)
//         FCylinder[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkCylinders::Draw(MkPaint *pb)
// {
//     for (int i = 0; i < FSize; i++)
//         FCylinder[i].Draw(pb);
// }
// #endif

//---------------------------------------------------------------------------
