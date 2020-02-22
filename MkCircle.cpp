//---------------------------------------------------------------------------
// This module is general purposed simple graphic class to store, draw,
// manipulate object. It is well suited to VCL component, but not restricted.
// It forms the base for the higher level class, such as tunnel component.
//
// Copyright (c) 1999 Myung Kyu Song, ESCO Consultant Co., Ltd.
#include "MkCircle.hpp"

MkCircle NullCircle(0);

//---------------------------------------------------------------------------
MkCircle::MkCircle()
{
    FCP.X = 0;
    FCP.Y = 0;
#ifdef __BCPLUSPLUS__
    Color = clBlack;
#endif
    FRadius = 0;
    className = "MkCircle";
}

MkCircle::MkCircle(int)
{
    FCP.X = 0;
    FCP.Y = 0;
    FRadius = 0;
    className = "MkCircle";
}

MkCircle::MkCircle(double cx, double cy, double radius)
{
    FCP.X = cx;
    FCP.Y = cy;
    FRadius = radius;
    CalArea();
    className = "MkCircle";
}

MkCircle::MkCircle(MkPoint cp, double radius)
{
    FCP = cp;
    FRadius = radius;
    CalArea();
    className = "MkCircle";
}

#ifdef __BCPLUSPLUS__
MkCircle::MkCircle(double cx, double cy, double radius, TColor C)
{
    FCP.X = cx;
    FCP.Y = cy;
    FRadius = radius;
    Color = C;
    CalArea();
}

MkCircle::MkCircle(MkPoint cp, double radius, TColor C)
{
    FCP = cp;
    FRadius = radius;
    Color = C;
    CalArea();
}
#endif

void MkCircle::SetCircle(double cx, double cy, double radius)
{
    FCP.X = cx;
    FCP.Y = cy;
    FRadius = radius;
    CalArea();
}

void MkCircle::SetCircle(MkPoint cp, double radius)
{
    FCP = cp;
    FRadius = radius;
    CalArea();
}
#ifdef __BCPLUSPLUS__
void MkCircle::SetCircle(double cx, double cy, double radius, TColor C)
{
    FCP.X = cx;
    FCP.Y = cy;
    FRadius = radius;
    Color = C;
    CalArea();
}

void MkCircle::SetCircle(MkPoint cp, double radius, TColor C)
{
    FCP = cp;
    FRadius = radius;
    Color = C;
    CalArea();
}
#endif

void MkCircle::CalArea()
{
    FCircleArea = M_PI * FRadius * FRadius;
}

void MkCircle::SetCenter(double cx, double cy)
{
    FCP.X = cx;
    FCP.Y = cy;
}

void MkCircle::SetCenter(MkPoint cp)
{
    FCP = cp;
}

void MkCircle::SetRadius(double radius)
{
    FRadius = radius;
    CalArea();
}

bool MkCircle::IsInSurface(MkPoint &pnt, double thick)
{
    double d;
    d = CalDist(FCP, pnt);
    if (FRadius - thick < d && d < FRadius + thick)
        return true;
    else
        return false;
}

bool MkCircle::IsInSpace(MkPoint &pnt)
{
    double d;
    d = CalDist(FCP, pnt);
    if (d < FRadius)
        return true;
    else
        return false;
}

MkPoint &MkCircle::operator[](int i)
{
    if (i == 0)
        return FCP;
    else
        return NullPoint;
}

MkCircle &MkCircle::operator=(MkCircle &rc)
{
    this->MkShape::operator=((MkShape &)rc);
    FCP.X = rc.FCP.X;
    FCP.Y = rc.FCP.Y;
    FRadius = rc.FRadius;
#ifdef __BCPLUSPLUS__
    Color = rc.Color;
#endif
    CalArea();
    return (*this);
}

bool MkCircle::operator&&(MkLine &rl)
{
    double xn1, yn1, xn2, yn2;
    double d;
    d = fabs(rl.GetA() * FCP.X + rl.GetB() * FCP.Y - rl.GetC()) / sqrt(rl.GetA() * rl.GetA() + rl.GetB() * rl.GetB());

    if (d > FRadius)
        return false;
    double r2 = FRadius * FRadius;
    if (rl.GetB() < FTOL)
    {
        double tmp;
        double B2, A2, C2, AC, Cx2, Cy2;
        B2 = rl.GetB() * rl.GetB();
        C2 = rl.GetC() * rl.GetC();
        A2 = rl.GetA() * rl.GetA();
        Cx2 = FCP.X * FCP.X;
        Cy2 = FCP.Y * FCP.Y;

        AC = (-2 * FCP.X + 2 * FCP.Y * rl.GetA() / rl.GetB() - 2 * rl.GetA() * rl.GetC() / B2);

        tmp = (AC * AC - 4 * (1 + A2 / B2) * (Cx2 + Cy2 - r2 + C2 / B2 - 2 * FCP.Y * rl.GetC() / rl.GetB()));
        if (tmp < 0)
            return false;
        tmp = sqrt(tmp);

        xn1 = (-AC + tmp) / (2 * (1 + A2 / B2));
        yn1 = -rl.GetA() / rl.GetB() * xn1 + rl.GetC() / rl.GetB();
        xn2 = (-AC - tmp) / (2 * (1 + A2 / B2));
        yn2 = -rl.GetA() / rl.GetB() * xn2 + rl.GetC() / rl.GetB();
    }
    else
    {
        xn1 = rl.GetC() / rl.GetA();
        xn2 = rl.GetC() / rl.GetA();
        yn1 = FCP.Y + sqrt(r2 - (xn1 - FCP.X) * (xn1 - FCP.X));
        yn2 = FCP.Y - sqrt(r2 - (xn2 - FCP.X) * (xn2 - FCP.X));
    }

    if (rl.IsIn(xn1, yn1))
    {
        if (rl.IsIn(xn2, yn2))
            return true;
        else
            return true;
    }

    else if (rl.IsIn(xn2, yn2))
        return true;
    else
        return false;
}

MkPoints MkCircle::operator&(MkLine &rl)
{
    double xn1, yn1, xn2, yn2;

    double d;
    d = fabs(rl.GetA() * FCP.X + rl.GetB() * FCP.Y - rl.GetC()) / sqrt(rl.GetA() * rl.GetA() + rl.GetB() * rl.GetB());

    if (d > FRadius)
        return NullPoints;
    double r2 = FRadius * FRadius;
    if (rl.GetB() < FTOL)
    {
        double tmp;
        double B2, A2, C2, AC, Cx2, Cy2;
        B2 = rl.GetB() * rl.GetB();
        C2 = rl.GetC() * rl.GetC();
        A2 = rl.GetA() * rl.GetA();
        Cx2 = FCP.X * FCP.X;
        Cy2 = FCP.Y * FCP.Y;

        AC = (-2 * FCP.X + 2 * FCP.Y * rl.GetA() / rl.GetB() - 2 * rl.GetA() * rl.GetC() / B2);

        tmp = (AC * AC - 4 * (1 + A2 / B2) * (Cx2 + Cy2 - r2 + C2 / B2 - 2 * FCP.Y * rl.GetC() / rl.GetB()));
        if (tmp < 0)
            return NullPoints;
        tmp = sqrt(tmp);

        xn1 = (-AC + tmp) / (2 * (1 + A2 / B2));
        yn1 = -rl.GetA() / rl.GetB() * xn1 + rl.GetC() / rl.GetB();
        xn2 = (-AC - tmp) / (2 * (1 + A2 / B2));
        yn2 = -rl.GetA() / rl.GetB() * xn2 + rl.GetC() / rl.GetB();
    }
    else
    {
        xn1 = rl.GetC() / rl.GetA();
        xn2 = rl.GetC() / rl.GetA();
        yn1 = FCP.Y + sqrt(r2 - (xn1 - FCP.X) * (xn1 - FCP.X));
        yn2 = FCP.Y - sqrt(r2 - (xn2 - FCP.X) * (xn2 - FCP.X));
    }

    if (rl.IsIn(xn1, yn1))
    {

        if (rl.IsIn(xn2, yn2))
        {
            MkPoints rps(2);
            rps[0].SetPoint(xn1, yn1);
            rps[1].SetPoint(xn2, yn2);
            return rps;
        }
        else
        {
            MkPoints rps(1);
            rps[0].SetPoint(xn1, yn1);
            return rps;
        }
    }

    else if (rl.IsIn(xn2, yn2))
    {
        MkPoints rps(1);
        rps[0].SetPoint(xn2, yn2);
        return rps;
    }
    else
        return NullPoints;
}

MkPoints &MkCircle::operator&(MkPoint &rp) //�������ϱ�
{
    double a, b, c;
    double xp, yp;
    double x0, y0;
    double A, A2;
    double r, r2;
    double g, g2;

    xp = rp.X - FCP.X;
    yp = rp.Y - FCP.Y;
    x0 = 0;
    y0 = 0;
    r = FRadius;
    r2 = r * r;

    if (MkLine(FCP, rp).GetLength() < FRadius)
        return NullPoints;

    FRealPoints.Initialize(2);
    if (fabs(x0 - xp) < FTOL)
    {
        double dy = yp - y0;
        FRealPoints[0].Y = FRealPoints[1].Y = -r2 / dy + FCP.Y;
        FRealPoints[0].X = x0 + sqrt(r2 - (r2 / dy + y0) * (r2 / dy + y0)) + FCP.X;
        FRealPoints[1].X = x0 + sqrt(r2 - (r2 / dy + y0) * (r2 / dy + y0)) + FCP.X;
    }

    else if (fabs(y0 - yp) < FTOL)
    {
        double dx = xp - x0;
        FRealPoints[0].X = FRealPoints[1].X = -r2 / dx + FCP.X;
        FRealPoints[0].Y = y0 + sqrt(r2 - (r2 / dx + x0) * (r2 / dx + x0)) + FCP.Y;
        FRealPoints[1].Y = y0 + sqrt(r2 - (r2 / dx + x0) * (r2 / dx + x0)) + FCP.Y;
    }
    else
    {
        g = (xp - x0) / (yp - y0);
        g2 = g * g;
        A = r2 / (xp - x0) - y0 / g;
        A2 = A * A;

        a = (1 + g2);
        b = -2 * (x0 + A * g2);
        c = x0 * x0 + A2 * g2 - r2;
        if (b * b - 4 * a * c <= FTOL)
        {
            FRealPoints.Clear();
            return FRealPoints;
        }
        FRealPoints[0].X = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
        FRealPoints[0].Y = -FRealPoints[0].X * g + r2 / (yp - y0);
        FRealPoints[1].X = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
        FRealPoints[1].Y = -FRealPoints[1].X * g + r2 / (yp - y0);
    }
    FRealPoints[0] += FCP;
    FRealPoints[1] += FCP;

    return FRealPoints;
}

bool MkCircle::operator==(MkCircle &c)
{
    return FCP == c.FCP && fabs(FRadius - c.FRadius) < EPS &&
           fabs(FCircleArea - c.FCircleArea) < EPS && FRealPoints == c.FRealPoints;
}

bool MkCircle::operator!=(MkCircle &c)
{
    return !operator==(c);
}

#ifdef __BCPLUSPLUS__
void MkCircle::Draw(TObject *Sender)
{
    TColor C;
    if (String(Sender->ClassName()) == String("MkPaintBox"))
    {
        MkPaintBox *pb = (MkPaintBox *)Sender;
        MkPoint cp = GetCenter();
        C = pb->Canvas->Pen->Color;
        pb->Canvas->Pen->Color = Color;
        pb->Circle2D(cp.X, cp.Y, GetRadius());
        pb->Canvas->Pen->Color = C;
    }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkCircle::Draw(MkPaint *pb)
{
    MkPoint cp;
    cp = GetCenter();
    pb->Circle2D(cp.X, cp.Y, GetRadius());
}
#endif

//---------------------------------------------------------------------------
// MkCircles::MkCircles(int size)
// {
//     //  try {
//     if (size <= 0)
//     {
//         MkDebug("::MkCircles - MkCircles(int size)");
//         return;
//     }

//     FSize = size;
//     FCircle = new MkCircle[FSize];
//     //  }
//     //  catch () {
//     //    ShowMessage(AnsiString(E.ClassName())+ E.Message);
//     //  }
// }

// void MkCircles::Initialize(int size)
// {
//     Clear();
//     FSize = size;
//     if (FSize == 0)
//     {
//         FCircle = NULL;
//         return;
//     }
//     FCircle = new MkCircle[FSize];
// }

// void MkCircles::Clear()
// {
//     FSize = 0;
//     if (FCircle)
//         delete[] FCircle;
//     FCircle = NULL;
// }

// MkCircle &MkCircles::operator[](int i)
// {
//     if (i >= 0 && i < FSize)
//         return FCircle[i];
//     else
//         return NullCircle;
// }

// MkCircles &MkCircles::operator=(MkCircles &circles)
// {
//     int i;
//     FSize = circles.FSize;
//     this->FCircle = new MkCircle[FSize];

//     for (i = 0; i < FSize; i++)
//         this->FCircle[i] = circles.FCircle[i];

//     return *this;
// }

// #ifdef __BCPLUSPLUS__
// void MkCircles::Draw(TObject *Sender)
// {
//     for (int i = 0; i < FSize; i++)
//         FCircle[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkCircles::Draw(MkPaint *pb)
// {
//     for (int i = 0; i < FSize; i++)
//         FCircle[i].Draw(pb);
// }
// #endif
