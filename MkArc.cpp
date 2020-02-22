//---------------------------------------------------------------------------
// This module is general purposed simple graphic class to store, draw,
// manipulate object. It is well suited to VCL component, but not restricted.
// It forms the base for the higher level class, such as tunnel component.
//
// Copyright (c) 1999 Myung Kyu Song, ESCO Consultant Co., Ltd.
#include "MkArc.hpp"

MkArc NullArc(0);

//---------------------------------------------------------------------------
MkArc::MkArc(int)
{
  FCP.X = 0;
  FCP.Y = 0;
  FRadius = 0;
  className = "MkArc";
}

MkArc::MkArc(double cx, double cy, double radius, double start_ang, double end_ang)
    : MkCircle(cx, cy, radius)
{
  FStartAng = start_ang;
  FEndAng = end_ang;
  CalStartPoint();
  CalEndPoint();
  className = "MkArc";
}

MkArc::MkArc(MkPoint cp, double radius, double start_ang, double end_ang)
    : MkCircle(cp, radius)
{
  FStartAng = start_ang;
  FEndAng = end_ang;
  CalStartPoint();
  CalEndPoint();
  className = "MkArc";
}

#ifdef __BCPLUSPLUS__
MkArc::MkArc(double cx, double cy, double radius, double start_ang, double end_ang, TColor C)
    : MkCircle(cx, cy, radius)
{
  FStartAng = start_ang;
  FEndAng = end_ang;
  Color = C;
  CalStartPoint();
  CalEndPoint();
}

MkArc::MkArc(MkPoint cp, double radius, double start_ang, double end_ang, TColor C)
    : MkCircle(cp, radius)
{
  FStartAng = start_ang;
  FEndAng = end_ang;
  Color = C;
  CalStartPoint();
  CalEndPoint();
}
#endif

MkArc::MkArc() : MkCircle()
{
  FStartAng = 0;
  FEndAng = 0;
  className = "MkArc";
}

MkArc::MkArc(MkPoint p1, MkPoint p2, MkLine line)
{
  MkLine rl1 = MkLine(p1, p2);
  MkLine rl2 = !rl1;
  MkPoint rp1 = line & rl2;
  int cnt = 0;
  while (rp1 == NullPoint && cnt < 20)
  {
    line.Extend(2);
    rp1 = line & rl2;
    cnt++;
  }
  MkLine rl3 = MkLine(rp1, p1);
  MkLine rl4 = MkLine(rp1, p2);

  double len1 = rl3.GetLength();
  double len2 = rl4.GetLength();
  double theta1 = rl3.GetTheta();
  double theta2 = rl4.GetTheta();

  if (fabs(len1 - len2) < EPS)
    *this = NullArc;
  if (cnt == 20 && rp1 == NullPoint)
    *this = NullArc;
  else
  {
    SetCenter(rp1);
    SetRadius(len1);
    SetStartAng(theta1);
    SetEndAng(theta2);
    CalStartPoint();
    CalEndPoint();
  }
  className = "MkArc";
}

void MkArc::SetCenter(double cx, double cy)
{
  MkCircle::SetCenter(cx, cy);
  CalStartPoint();
  CalEndPoint();
}

void MkArc::SetCenter(MkPoint cp)
{
  MkCircle::SetCenter(cp);
  CalStartPoint();
  CalEndPoint();
}

void MkArc::SetRadius(double radius)
{
  MkCircle::SetRadius(radius);
  CalStartPoint();
  CalEndPoint();
}

void MkArc::SetStartAng(double start_ang)
{
  FStartAng = start_ang;
  CalStartPoint();
  CalEndPoint();
}

void MkArc::SetEndAng(double end_ang)
{
  FEndAng = end_ang;
  CalStartPoint();
  CalEndPoint();
}

double MkArc::GetStartAng()
{
  return FStartAng;
}

double MkArc::GetEndAng()
{
  return FEndAng;
}

double MkArc::GetArea()
{
  if (FStartAng > FEndAng)
    FEndAng += 360;
  return MkCircle::GetArea() * fabs(FStartAng - FEndAng) / 360.0;
}

double MkArc::GetTriArea()
{
  MkTriangle rt(GetCenter(), FStartPoint, FEndPoint);
  return rt.GetArea();
}

bool MkArc::isWithInArc(MkPoint rp)
{
  MkLine rl(FCP, rp);
  return ((*this)(0) * rl > 0 && (*this)(1) * rl < 0) && (GetRadius() > rl.GetLength());
}

bool MkArc::isWithInAng(MkPoint rp)
{
  MkLine rl(FCP, rp);
  return ((*this)(0) * rl > 0 && (*this)(1) * rl < 0);
}

double MkArc::CrossProduct()
{
  MkLine rl1((*this)[0], (*this)[1]);
  MkLine rl2((*this)[0], (*this)[2]);
  return rl1 * rl2;
}
double MkArc::GetCrownArea()
{
  return GetArea() - GetTriArea();
}

MkLine MkArc::operator()(int i)
{
  if (i == 0)
    return MkLine((*this)[0], (*this)[1]);
  else if (i == 1)
    return MkLine((*this)[0], (*this)[2]);
  else
  {
    MkDebug("MkArc::operator() Bad Index");
    throw std::exception();
    //return NullLine;
  }
}

MkPoint &MkArc::operator[](int i)
{
  if (i == 0)
    return FCP;
  else if (i == 1)
    return FStartPoint;
  else if (i == 2)
    return FEndPoint;
  else
  {
    MkDebug("MkArc::operator[] Bad Index");
    return NullPoint;
  }
}

void MkArc::CalStartPoint()
{
  MkPoint cp;
  double x, y;

  cp = GetCenter();
  x = GetRadius() * cos(FStartAng * M_PI / 180.0);
  y = GetRadius() * sin(FStartAng * M_PI / 180.0);
  FStartPoint.X = cp.X + x;
  FStartPoint.Y = cp.Y + y;
  return;
}

void MkArc::CalEndPoint()
{
  MkPoint cp;
  double x, y;

  cp = GetCenter();
  x = FRadius * cos(FEndAng * M_PI / 180.0);
  y = FRadius * sin(FEndAng * M_PI / 180.0);
  FEndPoint.X = cp.X + x;
  FEndPoint.Y = cp.Y + y;
  return;
}

void MkArc::CalStartAngle()
{
  MkLine rl((*this)[0], (*this)[1]);
  FStartAng = rl.GetTheta();
  return;
}

void MkArc::CalEndAngle()
{
  MkLine rl((*this)[0], (*this)[2]);
  FEndAng = rl.GetTheta();
  return;
}

void MkArc::ReCalcAng()
{
  CalStartAngle();
  CalEndAngle();
}

void MkArc::ReCalcPoint()
{
  CalStartPoint();
  CalEndPoint();
}

MkArc &MkArc::operator=(MkArc &rc)
{
  this->MkShape::operator=((MkShape &)rc);
  FCP.X = rc.FCP.X;
  FCP.Y = rc.FCP.Y;
  FRadius = rc.FRadius;
#ifdef __BCPLUSPLUS__
  Color = rc.Color;
#endif
  CalArea();
  FStartAng = rc.FStartAng;
  FEndAng = rc.FEndAng;
  CalStartPoint();
  CalEndPoint();
  return (*this);
}

#ifdef __BCPLUSPLUS__
void MkArc::Draw(TObject *Sender)
{
  TColor C;
  if (String(Sender->ClassName()) == String("MkPaintBox"))
  {
    MkPaintBox *pb = (MkPaintBox *)Sender;
    C = pb->Canvas->Pen->Color;
    pb->Canvas->Pen->Color = Color;

    pb->Arc2D((*this)[0].X, (*this)[0].Y, GetRadius(), GetStartAng(), GetEndAng());
    pb->MoveTo3D((*this)[0].X, (*this)[0].Y, 0);
    pb->LineTo3D((*this)[1].X, (*this)[1].Y, 0);
    pb->MoveTo3D((*this)[0].X, (*this)[0].Y, 0);
    pb->LineTo3D((*this)[2].X, (*this)[2].Y, 0);
    pb->Canvas->Pen->Color = C;
  }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkArc::Draw(MkPaint *pb)
{
  pb->Arc2D((*this)[0].X, (*this)[0].Y, GetRadius(), GetStartAng(), GetEndAng());
  pb->MoveTo3D((*this)[0].X, (*this)[0].Y, 0);
  pb->LineTo3D((*this)[1].X, (*this)[1].Y, 0);
  pb->MoveTo3D((*this)[0].X, (*this)[0].Y, 0);
  pb->LineTo3D((*this)[2].X, (*this)[2].Y, 0);
}
#endif

//---------------------------------------------------------------------------
