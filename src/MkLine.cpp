//---------------------------------------------------------------------------
// This module is general purposed simple graphic class to store, draw,
// manipulate object. It is well suited to VCL component, but not restricted.
// It forms the base for the higher level class, such as tunnel component.
//
// Copyright (c) 1999 Myung Kyu Song, ESCO Consultant Co., Ltd.
#include "MkLine.hpp"

extern MkPoint NullPoint;
MkLine NullLine(NullPoint, NullPoint);
MkLines NullLines;
//---------------------------------------------------------------------------
MkLine::MkLine(MkPoint sp, MkPoint ep)
{
   StartPoint = sp;
   EndPoint = ep;
   isSelected = false;
   isFinite = true;
   Tol = EPS;
   CalCoeff();
   className = "MkLine";
}

MkLine::MkLine(double sx, double sy, double ex, double ey)
{
   StartPoint.X = sx;
   StartPoint.Y = sy;
   EndPoint.X = ex;
   EndPoint.Y = ey;
   isSelected = false;
   isFinite = true;
   Tol = EPS;
   CalCoeff();
   className = "MkLine";
}

void MkLine::SetLine(MkPoint sp, MkPoint ep)
{
   StartPoint = sp;
   EndPoint = ep;
   isFinite = true;
   Tol = EPS;
   CalCoeff();
}

void MkLine::SetLine(double sx, double sy, double ex, double ey)
{
   StartPoint.X = sx;
   StartPoint.Y = sy;
   EndPoint.X = ex;
   EndPoint.Y = ey;
   isFinite = true;
   Tol = EPS;
   CalCoeff();
}

MkLine::MkLine()
{
   A = 0;
   B = 0;
   C = 0, Theta = 0;
   Length = 0;
   StartPoint.X = 0;
   StartPoint.Y = 0;
   EndPoint.X = 0;
   EndPoint.Y = 0;
   EndPoint.Z = 0;
   EndPoint.Z = 0;
   Tol = EPS;
   isFinite = false;
   className = "MkLine";
}

void MkLine::CalLength()
{
   if (StartPoint == EndPoint)
      Length = 0;
   else
      Length = sqrt(pow(StartPoint.X - EndPoint.X, 2) + pow(StartPoint.Y - EndPoint.Y, 2) + pow(StartPoint.Z - EndPoint.Z, 2));
}

void MkLine::CalTheta()
{
   CalLength();
   if (Length < 1.0e-6)
   {
      Theta = 0;
      //      ShowMessage("MkLine::Length is below the limit, Can't calculate theta.");
      MkDebug("MkLine::Length is below the limit, Can't calculate theta.");
      return;
   }
   if (fabs(DeltaX()) < 1.0e-3)
   {
      //if (fabs(NormDeltaX()) < 1.0e-5 && fabs(NormDeltaY()) < 1.0e-5) {
      //         ShowMessage("MkLine::CalTheta Two points are the same, Can't calculate theta.");
      //MkDebug("MkLine::CalTheta Two points are the same, Can't calculate theta.");
      //Theta = 0;
      //return;
      //}
      //      if (NormDeltaX() < 1.0e-6)
      Theta = DeltaY() > 0 ? 90 : 270;
   }
   else
      Theta = atan2(DeltaY(), DeltaX()) * 180 / M_PI;
   AdjustTheta();
}

void MkLine::CalCoeff() // Ax + By = C for 2d, (x-x0)/L = (y-y0)/M = (z-z0)/N = t for 3d
{
   double len;
   L = M = N = 0;

   L = EndPoint.X - StartPoint.X;
   M = EndPoint.Y - StartPoint.Y;
   N = EndPoint.Z - StartPoint.Z;

   len = sqrt(L * L + M * M + N * N);
   if (len > 0.001)
   {
      L /= len;
      M /= len;
      N /= len;
   }

   if (fabs(DeltaX()) < 1.0e-7 && fabs(DeltaY()) < 1.0e-7)
   {
      A = B = C = 0;
      // Can't define the line, it's a point!
   }
   else if (fabs(DeltaX()) < 1.0e-7)
   {
      A = 1;
      B = 0;
      C = (StartPoint.X + EndPoint.X) / 2.0;
   }
   else if (fabs(DeltaY()) < 1.0e-7)
   {
      A = 0;
      B = 1;
      C = (StartPoint.Y + EndPoint.Y) / 2.0;
   }
   else if (!(StartPoint.X || StartPoint.Y))
   {
      A = EndPoint.Y;
      B = -EndPoint.X;
      C = 0;
   }
   else if (!(EndPoint.X || EndPoint.Y))
   {
      A = StartPoint.Y;
      B = -StartPoint.X;
      C = 0;
   }
   else if (fabs(StartPoint.X * EndPoint.Y - EndPoint.X * StartPoint.Y) < 1.0e-7)
   {
      A = -DeltaY();
      B = DeltaX();
      C = 0;
   }
   else
   {
      A = (EndPoint.Y - StartPoint.Y) / (StartPoint.X * EndPoint.Y - EndPoint.X * StartPoint.Y);
      B = (EndPoint.X - StartPoint.X) / (EndPoint.X * StartPoint.Y - EndPoint.Y * StartPoint.X);
      C = 1;
   }
}

double MkLine::CalDist3D(MkPoint rp)
{
   MkPoint mp, sp, ep;
   double dist;

   sp = StartPoint;
   ep = EndPoint;
   mp = ep - sp;
   mp.Normalize();

   double atheta;
   atheta = asin(mp.X);
   atheta = fabs(sin(atheta) - mp.X) < 0.001 ? atheta * 180 / M_PI : (atheta + M_PI) * 180 / M_PI;

   mp.X = rp.X;
   mp.Y = rp.Y;
   mp.Z = rp.Z;

   mp.Translate(-sp.X, -sp.Y, -sp.Z);
   mp.Rotate((90 - atheta), 0, 0);

   return dist = sqrt(mp.Y * mp.Y + mp.Z * mp.Z);
}

double MkLine::CalDist3D(double x, double y, double z)
{
   MkPoint mp, sp, ep, rp(x, y, z);
   double dist;

   sp = StartPoint;
   ep = EndPoint;
   mp = ep - sp;
   mp.Normalize();

   double atheta;
   atheta = asin(mp.X);
   atheta = fabs(sin(atheta) - mp.X) < 0.001 ? atheta * 180 / M_PI : (atheta + M_PI) * 180 / M_PI;

   mp.X = rp.X;
   mp.Y = rp.Y;
   mp.Z = rp.Z;

   mp.Translate(-sp.X, -sp.Y, -sp.Z);
   mp.Rotate((90 - atheta), 0, 0);

   return dist = sqrt(mp.Y * mp.Y + mp.Z * mp.Z);
}

MkPoint &MkLine::GetNearestPnt(MkPoint rp) // Test is qued...^^
{
   static MkPoint mp;
   MkPoint sp, ep;

   sp = StartPoint;
   ep = EndPoint;
   mp = ep - sp;
   mp.Normalize();

   double atheta;
   atheta = asin(mp.X);
   atheta = fabs(sin(atheta) - mp.X) < 0.001 ? atheta * 180 / M_PI : (atheta + M_PI) * 180 / M_PI;

   mp.X = rp.X;
   mp.Y = rp.Y;
   mp.Z = rp.Z;

   mp.Translate(-sp.X, -sp.Y, -sp.Z);
   mp.Rotate((90 - atheta), 0, 0);

   ep.Translate(-sp.X, -sp.Y, -sp.Z);
   ep.Rotate((90 - atheta), 0, 0);

   if (isFinite)
   {
      if (mp.X < -0.001)
         return StartPoint;
      else if (mp.X > ep.X + 0.001)
         return EndPoint;
   }

   mp.Y = 0;
   mp.Z = 0;

   mp.Rotate((atheta - 90), 0, 0);
   mp.Translate(sp.X, sp.Y, sp.Z);

   return mp;
}

MkPoint &MkLine::GetNearestPnt(double x, double y, double z)
{
   MkPoint rp(x, y, z);
   static MkPoint mp;
   mp = GetNearestPnt(rp);
   return mp;
}

bool MkLine::GetIntParam(MkLine &rl, double &t) // full 3d & check it out.  test is qued...^^
{
   double l1, m1, n1;
   double l2, m2, n2;
   double x1, y1, z1;
   double x2, y2, z2;
   double t1max, t2max;
   double len;
   MkPoint sp1, ep1, sp2, ep2;
   MkMatrix<double> mat(2, 2);
   MkVector<double> vec(2);

   t = 0;
   if (!IsInSamePlane(rl))
      return false;

   sp1 = StartPoint;
   ep1 = EndPoint;
   sp2 = rl.StartPoint;
   ep2 = rl.EndPoint;

   x1 = sp1.X;
   y1 = sp1.Y;
   z1 = sp1.Z;

   x2 = sp2.X;
   y2 = sp2.Y;
   z2 = sp2.Z;

   l1 = ep1.X - sp1.X;
   m1 = ep1.Y - sp1.Y;
   n1 = ep1.Z - sp1.Z;

   t1max = len = sqrt(l1 * l1 + m1 * m1 + n1 * n1);
   l1 /= len;
   m1 /= len;
   n1 /= len;

   l2 = ep2.X - sp2.X;
   m2 = ep2.Y - sp2.Y;
   n2 = ep2.Z - sp2.Z;

   t2max = len = sqrt(l2 * l2 + m2 * m2 + n2 * n2);
   l2 /= len;
   m2 /= len;
   n2 /= len;

   if (fabs(l1 - l2) < EPS && fabs(m1 - m2) < EPS && fabs(n1 - n2) < EPS)
      return false;
   else if (fabs(l1) < EPS && fabs(l2) < EPS)
   {
      mat(0, 0) = m1;
      mat(0, 1) = n1;
      mat(1, 0) = -m2;
      mat(1, 1) = -n2;
      vec(0) = y2 - y1;
      vec(1) = z2 - z1;
      mat.Solve(vec, stLUD);
   }
   else if (fabs(m1) < EPS && fabs(m2) < EPS)
   {
      mat(0, 0) = l1;
      mat(0, 1) = n1;
      mat(1, 0) = -l2;
      mat(1, 1) = -n2;
      vec(0) = x2 - x1;
      vec(1) = z2 - z1;
      mat.Solve(vec, stLUD);
   }
   else
   {
      mat(0, 0) = l1;
      mat(0, 1) = -l2;
      mat(1, 0) = m1;
      mat(1, 1) = -m2;
      vec(0) = x2 - x1;
      vec(1) = y2 - y1;
      mat.Solve(vec, stLUD);
   }

   t = vec(0) / GetLength();
   return true;
}

bool MkLine::GetIntParam(MkPoint p, double &t) // full 3d & check it out.  test is qued...^^
{
   MkPoint sp, ep;
   sp = StartPoint;
   ep = EndPoint;

   t = MAX_VAL;
   if (GetLength() < EPS || CalDist(p) > EPS)
      return false;

   if (fabs(ep.X - sp.X) < EPS && fabs(p.X - sp.X) > EPS)
      return false;
   if (fabs(ep.Y - sp.Y) < EPS && fabs(p.Y - sp.Y) > EPS)
      return false;
   if (fabs(ep.Z - sp.Z) < EPS && fabs(p.Z - sp.Z) > EPS)
      return false;

   if (fabs(ep.X - sp.X) > EPS)
   {
      t = (p.X - sp.X) / (ep.X - sp.X);
      return true;
   }
   if (fabs(ep.Y - sp.Y) > EPS)
   {
      t = (p.Y - sp.Y) / (ep.Y - sp.Y);
      return true;
   }
   if (fabs(ep.Z - sp.Z) > EPS)
   {
      t = (p.Z - sp.Z) / (ep.Z - sp.Z);
      return true;
   }

   return false;
}

MkPoint &MkLine::GetIntPoint(MkLine &rl)
{
   double t;
   bool flag;
   if (!(flag = GetIntParam(rl, t)))
      return NullPoint;
   if (t < 0.0 || t > 1.0)
      return NullPoint;
   return GetDivision(t);
}

double MkLine::GetLength()
{
   CalLength();
   return Length;
}

double MkLine::GetTheta()
{
   CalTheta();
   return Theta;
}

MkPoint &MkLine::GetDivision(double d)
{
   double x, y, z;
   static MkPoint mp;
   if (GetLength() == 0)
      return (*this)[0];

   x = (*this)[0].X + d * ((*this)[1].X - (*this)[0].X);
   y = (*this)[0].Y + d * ((*this)[1].Y - (*this)[0].Y);
   z = (*this)[0].Z + d * ((*this)[1].Z - (*this)[0].Z);
   mp.SetPoint(x, y, z);
   return mp;
}

void MkLine::Extend(double f)
{
   MkPoint P1, P2, MP;
   MP = GetMiddlePoint();
   P1 = MkLine(MP, StartPoint).GetDivision(f);
   P2 = MkLine(MP, EndPoint).GetDivision(f);
   StartPoint = P1;
   EndPoint = P2;
}
void MkLine::Clear()
{
   StartPoint = NullPoint;
   EndPoint = NullPoint;
   Theta = Length = 0;
   A = B = C = 0; // A*X + B*Y = C  (C = 0 or 1)
   L = M = N = 0; // direction cosine
   isSelected = false;
   isFinite = false;

   Tol = EPS;
}

MkPoint &MkLine::operator[](int i)
{
   if (i == 0)
      return StartPoint;
   else if (i == 1)
      return EndPoint;
   else
   {
#ifdef __BCPLUSPLUS__
      ShowMessage("MkPoint::operator, error index");
#else
      MkDebug("MkPoint::operator, error index");
#endif
      return NullPoint;
   }
}

bool MkLine::operator&&(MkLine rline1)
{
   return IsIntersect(rline1);
}

bool MkLine::operator==(MkLine rline1)
{
   bool flag = true;
   for (int i = 0; i < 2 && flag; i++)
      flag = flag && ((*this)[i] == rline1[i]);
   return flag;
}

bool MkLine::operator!=(MkLine rline1)
{
   return !(*this == rline1);
}

MkPoint &MkLine::operator&(MkLine &line) // have to revisit...***
{
   //   if (!IsIntersect(line)) return NullPoint;
   static MkPoint pnt;
   pnt.X = (line.B * C - B * line.C) /
           (line.B * A - B * line.A);
   pnt.Y = (line.A * C - A * line.C) /
           (line.A * B - A * line.B);
   return pnt;

   if (fabs(L) < EPS && fabs(M) < EPS && fabs(line.L) < EPS && fabs(line.M) < EPS)
   {
      if (fabs(StartPoint.X - line.StartPoint.X) < EPS &&
          fabs(StartPoint.Y - line.StartPoint.Y) < EPS)
         return StartPoint;
      else
         NullPoint;
   }
   else
   {
      MkMatrix<double> mat(2, 2);
      MkVector<double> t(2);
      mat(0, 0) = L;
      mat(0, 1) = -line.L;
      mat(1, 0) = M;
      mat(1, 1) = -line.M;
      t[0] = line.StartPoint.X - StartPoint.Y;
      t[1] = line.StartPoint.Y - StartPoint.Y;
      mat.Solve(t);
      if (fabs(StartPoint.Z - line.StartPoint.Z + N * t[0] - line.N * t[1]) < EPS)
      {
         pnt.SetPoint(StartPoint.X + L * t[0], StartPoint.Y + M * t[0], StartPoint.Z + N * t[0]);
         return pnt;
      }
      else
         NullPoint;
   }
}

double MkLine::operator+(MkLine rl) //dot product
{
   double a, b, c, d, l;
   a = (*this)[1].X - (*this)[0].X;
   b = (*this)[1].Y - (*this)[0].Y;
   if (a * a + b * b < EPS)
      return 0;
   l = sqrt(a * a + b * b);
   a = a / l;
   b = b / l;
   c = rl[1].X - rl[0].X;
   d = rl[1].Y - rl[0].Y;
   if (c * c + d * d < EPS)
      return 0;
   l = sqrt(c * c + d * d);
   c = c / l;
   d = d / l;

   return a * d + b * c;
}

double MkLine::operator*(MkLine rl) //cross product
{
   double a, b, c, d, l;
   a = (*this)[1].X - (*this)[0].X;
   b = (*this)[1].Y - (*this)[0].Y;
   if (a * a + b * b < EPS)
      return 0;
   l = sqrt(a * a + b * b);
   a = a / l;
   b = b / l;
   c = rl[1].X - rl[0].X;
   d = rl[1].Y - rl[0].Y;
   if (c * c + d * d < EPS)
      return 0;
   l = sqrt(c * c + d * d);
   c = c / l;
   d = d / l;
   return a * d - b * c;
}

double MkLine::operator*(MkPoint rp) //cross product
{
   double a, b, c, d, l;

   MkLine rl((*this)[0], rp);

   a = (*this)[1].X - (*this)[0].X;
   b = (*this)[1].Y - (*this)[0].Y;
   if (a * a + b * b < EPS)
      return 0;
   l = sqrt(a * a + b * b);
   a = a / l;
   b = b / l;
   c = rl[1].X - rl[0].X;
   d = rl[1].Y - rl[0].Y;
   if (c * c + d * d < EPS)
      return 0;
   l = sqrt(c * c + d * d);
   c = c / l;
   d = d / l;
   return a * d - b * c;
}

MkLine &MkLine::operator*(MkMatrix4<double> &rm)
{
   StartPoint = StartPoint * rm;
   EndPoint = EndPoint * rm;
   CalCoeff();
   CalTheta();
   return *this;
}

double MkLine::operator+=(MkPoint rp)
{
   StartPoint = StartPoint + rp;
   EndPoint = EndPoint + rp;
   CalCoeff();
   return 1;
}

double MkLine::operator-=(MkPoint rp)
{
   StartPoint = StartPoint - rp;
   EndPoint = EndPoint - rp;
   CalCoeff();
   return 1;
}

bool MkLine::IsIntersect(MkLine &rl)
{
   if (!IsInSamePlane(rl))
      return false;
   //   if (!isFinite && !rl.isFinite) return true;

   MkVector<double> v12, v13, v14, v31, v32, v34, v142, v123, v314, v342;
   MkPoint p1, p2, p3, p4;
   MkLine l13, l14, l31, l32;

   p1 = StartPoint;
   p2 = EndPoint;
   p3 = rl.StartPoint;
   p4 = rl.EndPoint;

   if (p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4)
      return true;

   l13.SetLine(p1, p3);
   l14.SetLine(p1, p4);
   l31.SetLine(p3, p1);
   l32.SetLine(p3, p2);

   CalCoeff();
   rl.CalCoeff();

   v12 = GetVector();
   v34 = rl.GetVector();

   if (v12 == v34)
      return false;

   v13 = l13.GetVector();
   v14 = l14.GetVector();
   v31 = l31.GetVector();
   v32 = l32.GetVector();

   v14.Cross(v12, v142);
   v142.Unify();
   v12.Cross(v13, v123);
   v123.Unify();
   v31.Cross(v34, v314);
   v314.Unify();
   v34.Cross(v32, v342);
   v342.Unify();

   if (!isFinite && !rl.isFinite)
   {
      /*
      MkDouble t;
      int i,j,cnt=0;
      bool flag=true;

      if (fabs(L)>EPS*1000) cnt++;
      if (fabs(M)>EPS*1000) cnt++;
      if (fabs(N)>EPS*1000) cnt++;

      t.Initialize(cnt);

      cnt=0;
      if (fabs(L)>EPS*1000) {t(cnt) = (p3.X-p1.X)/L;cnt++;} else if (fabs(p3.X-p1.X)>EPS*1000){flag=false;}
      if (fabs(M)>EPS*1000) {t(cnt) = (p3.Y-p1.Y)/M;cnt++;} else if (fabs(p3.Y-p1.Y)>EPS*1000){flag=false;}
      if (fabs(N)>EPS*1000) {t(cnt) = (p3.Z-p1.Z)/N;cnt++;} else if (fabs(p3.Z-p1.Z)>EPS*1000){flag=false;}

      for (i=0;i<t.getSzX()-1;i++) {
        for(j=i+1;j<t.getSzX();j++) {
          flag = flag && fabs(t(i)-t(j))<EPS*1000;
        }
      }
      return fabs(v12*v34) - 1 < EPS && flag;
*/
      return true;
   }
   else if (!isFinite && rl.isFinite)
   {
      return v142 * v123 > 0;
   }

   else if (!rl.isFinite && isFinite)
   {
      return v314 * v342 > 0;
   }

   else
   {
      return v142 * v123 > 0 && v314 * v342 > 0;
   }
}

bool MkLine::IsIn(MkPoint rp)
{
   bool flag;
   double sx, sy, ex, ey;
   sx = min((*this)[0].X, (*this)[1].X) - Tol;
   sy = min((*this)[0].Y, (*this)[1].Y) - Tol;
   ex = max((*this)[0].X, (*this)[1].X) + Tol;
   ey = max((*this)[0].Y, (*this)[1].Y) + Tol;

   if (is_eq((*this)[0].X, (*this)[1].X))
      flag = (rp.Y >= sy && rp.Y <= ey) ? true : false;
   else if (is_eq((*this)[0].Y, (*this)[1].Y))
      flag = (rp.X >= sx && rp.X <= ex) ? true : false;
   else if (rp.X >= sx && rp.X <= ex &&
            rp.Y >= sy && rp.Y <= ey)
      flag = true;
   else
      flag = false;
   return flag;
}

bool MkLine::IsIn(double x, double y)
{
   bool flag;
   double sx, sy, ex, ey;
   char str[256];

   sprintf(str, "Tol of line is %f\n", Tol);
   MkDebug(str);

   sx = min((*this)[0].X, (*this)[1].X) - Tol;
   sy = min((*this)[0].Y, (*this)[1].Y) - Tol;
   ex = max((*this)[0].X, (*this)[1].X) + Tol;
   ey = max((*this)[0].Y, (*this)[1].Y) + Tol;

   if (is_eq((*this)[0].X, (*this)[1].X))
      flag = (y >= sy && y <= ey) ? true : false;
   else if (is_eq((*this)[0].Y, (*this)[1].Y))
      flag = (x >= sx && x <= ex) ? true : false;
   else if (x >= sx && x <= ex &&
            y >= sy && y <= ey)
      flag = true;
   else
      flag = false;
   return flag;
}

bool MkLine::IsInLine(MkPoint rp)
{
   double dist;
   CalCoeff();
   CalLength();
   dist = CalDist(rp);
   if (IsIn(rp) && (dist < Tol))
      return true;
   else
      return false;
}

bool MkLine::IsInLine(double x, double y)
{
   double dist;
   CalCoeff();
   CalLength();
   dist = CalDist(x, y);
   if (IsIn(x, y) && (dist < Tol))
      return true;
   else
      return false;
}

bool MkLine::IsInSamePlane(MkLine rl)
{
   MkPoint sp1, ep1, sp2, ep2;
   MkVector<double> v1, v2, norm;
   MkVector<double> v12(3), v21(3), norm12, norm21;
   MkLine l12, l21;

   v1 = GetVector();
   v2 = rl.GetVector();

   if (v1 == v2)
      return true;

   v1.Cross(v2, norm);

   sp1 = StartPoint;
   ep1 = EndPoint;
   sp2 = rl.StartPoint;
   ep2 = rl.EndPoint;

   if (sp1 == sp2)
      l12.SetLine(sp1, ep2);
   else
      l12.SetLine(sp1, sp2);

   if (sp2 == sp1)
      l21.SetLine(sp2, ep1);
   l21.SetLine(sp2, sp1);

   v12 = l12.GetVector();
   v21 = l21.GetVector();

   if ((v1 == v12) || (v1 == v21) || (v2 == v12) || (v2 == v21))
      return true;

   v1.Cross(v12, norm12);
   v21.Cross(v2, norm21);

   norm12.Normalize();
   norm21.Normalize();

   return norm12 == norm21;
   /*
   double l1,m1,n1;
   double l2,m2,n2;
   double l3,m3,n3;
   double l4,m4,n4;
   double len;
   MkPoint sp1,ep1,sp2,ep2;

   sp1 = StartPoint;
   ep1 = EndPoint;
   sp2 = rl.StartPoint;
   ep2 = rl.EndPoint;

   l1 = ep1.X - sp1.X;
   m1 = ep1.Y - sp1.Y;
   n1 = ep1.Z - sp1.Z;

   len = sqrt(l1*l1+m1*m1+n1*n1);

   if(len<0.0001) 
     return IsInLine(sp1);

   l1 /= len;m1 /= len;n1 /= len;

   l2 = sp2.X - sp1.X;
   m2 = sp2.Y - sp1.Y;
   n2 = sp2.Z - sp1.Z;

   len = sqrt(l2*l2+m2*m2+n2*n2);

   if(len<0.0001)
     return rl.IsInLine(sp2);

   l2 /= len;m2 /= len;n2 /= len;

   if (fabs(l1-l2) < 0.01 && fabs(m1-m2) < 0.01 && fabs(n1-n2) < 0.01) return true;

   l3 = m1*n2 - m2*n1;
   m3 = n1*l2 - n2*l1;
   n3 = l1*m2 - l2*m1;

   len = sqrt(l3*l3+m3*m3+n3*n3);

   if(len<0.0001) return true;

   l3 /= len;m3 /= len;n3 /= len;

   if (fabs((l3*sp2.X+m3*sp2.Y+n3*sp2.Z)-(l3*sp1.X+m3*sp1.Y+n3*sp1.Z))<0.001) return true;

   l2 = ep2.X - sp1.X;
   m2 = ep2.Y - sp1.Y;
   n2 = ep2.Z - sp1.Z;

   len = sqrt(l2*l2+m2*m2+n2*n2);

   if(len<0.0001) return true;
   
   l2 /= len;m2 /= len;n2 /= len;

   l4 = m1*n2 - m2*n1;
   m4 = n1*l2 - n2*l1;
   n4 = l1*m2 - l2*m1;

   len = sqrt(l4*l4+m4*m4+n4*n4);
   if(len<0) return true;

   l4 /= len;m4 /= len;n4 /= len;

   if (fabs(l3-l4) < 0.01 && fabs(m3-m4) < 0.01 && fabs(n3-n4) < 0.01) return true;
   else return false;*/
}

MkLine MkLine::Translate(MkPoint rp)
{
   StartPoint.Translate(rp);
   EndPoint.Translate(rp);
   CalCoeff();
   CalTheta();
   return *this;
}

MkLine MkLine::Translate(double x, double y, double z)
{
   StartPoint.Translate(x, y, z);
   EndPoint.Translate(x, y, z);
   CalCoeff();
   CalTheta();
   return *this;
}

MkLine MkLine::Rotate(double alpha, double beta, double gamma)
{
   StartPoint.Rotate(alpha, beta, gamma);
   EndPoint.Rotate(alpha, beta, gamma);
   CalCoeff();
   CalTheta();
   return *this;
}

MkLine MkLine::RotateInX(double ang)
{
   StartPoint.RotateInX(ang);
   EndPoint.RotateInX(ang);
   CalCoeff();
   CalTheta();
   return *this;
}

MkLine MkLine::RotateInY(double ang)
{
   StartPoint.RotateInY(ang);
   EndPoint.RotateInY(ang);
   CalCoeff();
   CalTheta();
   return *this;
}

MkLine MkLine::RotateInZ(double ang)
{
   StartPoint.RotateInZ(ang);
   EndPoint.RotateInZ(ang);
   CalCoeff();
   CalTheta();
   return *this;
}

MkLine MkLine::RotateInA(double ang, double l, double m, double n)
{
   StartPoint.RotateInA(ang, l, m, n);
   EndPoint.RotateInA(ang, l, m, n);
   CalCoeff();
   CalTheta();
   return *this;
}

MkLine MkLine::Scale(double sx, double sy, double sz)
{
   StartPoint.Scale(sx, sy, sz);
   EndPoint.Scale(sx, sy, sz);
   CalCoeff();
   CalTheta();
   return *this;
}

#if !defined(_MSC_VER) && !defined(_WINDOWS_) || defined(__BCPLUSPLUS__)
MkLine &MkLine::operator=(const MkLine &rl)
{
   MkShape::operator=((MkShape &)rl);
   StartPoint = rl.StartPoint;
   EndPoint = rl.EndPoint;
   Theta = rl.Theta;
   Length = rl.Length;

   A = rl.A;
   B = rl.B;
   C = rl.C;
   Tol = rl.Tol;
   return (*this);
}
#endif

MkLine &MkLine::operator=(MkLine &rl)
{
   MkShape::operator=((MkShape &)rl);
   StartPoint = rl.StartPoint;
   EndPoint = rl.EndPoint;
   Theta = rl.Theta;
   Length = rl.Length;

   A = rl.A;
   B = rl.B;
   C = rl.C;
   Tol = rl.Tol;
   return (*this);
}

MkLine MkLine::operator!()
{
   CalLength();
   CalTheta();
   CalCoeff();

   MkLine norm;
   norm[0].X = ((*this)[0].X + (*this)[1].X) / 2.0;
   norm[0].Y = ((*this)[0].Y + (*this)[1].Y) / 2.0;

   if (fabs(A) < 1.0e-5)
   {
      norm.B = 0;
      norm.A = 1;
   }
   else if (fabs(B) < 1.0e-5)
   {
      norm.A = 0;
      norm.B = 1;
   }
   else
   {
      norm.A = -B;
      norm.B = A;
   }
   norm.C = norm.A * norm[0].X + norm.B * norm[0].Y;

   if (fabs(norm.A) < 1.0e-5)
   {
      norm[0].X = -100;
      norm[0].Y = norm.C / norm.B;
      norm[1].X = 100;
      norm[1].Y = norm.C / norm.B;
   }
   else if (fabs(norm.B) < 1.0e-5)
   {
      norm[0].Y = -100;
      norm[0].X = norm.C / norm.A;
      norm[1].Y = 100;
      norm[1].X = norm.C / norm.A;
   }
   else if (fabs(norm.A / norm.B) > 1)
   {
      norm[0].Y = -100;
      norm[0].X = -norm.B / norm.A * norm[0].Y + norm.C / norm.A;
      norm[1].Y = 100;
      norm[1].X = -norm.B / norm.A * norm[1].Y + norm.C / norm.A;
   }
   else if (fabs(norm.A / norm.B) < 1)
   {
      norm[0].X = -100;
      norm[0].Y = -norm.A / norm.B * norm[0].X + norm.C / norm.B;
      norm[1].X = 100;
      norm[1].Y = -norm.A / norm.B * norm[1].X + norm.C / norm.B;
   }

   return norm;
}

void MkLine::Out()
{
   char str[256];
   sprintf(str, "(%10.5f,%10.5f,%10.5f)-(%10.5f,%10.5f,%10.5f)\n", StartPoint.X, StartPoint.Y, StartPoint.Z, EndPoint.X, EndPoint.Y, EndPoint.Z);
   MkDebug(str);
}

void MkLine::Draw()
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glPushMatrix();
   //  glTranslatef(-(*this)[0].X,-(*this)[0].Y,-(*this)[0].Z);
   glColor3f(1, 1, 0);
   glLineWidth(10.0);
   glutSolidSphere(1.0, 16, 16);
   glBegin(GL_POLYGON);
   glVertex3f((*this)[0].X, (*this)[0].Y, (*this)[0].Z);
   glVertex3f((*this)[1].X, (*this)[1].Y, (*this)[1].Z);
   glVertex3f((*this)[0].X, (*this)[0].Y, (*this)[0].Z);
   glEnd();
   glFlush();
   glPopMatrix();
}

#ifdef __BCPLUSPLUS__
void MkLine::Draw(TObject *Sender)
{
   TColor C;
   TPenStyle PS;
   if (String(Sender->ClassName()) == String("MkPaintBox"))
   {
      MkPaintBox *pb = (MkPaintBox *)Sender;
      C = pb->Canvas->Pen->Color;
      PS = pb->Canvas->Pen->Style;
      pb->Canvas->Pen->Color = Color;
      //      pb->Canvas->Pen->Style = PenStyle;

      pb->MoveTo3D((*this)[0].X, (*this)[0].Y, (*this)[0].Z);
      pb->LineTo3D((*this)[1].X, (*this)[1].Y, (*this)[1].Z);
      pb->Canvas->Pen->Color = C;
      pb->Canvas->Pen->Style = PS;
   }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkLine::Draw(MkPaint *pb)
{
   pb->MoveTo3D((*this)[0].X, (*this)[0].Y, (*this)[0].Z);
   pb->LineTo3D((*this)[1].X, (*this)[1].Y, (*this)[1].Z);
}
#endif

//---------------------------------------------------------------------------
// MkLines::MkLines(int size, MkLine *rl)
// {
//    if (size < 0)
//    {
// #ifdef __BCPLUSPLUS__
//       ShowMessage("::MkLines - MkLines(int size)");
// #else
//       MkDebug("::MkLines - MkLines(int size)");
// #endif
//       return;
//    }

//    FSize = size;
//    if (FSize == 0)
//    {
//       FRealLine = NULL;
//       return;
//    }

//    FRealLine = new MkLine[FSize];
//    for (int i = 0; i < FSize; i++)
//       FRealLine[i] = rl[i];
// }

// MkLines::MkLines(int size)
// {
//    if (size < 0)
//    {
// #ifdef __BCPLUSPLUS__
//       ShowMessage("::MkLines - MkLines(int size)");
// #else
//       MkDebug("::MkLines - MkLines(int size)");
// #endif

//       return;
//    }

//    FSize = size;
//    if (FSize == 0)
//    {
//       FRealLine = NULL;
//       return;
//    }

//    FRealLine = new MkLine[FSize];
// }

// MkLines::~MkLines()
// {
//    FSize = 0;
//    if (FRealLine)
//    {
//       delete[](MkLine *) FRealLine;
//       FRealLine = NULL;
//    }
// }

// void MkLines::Initialize(int size)
// {
//    if (size < 0)
//    {
// #ifdef __BCPLUSPLUS__
//       ShowMessage("::MkLines - Initialize(int size)");
// #else
//       MkDebug("::MkLines - Initialize(int size)");
// #endif

//       return;
//    }
//    if (FSize == size)
//       return;

//    FSize = size;
//    if (FSize == 0)
//    {
//       if (FRealLine != NULL)
//          delete[](MkLine *) FRealLine;
//       FRealLine = NULL;
//       return;
//    }

//    if (FRealLine != NULL)
//       delete[](MkLine *) FRealLine;
//    FRealLine = new MkLine[FSize];
// }

// void MkLines::Initialize(int size, MkLine *rl)
// {
//    if (size < 0)
//    {
// #ifdef __BCPLUSPLUS__
//       ShowMessage("::MkLines - Initialize(int size)");
// #else
//       MkDebug("::MkLines - Initialize(int size)");
// #endif

//       return;
//    }

//    FSize = size;
//    if (FSize == 0)
//    {
//       if (FRealLine != NULL)
//          delete[](MkLine *) FRealLine;
//       FRealLine = NULL;
//       return;
//    }

//    if (FRealLine != NULL)
//       delete[](MkLine *) FRealLine;
//    FRealLine = new MkLine[FSize];
//    for (int i = 0; i < FSize; i++)
//       FRealLine[i] = rl[i];
// }

// void MkLines::Grow(int size) // Grow allocate extra memory
// {
//    if (size <= 0)
//    {
//       return;
//    }

//    MkLine *pLine = new MkLine[FSize + size];

//    for (int i = 0; i < FSize; i++)
//       pLine[i] = FRealLine[i];

//    if (FRealLine != NULL)
//    {
//       delete[](MkLine *) FRealLine;
//       FRealLine = NULL;
//    }

//    FRealLine = pLine;
//    FSize = FSize + size;
// }

// void MkLines::DeleteSelected()
// {
//    int i, c, del_line_num;
//    del_line_num = 0;
//    for (i = 0; i < FSize; i++)
//    {
//       if (FRealLine[i].GetSelected())
//          del_line_num++;
//    }

//    MkLine *pLine = new MkLine[FSize - del_line_num];
//    for (i = 0, c = 0; i < FSize; i++)
//    {
//       if (!FRealLine[i].GetSelected())
//          pLine[c++] = FRealLine[i];
//    }
//    if (FRealLine != NULL)
//    {
//       delete[](MkLine *) FRealLine;
//       FRealLine = NULL;
//    }
//    FRealLine = pLine;
//    FSize = FSize - del_line_num;
// }

// bool MkLines::Clear()
// {
//    FSize = 0;
//    if (FRealLine)
//    {
//       delete[](MkLine *) FRealLine;
//       FRealLine = NULL;
//    }
//    return true;
// }

// MkLines &MkLines::RotateInX(double ang)
// {
//    for (int i = 0; i < FSize; i++)
//       FRealLine[i].RotateInX(ang);
//    return *this;
// }

// MkLines &MkLines::RotateInY(double ang)
// {
//    for (int i = 0; i < FSize; i++)
//       FRealLine[i].RotateInY(ang);
//    return *this;
// }

// MkLines &MkLines::RotateInZ(double ang)
// {
//    for (int i = 0; i < FSize; i++)
//       FRealLine[i].RotateInZ(ang);
//    return *this;
// }

// MkLines &MkLines::RotateInA(double ang, double l, double m, double n)
// {
//    for (int i = 0; i < FSize; i++)
//       FRealLine[i].RotateInA(ang, l, m, n);
//    return *this;
// }

// MkLine &MkLines::operator[](int i)
// {
//    if (FSize == 0)
//       return NullLine;
//    else if (i >= 0 && i < FSize)
//       return FRealLine[i];
//    else
//       return NullLine;
// }

// MkLines &MkLines::operator*(MkMatrix4<double> &rm)
// {
//    for (int i = 0; i < FSize; i++)
//       FRealLine[i] = FRealLine[i] * rm;
//    return *this;
// }

// MkLines &MkLines::operator=(MkLines &rls)
// {
//    int i;

//    Clear();
//    FSize = rls.FSize;
//    if (FSize == 0)
//    {
//       FRealLine = NULL;
//       return *this;
//    }
//    FRealLine = new MkLine[FSize];

//    for (i = 0; i < FSize; i++)
//       FRealLine[i] = rls.FRealLine[i];

//    return *this;
// }

// bool MkLines::operator==(MkLines &lines)
// {
//    int i;

//    if (FSize != lines.FSize)
//       return false;
//    for (i = 0; i < FSize; i++)
//       if (FRealLine[i] != lines.FRealLine[i])
//          return false;

//    return true;
// }

// void MkLines::Draw()
// {
//    for (int i = 0; i < GetSize(); i++)
//       FRealLine[i].Draw();
// }

// #ifdef __BCPLUSPLUS__
// void MkLines::Draw(TObject *Sender)
// {
//    for (int i = 0; i < GetSize(); i++)
//       FRealLine[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkLines::Draw(MkPaint *pb)
// {
//    for (int i = 0; i < GetSize(); i++)
//       FRealLine[i].Draw(pb);
// }
// #endif

//---------------------------------------------------------------------------

//template class MkContainer<MkLine>;