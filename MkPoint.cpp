//---------------------------------------------------------------------------
// This module is general purposed simple graphic class to store, draw,
// manipulate object. It is well suited to VCL component, but not restricted.
// It forms the base for the higher level class, such as tunnel component.
//
// Copyright (c) 1999 Myung Kyu Song, ESCO Consultant Co., Ltd.
#include "MkPoint.hpp"

// #ifdef __BCPLUSPLUS__
// #include <Extctrls.hpp>
// #include "MkPaintBox.h"
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// #include <Windows.h>
// #include "MkPaint.h"
// #endif

//-----------------------------------------------------------------------------

#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif

MkPoint NullPoint(0, 0, 0);
MkPoints NullPoints(0);

MkPoint &MkPoint::operator=(MkPoint &rp)
{
   X = rp.X;
   Y = rp.Y;
   Z = rp.Z;
#if defined(__GL_H__)
   Color = rp.Color;
#endif
   return (*this);
}

#if !defined(_MSC_VER) && !defined(_WINDOWS_) || defined(__BCPLUSPLUS__)

MkPoint &MkPoint::operator=(const MkPoint &rp)
{
   X = rp.X;
   Y = rp.Y;
   Z = rp.Z;
   return (*this);
}
#endif

double MkPoint::GetAng()
{
   double theta;
   double len = X * X + Y * Y;
   if (len > 0.001)
   {
      len = sqrt(len);
   }
   X = X / len;
   Y = Y / len;

   if (fabs(X) < 0.0001)
   {
      return theta = Y > 0 ? 90 : 270;
   }
   else
   {
      theta = acos(X);
      if (fabs(Y - sin(theta)) < 0.0001)
         return theta * 180 / M_PI;
      else
         return (2 * M_PI - theta) * 180 / M_PI;
   }
}

void MkPoint::GetAng(double &alpha, double &beta, double &gamma)
{
   double len;

   double l, m, n;
   l = X;
   m = Y;
   n = Z;
   len = sqrt(l * l + m * m + n * n);
   l /= len;
   m /= len;
   n /= len;
   alpha = 0;
   beta = acos(n);

   if (fabs(l) > 0.0001)
   {
      gamma = acos(sin(beta) / l);
      if (l != sin(beta) * cos(gamma))
         gamma = -gamma;
      beta = beta * 180 / M_PI;
      gamma = gamma * 180 / M_PI;
   }
   else if (fabs(m) > 0.0001)
   {
      gamma = asin(sin(beta) / m);
      if (l != sin(beta) * sin(gamma))
         gamma = M_PI - gamma;
      beta = beta * 180 / M_PI;
      gamma = gamma * 180 / M_PI;
   }
}

bool MkPoint::operator==(MkPoint &rp)
{
   char str[256];
   sprintf(str, "(%f,%f,%f),(%f,%f,%f)\n", X, Y, Z, rp.X, rp.Y, rp.Z);
   //MkDebug("MkPoint::operator==() is called\n");
   //MkDebug(str);
   return (fabs(X - rp.X) < EPS) && (fabs(Y - rp.Y) < EPS) && (fabs(Z - rp.Z) < EPS);
}

bool MkPoint::operator!=(MkPoint &rp)
{
   return (fabs(X - rp.X) > EPS) || (fabs(Y - rp.Y) > EPS) || (fabs(Z - rp.Z) > EPS);
}

MkPoint &MkPoint::operator*=(MkMatrix4<double> &rm)
{
   static MkPoint rp;

   rp.X = X * rm(0, 0) + Y * rm(0, 1) + Z * rm(0, 2) + 1 * rm(0, 3);
   rp.Y = X * rm(1, 0) + Y * rm(1, 1) + Z * rm(1, 2) + 1 * rm(1, 3);
   rp.Z = X * rm(2, 0) + Y * rm(2, 1) + Z * rm(2, 2) + 1 * rm(2, 3);

   return *this = rp;
}

MkPoint &MkPoint::operator*=(double f)
{
   X *= f;
   Y *= f;
   Z *= f;
   return *this;
}

void MkPoint::Unify()
{
   double len = X * X + Y * Y + Z * Z;
   try
   {
      len = sqrt(len);
      if (len < 0.000001)
      {
         throw len;
      }

      X = X / len;
      Y = Y / len;
      Z = Z / len;
   }
   catch (double f)
   {
      MkDebug("MkPoint Unify length is too small %f", f);
      throw f;
   }
}

MkPoint &MkPoint::Translate(MkPoint &rp)
{
   X = X + rp.X;
   Y = Y + rp.Y;
   Z = Z + rp.Z;
   return *this;
}

MkPoint &MkPoint::Translate(double x, double y, double z)
{
   X = X + x;
   Y = Y + y;
   Z = Z + z;
   return *this;
}

MkPoint &MkPoint::Rotate(double alpha, double beta, double gamma)
{
   MkMatrix4<double> rm;
   rm.Identity();

   rm.Rotate(alpha, beta, gamma);
   return *this *= rm;
}

MkPoint &MkPoint::RotateInX(double ang)
{
   MkMatrix4<double> rm;
   rm.Identity();

   rm.RotateInX(ang);
   return *this *= rm;
}

MkPoint &MkPoint::RotateInY(double ang)
{
   MkMatrix4<double> rm;
   rm.Identity();

   rm.RotateInY(ang);
   return *this *= rm;
}

MkPoint &MkPoint::RotateInZ(double ang)
{
   MkMatrix4<double> rm;
   rm.Identity();

   rm.RotateInZ(ang);
   return *this *= rm;
}

MkPoint &MkPoint::RotateInA(double ang, double l, double m, double n)
{
   MkMatrix4<double> rm;
   rm.Identity();

   rm.RotateInA(ang, l, m, n);
   return *this *= rm;
}

MkPoint &MkPoint::Scale(double sx, double sy, double sz)
{
   X = X * sx;
   Y = Y * sy;
   Z = Z * sz;
   return *this;
}

void MkPoint::Normalize()
{
   double dist = X * X + Y * Y + Z * Z;
   if (dist < 0.0001)
      return;
   dist = sqrt(dist);
   X = X / dist;
   Y = Y / dist;
   Z = Z / dist;
   return;
}

#ifdef __BCPLUSPLUS__

void MkPoint::Draw(TObject *Sender)
{
   TColor C;
   double Offset;

   if (String(Sender->ClassName()) == String("MkPaintBox"))
   {
      MkPaintBox *pb = (MkPaintBox *)Sender;
      C = pb->Canvas->Pen->Color;
      pb->Canvas->Pen->Color = Color;
      Offset = pb->Offset(3);
      pb->Rectangle3D(X - Offset, Y - Offset, Z - Offset, X + Offset, Y + Offset, Z + Offset);
      pb->Canvas->Pen->Color = C;
   }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)

void MkPoint::Draw(MkPaint *pb)
{
   double sr;
   double DotRadius = 3;
   sr = pb->Offset(DotRadius);
   //  pb->Sphere(X,Y,Z,sr);
}
#endif

#if defined(__GL_H__)

void MkPoint::Draw()
{
}
#endif

void Swap(MkPoint &p1, MkPoint &p2)
{
   MkPoint tp;
   tp = p1;
   p1 = p2;
   p2 = tp;
}

//---------------------------------------------------------------------------

MkPoints::MkPoints(int size, MkPoint *rps)
{

   if (size < 0)
   {
      MkDebug("::MkPoints - MkPoints(int size)");
      ;
      return;
   }

   FSizeOfArray = FSize = size;
   if (FSize == 0)
   {
      FPoint.reset();
      return;
   }

   //   try
   //   {
   //     F.reset(new T[FSize]);
   //   }
   //   catch (std::bad_alloc &a)
   //   {
   //     MkDebug("MkContainer Contructor memory allocation std::bad_alloc thrown\n");
   //     throw Alloc(std::string("MkContainer<T>::bad_alloc error"));
   //   }

   try
   {
      FPoint.reset(new MkPoint[FSize]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkPoints::MkPoints bad_alloc thrown!!!\n");
      throw Alloc(a.what());
   }

   for (int i = 0; i < FSize; i++)
      (*this)[i] = rps[i];
}

MkPoints::MkPoints(int size)
{
   if (size < 0)
   {
      MkDebug("::MkPoints - MkPoints(int size)");
      throw Size(std::string("MkPoints::Constructur size is negative"), size);
   }

   FSize = FSizeOfArray = size;

   if (FSizeOfArray == 0)
   {
      FPoint.reset();
      return;
   }
   try
   {
      FPoint.reset(new MkPoint[FSizeOfArray]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkPoints::MkPoints bad_alloc thrown!!!\n");
      throw Alloc(a.what());
   }
}

MkPoints::~MkPoints()
{
   FSizeOfArray = FSize = 0;
   FPoint.reset();
}

void MkPoints::Initialize(int size)
{
   if (size < 0)
   {
      MkDebug("::MkPoints - Initialize(int size)");
      ;
      return;
   }
   if (FSizeOfArray == size)
      return;

   FSize = FSizeOfArray = size;

   if (FSizeOfArray == 0)
   {
      FPoint.reset();
      return;
   }

   try
   {
      FPoint.reset(new MkPoint[FSizeOfArray]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkPoints::MkPoints bad_alloc thrown!!!\n");
      throw Alloc(a.what());
   }
}

void MkPoints::Initialize(int size, MkPoint *rps)
{

   if (size < 0)
   {
      MkDebug("MkPoints::MkPoints - Initialize(int size) size <0 or rps is null");
      throw Size(std::string("MkPoints::Initialize negative size"), size);
   }
   if (rps == NULL)
   {
      MkDebug("MkPoints::MkPoints - Initialize(int size) size <0 or rps is null");
      throw Size(std::string("MkPoints::Initialize null points argument "), 0);
   }

   if (FSizeOfArray == size)
   {
      FSize = size;
      for (int i = 0; i < FSizeOfArray; i++)
         FPoint[i] = rps[i];
      return;
   }

   FSize = FSizeOfArray = size;
   if (FSizeOfArray == 0)
   {
      FPoint.reset();
      return;
   }

   try
   {
      FPoint.reset(new MkPoint[FSizeOfArray]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkPoints::MkPoints bad_alloc thrown!!!\n");
      throw Alloc(a.what());
   }

   for (int i = 0; i < FSizeOfArray; i++)
      FPoint[i] = rps[i];
}

void MkPoints::FindCenter() // find maximum and minimum for each axis and average it, the center
{
   if (FSize == 0)
   {
      FCenter = NullPoint;
      return;
   }

   double x1, y1, z1, x2, y2, z2;
   x1 = (*this)[0].X;
   y1 = (*this)[0].Y;
   z1 = (*this)[0].Z;
   x2 = (*this)[0].X;
   y2 = (*this)[0].Y;
   z2 = (*this)[0].Z;

   for (int i = 1; i < GetSize(); i++)
   {
      x1 = x1 < (*this)[i].X ? x1 : (*this)[i].X;
      y1 = y1 < (*this)[i].Y ? y1 : (*this)[i].Y;
      z1 = z1 < (*this)[i].Z ? z1 : (*this)[i].Z;
      x2 = x2 > (*this)[i].X ? x2 : (*this)[i].X;
      y2 = y2 > (*this)[i].Y ? y2 : (*this)[i].Y;
      z2 = z2 > (*this)[i].Z ? z2 : (*this)[i].Z;
   }
   FCenter.X = (x1 + x2) / 2.0;
   FCenter.Y = (y1 + y2) / 2.0;
   FCenter.Z = (z1 + z2) / 2.0;
}

int MkPoints::Grow(int delta)
{
   int i;
   boost::shared_array<MkPoint> rp;

   try
   {
      rp.reset(new MkPoint[FSizeOfArray + delta]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkPoint::Grow std::bad_alloc thrown");
      throw Alloc(a.what());
   }

   for (i = 0; i < FSize; i++)
      rp[i] = FPoint[i];
   for (i = FSize; i < FSizeOfArray + delta; i++)
      rp[i] = NullPoint;

   FPoint = rp; // it doesn't leak memory as it is smart pointers (shared_array)
   FSizeOfArray = FSizeOfArray + delta;
   return FSizeOfArray;
}

int MkPoints::Shrink(int delta)
{
   int i;
   boost::shared_array<MkPoint> rp;

   try
   {
      rp.reset(new MkPoint[FSizeOfArray - delta]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkPoint::Shrink std::bad_alloc thrown");
      throw Alloc(a.what());
   }

   for (i = 0; i < FSize; i++)
      rp[i] = FPoint[i];

   for (i = FSize; i < FSizeOfArray - delta; i++)
      rp[i] = NullPoint;

   FPoint = rp; // it doesn't leak memory as it is smart pointers (shared_array)
   FSizeOfArray = FSizeOfArray - delta;
   return FSizeOfArray;
}

bool MkPoints::Add(MkPoint point)
{
   if (FSize == FSizeOfArray)
   {
      try
      {
         Grow(FSize - FSizeOfArray + 10);
      }
      catch (Alloc &a)
      {
         MkDebug("MkPoints::Add single point, Alloc thrown while attempt Grow()");
         throw Alloc(a.what());
      }
   }
   if (FSize == FSizeOfArray)
      return false;

   FSize++;
   FPoint[FSize - 1] = point;
   return true;
}

bool MkPoints::Add(int index, MkPoint point)
{
   if (FSize == FSizeOfArray)

      try
      {
         Grow(FSize - FSizeOfArray + 10);
      }
      catch (Alloc &a)
      {
         MkDebug("MkPoints::Add multiple point, Alloc thrown while attempt Grow()");
         throw Alloc(a.what());
      }

   if (FSize == FSizeOfArray)
      return false;

   for (int i = FSize - 1; i >= index; i--)
      FPoint[i + 1] = FPoint[i];
   FSize++;
   FPoint[index] = point;
   return true;
}

bool MkPoints::Delete(MkPoint point)
{
   int i;
   for (i = 0; i < FSize; i++)
   {
      if (FPoint[i] == point)
         break;
   }
   if (i == FSize)
      return false;
   if (FPoint[i] == point)
   {
      for (int j = i; j < FSize - 1; j++)
         FPoint[j] = FPoint[j + 1];
   }
   FSize--;
   FPoint[FSize] = NullPoint;
   return true;
}

bool MkPoints::Delete(int index)
{
   for (int j = index; j < FSize - 1; j++)
      FPoint[j] = FPoint[j + 1];

   FSize--;
   FPoint[FSize] = NullPoint;
   return true;
}

bool MkPoints::Swap(int i, int j)
{
   MkPoint p;
   if (i >= FSize || j >= FSize)
      return false;
   p = FPoint[i];
   FPoint[i] = FPoint[j];
   FPoint[j] = p;
   return true;
}

bool MkPoints::Clear()
{
   FSizeOfArray = FSize = 0;
   FPoint.reset();
   return true;
}

MkPoint &MkPoints::operator[](int i)
{
   if (FSizeOfArray == 0)
      return NullPoint;

   if (i >= FSize && i < FSizeOfArray)
      FSize = i + 1;

   if (i >= 0 && i < FSize)
      return FPoint[i];
   else
      return NullPoint;
}

MkPoints &MkPoints::operator=(MkPoints &points)
{
   int i;

   Clear();
   FSize = points.FSize;
   FSizeOfArray = points.FSizeOfArray;
   try
   {
      FPoint.reset(new MkPoint[FSizeOfArray]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkPoints::operator= throw std::bad_alloc while reserving memory for copying points\n");
      throw Alloc(a.what());
   }
   for (i = 0; i < FSize; i++)
      FPoint[i] = points.FPoint[i];

   for (i = FSize; i < FSizeOfArray; i++)
      FPoint[i] = NullPoint;

   return *this;
}

bool MkPoints::operator==(MkPoints &points)
{
   int i;

   if (FSize != points.FSize)
      return false;
   for (i = 0; i < FSize; i++)
      if (this->FPoint[i] != points.FPoint[i])
         return false;

   return true;
}

MkPoints &MkPoints::operator*=(MkMatrix4<double> &rm)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i] = this->FPoint[i] * rm;
   return *this;
}

MkPoints &MkPoints::Translate(MkPoint rp)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].Translate(rp);
   return *this;
}

MkPoints &MkPoints::Translate(double x, double y, double z)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].Translate(x, y, z);
   return *this;
}

MkPoints &MkPoints::Rotate(double alpha, double beta, double gamma)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].Rotate(alpha, beta, gamma);
   return *this;
}

MkPoints &MkPoints::RotateInX(double ang)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].RotateInX(ang);
   return *this;
}

MkPoints &MkPoints::RotateInY(double ang)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].RotateInY(ang);
   return *this;
}

MkPoints &MkPoints::RotateInZ(double ang)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].RotateInZ(ang);
   return *this;
}

MkPoints &MkPoints::RotateInA(double ang, double l, double m, double n)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].RotateInA(ang, l, m, n);
   return *this;
}

MkPoints &MkPoints::Scale(double sx, double sy, double sz)
{
   for (int i = 0; i < FSize; i++)
      this->FPoint[i].Scale(sx, sy, sz);
   return *this;
}

#ifdef __BCPLUSPLUS__

void MkPoints::Draw(TObject *Sender)
{
   for (int i = 0; i < FSize; i++)
      (*this)[i].Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)

void MkPoints::Draw(MkPaint *pb)
{
   for (int i = 0; i < FSize; i++)
      (*this)[i].Draw(pb);
}
#endif
//---------------------------------------------------------------------------
