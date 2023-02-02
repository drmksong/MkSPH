//---------------------------------------------------------------------------
// ::new convienient float vector class
// This can be any dimension up to 3.
// Vector can be re-initialized even though
// previous data will be erased.
// But if you want to preserve it, you can do it simply...
// Author : M. K. Song. Seoul, Korea. 1999.2.20
//---------------------------------------------------------------------------

#include "MkFloat.h"


#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif

MkFloat::MkFloat(int s_x,int s_y,int s_z)
{
    if (s_x<=0 || s_y<=0 || s_z <=0) {
       printf("Three D::bad MkFloat size\n");
       exit(-10);// Size(s_x,s_y,s_z);
    }
    Zero = 0;
    FDimension = 3;
    FI = FJ = FK = 0;

    long sz = long(s_x)*long(s_y)*long(s_z);
    sz_x = long(s_x);
    sz_y = long(s_y);
    sz_z = long(s_z);

    try {F = new float[sz];}

    catch (std::bad_alloc){
       printf( "Not Enough Memory\n");
       exit(-10);// Alloc();
    }
    for(int i=0;i<sz;i++) F[i]=0;
}

MkFloat::MkFloat(int s_x,int s_y)
{
    if (s_x<=0 || s_y<=0) {
       printf( "Two D::bad MkFloat size\n");
       exit(-10);// Size(s_x,s_y);
    }
    Zero = 0;
    FDimension = 2;
    FI = FJ = FK = 0;

    long sz = long(s_x)*long(s_y);
    sz_x = long(s_x);
    sz_y = long(s_y);
    sz_z = 1;

    try {
        F = new float[sz];
    }
    catch (std::bad_alloc){
       printf( "Not Enough Memory\n");
       exit(-10);// Alloc();
    }
    for(int i=0;i<sz;i++) F[i]=0;
}

MkFloat::MkFloat(int s_x)
{
    if (s_x<=0) {
       printf( "One D::bad MkFloat size\n");
       exit(-10);// Size(s_x);
    }

    Zero = 0;
    FDimension = 1;
    FI = FJ = FK = 0;

    int sz = s_x;
    sz_x = long(s_x);
    sz_y = 1;
    sz_z = 1;

    try {
        F = new float[sz];
    }
    catch (std::bad_alloc){
       printf( "Not Enough Memory\n");
       exit(-10);// Alloc();
    }
    for(int i=0;i<sz;i++) F[i]=0;
}

MkFloat::MkFloat()
{
    Zero = 0;
    FDimension = 0;
    FI = FJ = FK = 0;

    sz_x = 0;
    sz_y = 0;
    sz_z = 0;

    F = NULL;
}

MkFloat::~MkFloat()
{
   if (F!=NULL) {
     delete[] F;
     F=NULL;
   }
}

void MkFloat::Clear()
{

    if (F!=NULL) {
       delete[] F;
       F = NULL;
    }
    Zero = 0;
    FDimension = 0;
    FI = FJ = FK = 0;

    sz_x = 0;
    sz_y = 0;
    sz_z = 0;
}

void MkFloat::Initialize(int s_x,int s_y,int s_z)
{
    if (s_x<=0 || s_y<=0 || s_z <=0) {
       printf("bad MkFloat size\n");
       exit(-10);// Size(s_x,s_y,s_z);
    }

    FDimension = 3;
    FI = FJ = FK = 0;

    long sz = long(s_x)*long(s_y)*long(s_z);
    sz_x = long(s_x);
    sz_y = long(s_y);
    sz_z = long(s_z);

    if (F != NULL) {
       delete[] F;
       F = NULL;
    }

    try {
        F = new float[sz];
    }
    catch (std::bad_alloc){
       printf( "Not Enough Memory\n");
       exit(-10);// Alloc();
    }
    for(int i=0;i<sz;i++) F[i]=0;
}

void MkFloat::Initialize(int s_x,int s_y)
{
    if (s_x<=0 || s_y<=0) {
       printf( "bad MkFloat size\n");
       exit(-10);// Size(s_x,s_y);
    }

    FDimension = 2;
    FI = FJ = FK = 0;

    long sz = long(s_x)*long(s_y);
    sz_x = long(s_x);
    sz_y = long(s_y);
    sz_z = 1;

    if (F != (float*)NULL) {
       delete[] F;
       F = NULL;
    }

    try {
        F = new float[sz];
    }
    catch (std::bad_alloc){
       printf( "Not Enough Memory\n");
       exit(-10);// Alloc();
    }
    for(int i=0;i<sz;i++) F[i]=0;
}

void MkFloat::Initialize(int s_x)
{
    if (s_x<=0) {
       printf( "bad MkFloat size\n");
       exit(-10);// Size(s_x);
    }
    FDimension = 1;
    FI = FJ = FK = 0;

    long sz = long(s_x);
    sz_x = long(s_x);
    sz_y = 1;
    sz_z = 1;

    if (F != (float*)NULL) {
       delete[] F;
       F = NULL;
    }

    try {
        F = new float[sz];
    }
    catch (std::bad_alloc){
       printf( "Not Enough Memory\n");
       exit(-10);// Alloc();
    }
    for(int i=0;i<sz;i++) F[i]=0;    
}

void MkFloat::CopyFrom(MkFloat &value)
{
    if (!(value.sz_x*value.sz_y*value.sz_z)) return;
    if(sz_x!=value.sz_x || sz_y!=value.sz_y || sz_z!=value.sz_z) {
      Initialize(value.sz_x,value.sz_y,value.sz_z);
      FDimension = value.FDimension;
    }
    long sz = long(sz_x)*long(sz_y)*long(sz_z);    
    for (long i=0;i<sz;i++)
        F[i] = value.F[i];
}

float &MkFloat::operator()(int i,int j,int k)
{
    if ( FDimension <= 0) {
       printf("Possibly Memory dosen't allocated!\n");
       return Zero;
    }

    if (i<0 || sz_x<=i) {
       printf("MkFloat index out of range\n");
       exit(-10);// Range("sz_x",i);
    }
    if (j<0 || sz_y<=j) {
       printf("MkFloat index out of range\n");
       exit(-10);// Range("sz_y",i);
    }
    if (k<0 || sz_z<=k) {
       printf("MkFloat index out of range\n");
       exit(-10);// Range("sz_z",i);
    }
    return F[long(i)+long(j)*sz_x+long(k)*sz_x*sz_y];
}

float &MkFloat::operator()(int i,int j)
{
    if ( FDimension <= 0) {
       printf("Possibly Memory dosen't allocated!\n");
       return Zero;
    }

    if (i<0 || sz_x<=i) {
       printf("MkFloat index out of range\n");
       exit(-10);// Range("sz_x",i);
    }
    if (j<0 || sz_y<=j) {
       printf("MkFloat index out of range\n");
       exit(-10);// Range("sz_y",i);
    }
    return F[long(i)+long(j)*sz_x];
}

float &MkFloat::operator()(int i)
{
    if ( FDimension <= 0) {
       printf("Possibly Memory dosen't allocated!\n");
       return Zero;
    }

    if (i<0 || sz_x<=i) {
       printf("MkFloat index out of range\n");
       exit(-10);// Range("sz_x",i);
    }
    return F[long(i)];
}

MkFloat & MkFloat::operator+=(MkFloat &a)
{
   if (a.sz_x!=sz_x || a.sz_y!=sz_y || a.sz_z!=sz_z) return *this;

   for (int i=0;i<sz_x;i++)
     for (int j=0;j<sz_y;j++)
       for (int k=0;k<sz_z;k++)
         (*this)(i,j,k) += a(i,j,k);

   return *this;
}

MkFloat & MkFloat::operator-=(MkFloat &a)
{
   if (a.sz_x!=sz_x || a.sz_y!=sz_y || a.sz_z!=sz_z) return *this;

   for (int i=0;i<sz_x;i++)
     for (int j=0;j<sz_y;j++)
       for (int k=0;k<sz_z;k++)
         (*this)(i,j,k) -= a(i,j,k);

   return *this;
}

MkFloat & MkFloat::operator*=(float b)
{
  for (int i=0;i<sz_x;i++)
     for (int j=0;j<sz_y;j++)
       for (int k=0;k<sz_z;k++)
         (*this)(i,j,k)*=b;

  return *this;
}

MkFloat & MkFloat::operator/=(float b)
{
  if(b<FTOL &&b>-FTOL) return *this;
  for (int i=0;i<sz_x;i++)
     for (int j=0;j<sz_y;j++)
       for (int k=0;k<sz_z;k++)
         (*this)(i,j,k)/=b;

  return *this;
}

bool MkFloat::operator==(MkFloat &a)
{
  bool flag=true;
  int i,j,k;
  flag=flag && (FDimension==a.FDimension);
  flag=flag && (sz_x==a.sz_x);
  flag=flag && (sz_y==a.sz_y);
  flag=flag && (sz_z==a.sz_z);
  if(!flag) {
    return flag;
  }
  for(i=0;i<sz_x;i++) {
    for(j=0;j<sz_y;j++) {
      for(k=0;k<sz_z;k++) {
	      flag=flag&&(fabs((*this)(i,j,j)-a(i,j,k))<EPS);
      }
    }
  }
  return flag;
}

bool MkFloat::operator!=(MkFloat &a)
{
  return !(*this==a);
}

MkFloat & operator+(MkFloat &a,MkFloat &b)
{
   static MkFloat c;
   c.Clear();
   if (a.sz_x!=b.sz_x || a.sz_y!=b.sz_y || a.sz_z!=b.sz_z) return c;

   c.Initialize(a.sz_x,a.sz_y,a.sz_z);

   for (int i=0;i<a.sz_x;i++)
     for (int j=0;j<a.sz_y;j++)
       for (int k=0;k<a.sz_z;k++)
         c(i,j,k) = a(i,j,k) + b(i,j,k);

   return c;
}

MkFloat & operator-(MkFloat &a,MkFloat &b)
{
   static MkFloat c;
   c.Clear();
   if (a.sz_x!=b.sz_x || a.sz_y!=b.sz_y || a.sz_z!=b.sz_z) return c;

   c.Initialize(a.sz_x,a.sz_y,a.sz_z);

   for (int i=0;i<a.sz_x;i++)
     for (int j=0;j<a.sz_y;j++)
       for (int k=0;k<a.sz_z;k++)
         c(i,j,k) = a(i,j,k) - b(i,j,k);

   return c;
}

MkFloat & operator*(MkFloat &a,float b)
{
  static MkFloat c;
  c.Initialize(a.sz_x,a.sz_y,a.sz_z);
  for (int i=0;i<a.sz_x;i++)
     for (int j=0;j<a.sz_y;j++)
       for (int k=0;k<a.sz_z;k++)
         c(i,j,k) = a(i,j,k)*b;

  return c;
}

MkFloat & operator/(MkFloat &a,float b)
{
  static MkFloat c;
  if(b<FTOL&&b>-FTOL) return a;
  c.Initialize(a.sz_x,a.sz_y,a.sz_z);
  for (int i=0;i<a.sz_x;i++)
     for (int j=0;j<a.sz_y;j++)
       for (int k=0;k<a.sz_z;k++)
         c(i,j,k) = a(i,j,k)/b;

  return c;
}
//---------------------------------------------------------------------------




