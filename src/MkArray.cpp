#include "MkArray.hpp"

//---------------------------------------------------------------------------
// ::new convienient float vector class
// This can be any dimension up to 3.
// Vector can be re-initialized even though
// previous data will be erased.
// But if you want to preserve it, you can do it simply...
// Author : M. K. Song. Seoul, Korea. 1999.2.2;
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif

template <class T>
MkArray<T>::MkArray(int s_x, int s_y, int s_z)
{

  if (s_x <= 0 || s_y <= 0 || s_z <= 0)
  {
    MkDebug("Three D::bad MkArray size\n");
    // exit(-10);
    throw Size(std::string("MkArray<T>"), s_x, s_y, s_z);
  }
  Zero = 0;
  FDimension = 3;
  FI = FJ = FK = 0;

  long sz = long(s_x) * long(s_y) * long(s_z);
  sz_x = long(s_x);
  sz_y = long(s_y);
  sz_z = long(s_z);

  try
  {
    F.reset(new T[sz]);
    //F = new T[sz];
  }

  catch (std::bad_alloc &a)
  {
    MkDebug("MkArray::MkArray(szx,szy,szz) Not Enough Memory %s\n", a.what());
    throw Alloc(std::string("MkArray<T>"));
  }
  for (int i = 0; i < sz; i++)
    F[i] = 0;
}

template <class T>
MkArray<T>::MkArray(int s_x, int s_y)
{
  if (s_x <= 0 || s_y <= 0)
  {
    MkDebug("Two D::bad MkArray size\n");
    // exit(-10);
    throw Size(std::string("MkArray<T>"), s_x, s_y);
  }
  Zero = 0;
  FDimension = 2;
  FI = FJ = FK = 0;

  long sz = long(s_x) * long(s_y);
  sz_x = long(s_x);
  sz_y = long(s_y);
  sz_z = 1;

  try
  {
    //F = new T[sz];
    F.reset(new T[sz]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkArray::MkArray(szx,szy) Not Enough Memory\n");
    throw Alloc(std::string("MkArray<T>"));
  }
  for (int i = 0; i < sz; i++)
    F[i] = 0;
}

template <class T>
MkArray<T>::MkArray(int s_x)
{
  if (s_x <= 0)
  {
    MkDebug("One D::bad MkArray size\n");
    //exit(-10);
    throw Size(std::string("MkArray<T>"), s_x);
  }

  Zero = 0;
  FDimension = 1;
  FI = FJ = FK = 0;

  int sz = s_x;
  sz_x = long(s_x);
  sz_y = 1;
  sz_z = 1;

  try
  {
    //F = new T[sz];
    F.reset(new T[sz]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkArray::MkArray(szx) Not Enough Memory %s\n", a.what());
    //exit(-10);
    throw Alloc(std::string("MkArray<T>"));
  }
  for (int i = 0; i < sz; i++)
    F[i] = 0;
}

template <class T>
MkArray<T>::MkArray()
{
  Zero = 0;
  FDimension = 0;
  FI = FJ = FK = 0;

  sz_x = 0;
  sz_y = 0;
  sz_z = 0;
}

template <class T>
MkArray<T>::~MkArray()
{
  Clear();
}

template <class T>
void MkArray<T>::Clear()
{
  Zero = 0;
  FDimension = 0;
  FI = FJ = FK = 0;

  sz_x = 0;
  sz_y = 0;
  sz_z = 0;
  F.reset();
}

template <class T>
void MkArray<T>::Initialize(int s_x, int s_y, int s_z)
{
  if (s_x <= 0 || s_y <= 0 || s_z <= 0)
  {
    MkDebug("bad MkArray size\n");
    throw Size(std::string("MkArray<T>"), s_x, s_y, s_z);
  }

  FDimension = 3;
  FI = FJ = FK = 0;

  long sz = long(s_x) * long(s_y) * long(s_z);
  sz_x = long(s_x);
  sz_y = long(s_y);
  sz_z = long(s_z);

  try
  {
    F.reset(new T[sz]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("Not Enough Memory\n");
    throw Alloc(std::string("MkArray<T>"));
  }

  for (int i = 0; i < sz; i++)
    F[i] = 0;
}

template <class T>
void MkArray<T>::Initialize(int s_x, int s_y)
{
  if (s_x <= 0 || s_y <= 0)
  {
    MkDebug("bad MkArray size\n");
    throw Size(std::string("MkArray<T>"), s_x, s_y);
  }

  FDimension = 2;
  FI = FJ = FK = 0;

  long sz = long(s_x) * long(s_y);
  sz_x = long(s_x);
  sz_y = long(s_y);
  sz_z = 1;

  try
  {
    F.reset(new T[sz]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("Not Enough Memory\n");
    throw Alloc(std::string("MkArray<T>"));
  }
  for (int i = 0; i < sz; i++)
    F[i] = 0;
}

template <class T>
void MkArray<T>::Initialize(int s_x)
{
  if (s_x <= 0)
  {
    MkDebug("bad MkArray size\n");
    throw Size(std::string("MkArray<T>"), s_x);
  }
  FDimension = 1;
  FI = FJ = FK = 0;

  long sz = long(s_x);
  sz_x = long(s_x);
  sz_y = 1;
  sz_z = 1;

  try
  {
    F.reset(new T[sz]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("Not Enough Memory\n");
    throw Alloc(std::string("MkArray<T>"));
  }
  for (int i = 0; i < sz; i++)
    F[i] = 0;
}

template <class T>
void MkArray<T>::CopyFrom(MkArray<T> &value)
{
  if (!(value.sz_x * value.sz_y * value.sz_z))
    return;
  if (sz_x != value.sz_x || sz_y != value.sz_y || sz_z != value.sz_z)
  {
    try
    {
      Initialize(value.sz_x, value.sz_y, value.sz_z);
      FDimension = value.FDimension;
    }
    catch (Alloc &a)
    {
      MkDebug("MkArray::CopyFrom thows Alloc()");
      throw Alloc(std::string("MkArray<T>"));
    }
  }
  long sz = long(sz_x) * long(sz_y) * long(sz_z);
  for (long i = 0; i < sz; i++)
    F[i] = value.F[i];
}

template <class T>
T &MkArray<T>::operator()(int i, int j, int k)
{
  if (FDimension <= 0)
  {
    MkDebug("Possibly Memory dosen't allocated!\n");
    throw Alloc(std::string("MkArray<T>"));
  }

  if (i < 0 || sz_x <= i)
  {
    MkDebug("MkArray index out of range\n");
    throw Range(std::string("sz_x"), i);
  }
  if (j < 0 || sz_y <= j)
  {
    MkDebug("MkArray index out of range\n");
    throw Range(std::string("sz_y"), j);
  }
  if (k < 0 || sz_z <= k)
  {
    MkDebug("MkArray index out of range\n");
    throw Range(std::string("sz_z"), k);
  }
  return F[long(i) + long(j) * sz_x + long(k) * sz_x * sz_y];
}

template <class T>
T &MkArray<T>::operator()(int i, int j)
{
  if (FDimension <= 0)
  {
    MkDebug("Possibly Memory dosen't allocated!\n");
    throw Alloc(std::string("MkArray<T>"));
  }

  if (i < 0 || sz_x <= i)
  {
    MkDebug("MkArray index out of range\n");
    throw Range("sz_x", i);
  }
  if (j < 0 || sz_y <= j)
  {
    MkDebug("MkArray index out of range\n");
    throw Range(std::string("sz_y"), j);
  }
  return F[long(i) + long(j) * sz_x];
}

template <class T>
T &MkArray<T>::operator()(int i)
{
  if (FDimension <= 0)
  {
    MkDebug("Possibly Memory dosen't allocated!\n");
    throw Alloc(std::string("MkArray<T>"));
  }

  if (i < 0 || sz_x <= i)
  {
    MkDebug("MkArray index out of range\n");
    throw Range(std::string("sz_x"), i);
  }
  return F[long(i)];
}

template <class T>
MkArray<T> &MkArray<T>::operator+=(MkArray &a)
{
  if (a.sz_x != sz_x || a.sz_y != sz_y || a.sz_z != sz_z)
    return *this;

  for (int i = 0; i < sz_x; i++)
    for (int j = 0; j < sz_y; j++)
      for (int k = 0; k < sz_z; k++)
        (*this)(i, j, k) += a(i, j, k);

  return *this;
}

template <class T>
MkArray<T> &MkArray<T>::operator-=(MkArray &a)
{
  if (a.sz_x != sz_x || a.sz_y != sz_y || a.sz_z != sz_z)
    return *this;

  for (int i = 0; i < sz_x; i++)
    for (int j = 0; j < sz_y; j++)
      for (int k = 0; k < sz_z; k++)
        (*this)(i, j, k) -= a(i, j, k);

  return *this;
}

template <class T>
MkArray<T> &MkArray<T>::operator*=(T b)
{
  for (int i = 0; i < sz_x; i++)
    for (int j = 0; j < sz_y; j++)
      for (int k = 0; k < sz_z; k++)
        (*this)(i, j, k) *= b;

  return *this;
}

template <class T>
MkArray<T> &MkArray<T>::operator/=(T b)
{
  if (b < FTOL && b > -FTOL)
    return *this;
  for (int i = 0; i < sz_x; i++)
    for (int j = 0; j < sz_y; j++)
      for (int k = 0; k < sz_z; k++)
        (*this)(i, j, k) /= b;

  return *this;
}

template <class T>
bool MkArray<T>::operator==(MkArray &a)
{
  bool flag = true;
  int i, j, k;
  flag = flag && (FDimension == a.FDimension);
  flag = flag && (sz_x == a.sz_x);
  flag = flag && (sz_y == a.sz_y);
  flag = flag && (sz_z == a.sz_z);
  if (!flag)
  {
    return flag;
  }
  for (i = 0; i < sz_x; i++)
  {
    for (j = 0; j < sz_y; j++)
    {
      for (k = 0; k < sz_z; k++)
      {
        flag = flag && (fabs((*this)(i, j, j) - a(i, j, k)) < EPS);
      }
    }
  }
  return flag;
}

template <class T>
bool MkArray<T>::operator!=(MkArray &a)
{
  return !(*this == a);
}

// template <class T>
// MkArray<T> &operator+(MkArray<T> &a, MkArray<T> &b)
// {
//   static MkArray<T> c;
//   c.Clear();
//   if (a.sz_x != b.sz_x || a.sz_y != b.sz_y || a.sz_z != b.sz_z)
//     return c;

//   c.Initialize(a.sz_x, a.sz_y, a.sz_z);

//   for (int i = 0; i < a.sz_x; i++)
//     for (int j = 0; j < a.sz_y; j++)
//       for (int k = 0; k < a.sz_z; k++)
//         c(i, j, k) = a(i, j, k) + b(i, j, k);

//   return c;
// }

// template <class T>
// MkArray<T> &operator-(MkArray<T> &a, MkArray<T> &b)
// {
//   static MkArray<T> c;
//   c.Clear();
//   if (a.sz_x != b.sz_x || a.sz_y != b.sz_y || a.sz_z != b.sz_z)
//     return c;

//   c.Initialize(a.sz_x, a.sz_y, a.sz_z);

//   for (int i = 0; i < a.sz_x; i++)
//     for (int j = 0; j < a.sz_y; j++)
//       for (int k = 0; k < a.sz_z; k++)
//         c(i, j, k) = a(i, j, k) - b(i, j, k);

//   return c;
// }

// template <class T>
// MkArray<T> &operator*(MkArray<T> &a, T b)
// {
//   static MkArray<T> c;
//   c.Initialize(a.sz_x, a.sz_y, a.sz_z);
//   for (int i = 0; i < a.sz_x; i++)
//     for (int j = 0; j < a.sz_y; j++)
//       for (int k = 0; k < a.sz_z; k++)
//         c(i, j, k) = a(i, j, k) * b;

//   return c;
// }

// template <class T>
// MkArray<T> &operator/(MkArray<T> &a, T b)
// {
//   static MkArray<T> c;
//   if (b < FTOL && b > -FTOL)
//     return a;
//   c.Initialize(a.sz_x, a.sz_y, a.sz_z);
//   for (int i = 0; i < a.sz_x; i++)
//     for (int j = 0; j < a.sz_y; j++)
//       for (int k = 0; k < a.sz_z; k++)
//         c(i, j, k) = a(i, j, k) / b;

//   return c;
// }
//---------------------------------------------------------------------------
