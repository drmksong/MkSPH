//---------------------------------------------------------------------------
#ifndef MkIntH
#define MkIntH

#undef _MSC_EXTENSIONS

#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "MkMisc.hpp"
//---------------------------------------------------------------------------
class MkInt  {
private:
      int FDimension;
      int FI,FJ,FK;
      int Zero;
      int  *F;
      long sz_x,sz_y,sz_z;

public:
       MkInt(int,int,int);
       MkInt(int,int);
       MkInt(int);
       MkInt();
       ~MkInt();
       void Clear();
       void Initialize(int s_x,int s_y,int s_z);
       void Initialize(int s_x,int s_y);
       void Initialize(int s_x);
       void CopyFrom(MkInt &value);

       int & operator()(int,int,int);
       int & operator()(int,int);
       int & operator()(int);
       int & operator[](int i){return operator()(i);}
       MkInt & operator=(MkInt &a){CopyFrom(a); return *this;}

       MkInt & operator+=(MkInt &a);
       MkInt & operator-=(MkInt &a);
       MkInt & operator*=(int a);
       MkInt & operator/=(int a);

       bool operator==(MkInt &a);
       bool operator!=(MkInt &a);

       friend MkInt & operator+(MkInt &a,MkInt &b);
       friend MkInt & operator-(MkInt &a,MkInt &b);
       friend MkInt & operator*(MkInt &a,int b);
       friend MkInt & operator/(MkInt &a,int b);

       int & elem( int i, int j , int k ) { return F[(i*sz_y+j)*sz_z+k]; }
                                                 // first k is changed
                                                 // second j is changed
                                                 // third k is changed
       long getSzX() { return sz_x ; };
       long getSzY() { return sz_y ; };
       long getSzZ() { return sz_z ; };

       class Alloc {};
       class Size {
         public:
         int X,Y,Z;
         Size(int x,int y,int z) : X(x), Y(y), Z(z) {}
         Size(int x,int y) : X(x), Y(y), Z(1) {}
         Size(int x) : X(x),Y(1),Z(1) {}
       };
       class Range {
         public:
         char *Str;
         int N;
         Range(char *str,int n) : Str(str), N(n) {}
       };
};

#endif
