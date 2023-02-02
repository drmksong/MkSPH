//---------------------------------------------------------------------------
#ifndef MkFloatH
#define MkFloatH

#undef _MSC_EXTENSIONS

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "MkMisc.hpp"
//---------------------------------------------------------------------------
class MkFloat  {
private:
      int FDimension;
      int FI,FJ,FK;
      float Zero;
      float *F;
      long sz_x,sz_y,sz_z;

public:
       MkFloat(int,int,int);
       MkFloat(int,int);
       MkFloat(int);
       MkFloat();
       ~MkFloat();
       void Clear();
       void Initialize(int s_x,int s_y,int s_z);
       void Initialize(int s_x,int s_y);
       void Initialize(int s_x);
       void CopyFrom(MkFloat &value);
       
       float & operator()(int,int,int);
       float & operator()(int,int);
       float & operator()(int);
       MkFloat & operator=(MkFloat &a){CopyFrom(a); return *this;}

       MkFloat & operator+=(MkFloat &a);
       MkFloat & operator-=(MkFloat &a);
       MkFloat & operator*=(float a);
       MkFloat & operator/=(float a);
       
       bool operator==(MkFloat &a);
       bool operator!=(MkFloat &a);

       friend MkFloat & operator+(MkFloat &a,MkFloat &b);
       friend MkFloat & operator-(MkFloat &a,MkFloat &b);
       friend MkFloat & operator*(MkFloat &a,float b);
       friend MkFloat & operator/(MkFloat &a,float b);

       long getSzX() { return sz_x;}
       long getSzY() { return sz_y;}
       long getSzZ() { return sz_z;}

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
