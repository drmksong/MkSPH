//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif

#pragma hdrstop
#include "MkSpline.h"

#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif

MkSpline NullSpline(0);
//---------------------------------------------------------------------------
MkSpline::MkSpline()
{
   YP1 = 0;
   YPn = 0;
   N = 0;
}
//---------------------------------------------------------------------------
MkSpline::MkSpline(int i)
{
   YP1 = 0;
   YPn = 0;
   N = 0;
}
//---------------------------------------------------------------------------
MkSpline::MkSpline(MkFloat &x, MkFloat &y)
{
     YP1 = 0;
     YPn = 0;
     N = 0;

     if (x.getSzX() != y.getSzX()) return;
     N = x.getSzX();
     X.CopyFrom(x);
     Y.CopyFrom(y);

     Y2.Initialize(N);

     CalcY2();
}
//---------------------------------------------------------------------------
MkSpline::MkSpline(MkFloat &x, MkFloat &y, float yp1, float ypn)
{
     if (x.getSzX() != y.getSzX()) return;
     N = x.getSzX();

     X.CopyFrom(x);
     Y.CopyFrom(y);

     YP1 = yp1;
     YPn = ypn;

     Y2.Initialize(N);

     CalcY2();
}
//---------------------------------------------------------------------------
void MkSpline::Set(MkFloat &x, MkFloat &y)
{
     YP1 = 0;
     YPn = 0;
     N = 0;

     if (x.getSzX() != y.getSzX()) return;
     N = x.getSzX();
     X.CopyFrom(x);
     Y.CopyFrom(y);
     Y2.Initialize(N);

     CalcY2();
}
//---------------------------------------------------------------------------
void MkSpline::Set(float yp1,float ypn)
{
     YP1 = yp1;
     YPn = ypn;

     CalcY2();
}
//---------------------------------------------------------------------------
void MkSpline::Set(MkFloat &x, MkFloat &y, float yp1, float ypn)
{
     if (x.getSzX() != y.getSzX()) return;
     N = x.getSzX();

     X.CopyFrom(x);
     Y.CopyFrom(y);
     YP1 = yp1;
     YPn = ypn;
     Y2.Initialize(N);

     CalcY2();
}
//---------------------------------------------------------------------------
void MkSpline::CalcY2()
{
	int i,k;
	float p,qn,sig,un;
     MkFloat u(N);

 	if (YP1 > 0.99e30)
		Y2(0)=u(0)=0.0;
	else {
		Y2(0) = -0.5;
		u(0)=(3.0/(X(1)-X(0)))*((Y(1)-Y(0))/(X(1)-X(0))-YP1);
	}
	for (i=1;i<N-1;i++) {
		sig=(X(i)-X(i-1))/(X(i+1)-X(i-1));
		p=sig*Y2(i-1)+2.0;
		Y2(i)=(sig-1.0)/p;
		u(i)=(Y(i+1)-Y(i))/(X(i+1)-X(i)) - (Y(i)-Y(i-1))/(X(i)-X(i-1));
		u(i)=(6.0*u(i)/(X(i+1)-X(i-1))-sig*u(i-1))/p;
	}
	if (YPn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(X(N-1)-X(N-2)))*(YPn-(Y(N-1)-Y(N-2))/(X(N-1)-X(N-2)));
	}
	Y2(N-1)=(un-qn*u(N-2))/(qn*Y2(N-2)+1.0);
	for (k=N-2;k>=0;k--)
		Y2(k)=Y2(k)*Y2(k+1)+u(k);
}
//---------------------------------------------------------------------------
MkSpline & MkSpline::operator=(MkSpline &spline)
{
   X.CopyFrom(spline.X);
   Y.CopyFrom(spline.Y);
   Y2.CopyFrom(spline.Y2);
   YP1 = spline.YP1;
   YPn = spline.YPn;
   N = spline.N;
   return *this;
}
//---------------------------------------------------------------------------
bool MkSpline::operator==(MkSpline &)
{
  return false; // not yet implemented
}
//---------------------------------------------------------------------------
bool MkSpline::operator!=(MkSpline &)
{
  return true; // not yet implemented
}
//---------------------------------------------------------------------------
float MkSpline::operator[](float x)
{
	int klo,khi,k;
	float h,b,a;
     if (!X.getSzX()) return 0;

	klo=0;
	khi=N-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (X(k) > x) khi=k;
		else klo=k;
	}
	h=X(khi)-X(klo);
	if (h == 0.0) MkDebug("Bad xa input to routine splint");
	a=(X(khi)-x)/h;
	b=(x-X(klo))/h;
	return a*Y(klo)+b*Y(khi)+((a*a*a-a)*Y2(klo)+(b*b*b-b)*Y2(khi))*(h*h)/6.0;
}
//---------------------------------------------------------------------------
float MkSpline::operator()(float x)
{
     return (*this)[x];
}
//---------------------------------------------------------------------------
MkSplines::MkSplines(int size)
{
    if (size < 0) {
      MkDebug("::MkSplines - MkSplines(int size)");;
      return;
    }

    FSize = size;
    if (FSize == 0) {
       FSpline = NULL;
       return;
    }

    FSpline = new MkSpline[FSize];
}

MkSplines::~MkSplines()
{
    FSize = 0;
    if (FSpline) {
       delete[] (MkSpline*)FSpline;
       FSpline = NULL;
    }
}

void MkSplines::Initialize(int size)
{
    if (size < 0) {
      MkDebug("::MkSplines - MkSplines(int size)");;
      return;
    }

    FSize = size;
    if (FSize == 0) {
       FSpline = NULL;
       return;
    }

    if (FSpline) {
       delete[] FSpline;
       FSpline = NULL;
    }
    FSpline = new MkSpline[FSize];
}

bool MkSplines::Clear()
{
    FSize = 0;
    if (FSpline) {
       delete[] FSpline;
       FSpline = NULL;
    }
    return true;
}
MkSpline & MkSplines::operator[](int i)
{
    if (FSize == 0) return NullSpline;
    else if (i >=0 && i < FSize) return FSpline[i];
    else return NullSpline;
}

MkSplines & MkSplines::operator=(MkSplines &splines)
{
    int i;

    Clear();
    FSize = splines.FSize;
    if (FSize == 0) {
       this->FSpline = NULL;
       return *this;
    }
    this->FSpline = new MkSpline[FSize];

    for (i=0;i<FSize;i++)
      this->FSpline[i] = splines.FSpline[i];

    return *this;

}

bool MkSplines::operator==(MkSplines &splines)
{
    int i;

    if (FSize != splines.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FSpline[i] != splines.FSpline[i]) return false;

    return true;
}


#ifdef __BCPLUSPLUS__
void MkSplines::Draw(TObject *)
{

}
#endif

