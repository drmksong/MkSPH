//---------------------------------------------------------------------------
#ifndef MkMiscH
#define MkMiscH
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef __BCPLUSPLUS__
#include <system.hpp>
#include <graphics.hpp>
#endif

#define FTOL 1.0e-6
#define EPS 1.0e-4
#define TINY 1.0e-20
#define PI 3.1415926535

//#define MPa2KPa 1.0e3
//#define KPa2MPa 1.0e-3
//#define MPa2Tonf 1.0e2
//#define Tonf2MPa 1.0e-2
//#define KPa2Tonf 1.0e-1
//#define Tonf2KPa 1.0e+1

#define MPa2KPa 1
#define KPa2MPa 1
#define MPa2Tonf 1
#define Tonf2MPa 1
#define KPa2Tonf 1
#define Tonf2KPa 1

#ifndef MAX_VAL
#define MAX_VAL 1.0e+6
#endif

#ifdef MKDEBUG
#ifdef __BCPLUSPLUS__
#define MkDebug ShowMessage
#else
#define MkDebug printf
#endif
#else
#define MkDebug dumprintf
#endif

enum VectType {vtRow, vtCol,vtNone};
enum MatType { mtNormal, mtLUDecomposed, mtBackSubstitute, mtInverted,
               mtTransposed, mtGauss, mtRow, mtCol};
enum SolveType {stLUD,stGauss,stHybrid};
enum MkSteelType {stHBEAM, stIBEAM, stCHANNEL, stANGLE, stSHEETPILE};

bool is_eq(double i,double j);

#ifndef min
int min(int x,int y);
double min(double x,double y);
#endif
#ifndef max
int max(int x,int y);
double max(double x,double y);
#endif

#if !defined(_MSC_VER) && !defined(_WINDOWS_)
void dull(const char *__format, ...);
#endif

void dumprintf(const char *__format, ...);

void swap(int &x, int &y);
void swap(double &x, double &y);
#ifdef __BCPLUSPLUS__
void swap(TColor a, TColor b);
void Swap(TObject *Sender, int i, int j);
#endif

int delta(int a, int b);
bool ExtractFileExt(char *ext,char *str);   // ext must have a memory
bool TrimLeft(char *&dest, char *src);
bool ExtractStr(char *des,int n,char *src);
bool ToLower(char *str);
bool ToUpper(char *str);

#ifdef __BCPLUSPLUS__
bool ToOnlyAlpha(AnsiString &dest, AnsiString &src);
bool RemoveAnd(AnsiString &dest, AnsiString &src);
#else
bool ToOnlyAlpha(char *&dest, char *src);
bool RemoveAnd(char *&dest, char *src);
#endif

bool CompSub(char *str,char *txt);
int  NumOfParam(char *str);
int  NumOfParen(char *str);
bool ExtractFromParen(char *str, int n, double &x, double &y);

#ifdef __BCPLUSPLUS__
bool IsNumber(AnsiString str);
AnsiString ShortSteelName(AnsiString str);
#else
char *ShortSteelName(char* str);
#endif

bool IsNumber(char *str);

double ShapeFun1(double x, double l);
double ShapeFun2(double x, double l);
double ShapeFun3(double x, double l);
double ShapeFun4(double x, double l);

double ShapeFunInteg1(double aj, double bj, double l, double lj1, double lj);
double ShapeFunInteg2(double aj, double bj, double l, double lj1, double lj);
double ShapeFunInteg3(double aj, double bj, double l, double lj1, double lj);
double ShapeFunInteg4(double aj, double bj, double l, double lj1, double lj);

struct MkOrient {
   double DipDir;
   double Dip;
   MkOrient(){DipDir = 0; Dip = 0;}
   MkOrient(double dip,double dd){Dip=dip;DipDir=dd;}
   MkOrient & operator=(MkOrient &jn){DipDir = jn.DipDir;Dip = jn.Dip;return *this;}
   bool operator==(MkOrient &fo){return (DipDir == fo.DipDir) && (Dip == fo.Dip);}
};

double sosu1(double a);
double sosu2(double a);
double sosu3(double a);
//---------------------------------------------------------------------------
#endif
