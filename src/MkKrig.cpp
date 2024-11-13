//---------------------------------------------------------------------------
#include "MkKrig.hpp"

MkKrig::MkKrig()
{
   FDataPoints.Clear();
   FX.Clear();
   FM.Clear();
   FR.Clear();
   FSZ.Clear();
   FPZ.Clear();
   FA.Clear();
//   FB.Clear();
//   FL.Clear();
//   FA_bar.Clear();
//   FM_bar.Clear();
//   FR_bar.Clear();
//   FP_bar.Clear();

   FS.Clear();
   FP.Clear();
}

MkKrig::MkKrig(MkDataPoints &rps)
{
   FDataPoints = rps;
   FI = FDataPoints.GetSize();
   FJ = 10;

   FX.Initialize(FI);
   FM.Initialize(FI);
   FR.Initialize(FI);
   FSZ.Initialize(FI);
   FPZ.Initialize(FJ);
   FA.Initialize(FJ);
//   FB.Initialize(FI);
//   FL.Initialize(FI);
//   FA_bar.Initialize(FI);
//   FM_bar.Initialize(FI);
//   FR_bar.Initialize(FI);
//   FP_bar.Clear();

   FS.Initialize(FI,FI);
   FP.Initialize(FI,FJ);

   SetupMatrix();
}


void MkKrig::SetDataPoints(MkDataPoints &rps)
{
   FDataPoints = rps;
   FI = FDataPoints.GetSize();
   FJ = 10;

   FX.Initialize(FI);
   FM.Initialize(FI);
   FR.Initialize(FI);
   FSZ.Initialize(FI);
   FPZ.Initialize(FJ);
   FA.Initialize(FJ);
//   FB.Initialize(FI);
//   FL.Initialize(FI);
//   FA_bar.Initialize(FI);
//   FM_bar.Initialize(FI);
//   FR_bar.Initialize(FI);
//   FP_bar.Clear();

   FS.Initialize(FI,FI);
   FP.Initialize(FI,FJ);

   SetupMatrix();
}

void MkKrig::SetupMatrix()
{
    isSetup = true;
    X();
    S();
    P();
    A();
    M();
    R();
    isSetup = false;
}

MkMatrix<float> &MkKrig::S()
{
    if(!isSetup) return FS;
    MkVector<float> v;
    for (int j=0;j<FI;j++) {
        v = S(FDataPoints[j]);
        for (int i=0;i<FI;i++)
            FS(i,j) = v(i);
    }
    FIS = FS;
    FIS.GetInvert();
    return FS;
}

MkVector<float> &MkKrig::S(MkPoint &rp)
{
    float x,y,z;
    x = rp.X;y = rp.Y;z = rp.Z;
    return S(x,y,z);
}

MkVector<float> &MkKrig::S(float x,float y,float z)
{
    float h;
    for (int i = 0 ; i < FI ;i++){
        MkPoint rp1 = FDataPoints[i];
        MkPoint rp2(x,y,z);
        h = CalDist(rp1,rp2);
        FSZ(i) = Sigma(h);
    }
    return FSZ;
}

MkMatrix<float> &MkKrig::P()
{
    if(!isSetup) return FP;
    MkVector<float> v;
    for (int i=0;i<FI;i++) {
        v = P(FDataPoints[i]);
        for (int j=0;j<FJ;j++)
            FP(i,j) = v(j);
    }
    return FP;
}

MkVector<float> &MkKrig::P(MkPoint &rp)
{
    float x,y,z;
    x = rp.X;y = rp.Y;z = rp.Z;
    return P(x,y,z);
}

MkVector<float> &MkKrig::P(float x,float y,float z)
{
    FPZ(0) = 1;
    FPZ(1) = x;
    FPZ(2) = y;
    FPZ(3) = z;
    FPZ(4) = x*y;
    FPZ(5) = x*z;
    FPZ(6) = y*z;
    FPZ(7) = x*x;
    FPZ(8) = y*y;
    FPZ(9) = z*z;
    return FPZ;
}

MkVector<float> & MkKrig::A()
{
    if(!isSetup) return FA;
    MkMatrix<float> pt,pt1,p,s1,m,s2;
    MkVector<float> x;

    p = P();
    pt = p;
    pt.GetTranspose();
    s1 = FIS;
    s2 = s1*p;
    pt1 = pt*s2;
    pt1.GetInvert();

    s1=FIS;
    pt = p;
    pt.GetTranspose();

    MkVector fa = pt1*pt*s1*X();
    FA = fa;
    return FA;
}

MkVector<float> &MkKrig::M()
{
    if(!isSetup) return FM;
    float a;
    for (int i=0;i<FI;i++) {
        a = M(FDataPoints[i]);
        FM(i) = a;
    }
    return FM;
}

float MkKrig::M(MkPoint &rp)
{
    float x,y,z;
    x = rp.X;y = rp.Y;z = rp.Z;
    return M(x,y,z);
}

float MkKrig::M(float x,float y,float z)
{

//    static MkVector v;
    float result;
    MkVector<float> v;
    v = P(x,y,z);
    result = FA*v;
    return result;   
}

MkVector<float> &MkKrig::R()
{
    if(!isSetup) return FR;
    FR = X();
    FR - M();
    return FR;
}

float MkKrig::R(MkPoint &rp)
{
    float x,y,z;
    x = rp.X;y = rp.Y;z = rp.Z;
    return R(x,y,z);
}

float MkKrig::R(float x,float y,float z)
{
    float result;
    MkVector<float> r,v1,s;
    MkMatrix<float> m;
    r = R();
    s = S(x,y,z);
    m = FIS;
    v1 = m*r;
    result = v1*s;
    return result;
}

MkVector<float> &MkKrig::X()
{
  if(!isSetup) return FX;
   for (int i=0;i < FI;i++)
       FX(i) = FDataPoints[i].Data;
   return FX;
}

float MkKrig::X(MkPoint &rp)
{
    float x,y,z;
    x = rp.X;y = rp.Y;z = rp.Z;
    return X(x,y,z);
}

float MkKrig::X(MkDataPoint &rdp)
{
    float x,y,z;
    x = rdp.X;y = rdp.Y;z = rdp.Z;
    return X(x,y,z);
}

float MkKrig::X(float x,float y,float z)
{
    return M(x,y,z)+R(x,y,z);
}

float MkKrig::Estimate(float x,float y,float z)
{
    return X(x,y,z);
}

float MkKrig::Estimate(MkPoint &rp)
{
    return X(rp);
}

float MkKrig::Estimate(MkDataPoint &rdp)
{
    return X(rdp);
}

//---------------------------------------------------------------------------
