//---------------------------------------------------------------------------
#ifndef MkKrigH
#define MkKrigH
#include "MkMatrix.h"
#include "MkPointData.h"

class MkKrig {
private:
    MkDataPoints FDataPoints;
    int FI,FJ;
    MkVector FX,FM,FR;
    MkVector FSZ,FPZ,FA;//,FB,FL,FA_bar,FM_bar,FR_bar,FP_bar;
    MkMatrix FS,FP,FIS;
    bool isSetup;
public:
    MkKrig();
    MkKrig(MkDataPoints &);
    void SetDataPoints(MkDataPoints &);
    void SetupMatrix();
    float Sigma(float h){return h < 100 ? 1 - sqrt(h)/10 : 0;};

    MkMatrix &S();
    MkVector &S(MkPoint &);
    MkVector &S(float x,float y,float z);

    MkMatrix &P();
    MkVector &P(MkPoint &);
    MkVector &P(float x,float y,float z);

    MkVector &M();
    float M(MkPoint &);
    float M(float x,float y,float z);

    MkVector &R();
    float R(MkPoint &);
    float R(float x,float y,float z);

    MkVector &X();
    float X(MkPoint &);
    float X(MkDataPoint &);
    float X(float x,float y,float z);

    MkVector &A();

    float Estimate(float x,float y,float z);
    float Estimate(MkPoint &);
    float Estimate(MkDataPoint &);

    float operator()(float x,float y,float z){return Estimate(x,y,z);}
    float operator()(MkPoint &rp){Estimate(rp);}
    float operator()(MkDataPoint &rp){Estimate(rp);}

};

//---------------------------------------------------------------------------
#endif
