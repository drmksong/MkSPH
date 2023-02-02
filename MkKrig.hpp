//---------------------------------------------------------------------------
#ifndef MkKrigH
#define MkKrigH
#include "MkMatrix.hpp"
#include "MkPointData.hpp"

class MkKrig {
private:
    MkDataPoints FDataPoints;
    int FI,FJ;
    MkVector<float> FX,FM,FR;
    MkVector<float> FSZ,FPZ,FA;//,FB,FL,FA_bar,FM_bar,FR_bar,FP_bar;
    MkMatrix<float> FS,FP,FIS;
    bool isSetup;
public:
    MkKrig();
    MkKrig(MkDataPoints &);
    void SetDataPoints(MkDataPoints &);
    void SetupMatrix();
    float Sigma(float h){return h < 100 ? 1 - sqrt(h)/10 : 0;};

    MkMatrix<float> &S();
    MkVector<float> &S(MkPoint &);
    MkVector<float> &S(float x,float y,float z);

    MkMatrix<float> &P();
    MkVector<float> &P(MkPoint &);
    MkVector<float> &P(float x,float y,float z);

    MkVector<float> &M();
    float M(MkPoint &);
    float M(float x,float y,float z);

    MkVector<float> &R();
    float R(MkPoint &);
    float R(float x,float y,float z);

    MkVector<float> &X();
    float X(MkPoint &);
    float X(MkDataPoint &);
    float X(float x,float y,float z);

    MkVector<float> &A();

    float Estimate(float x,float y,float z);
    float Estimate(MkPoint &);
    float Estimate(MkDataPoint &);

    float operator()(float x,float y,float z){return Estimate(x,y,z);}
    float operator()(MkPoint &rp){Estimate(rp);}
    float operator()(MkDataPoint &rp){Estimate(rp);}

};

//---------------------------------------------------------------------------
#endif
