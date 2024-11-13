//---------------------------------------------------------------------------
#ifndef MkStressH
#define MkStressH

#include <stdio.h>
//#include <conio.h>
#include <math.h>
#include <stdlib.h>
// #include "MkVect.hpp"
#include "MkMisc.hpp"
#include "MkMatrix.hpp"

#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif

#define M_PI 3.14159265358979323846

class MkStress
{
public:
    double FStress[3][3];
    double FMeanStress; // P
    double FDeviatoric[3][3];
    double FPrincipal[3];
    double LMN[3][3];
    double FI_1, FI_2, FI_3, FJ_2;
    double FPhi; //
public:
    MkStress();
    MkStress(double stress[3][3]);
    void Identity();
    void LoadIdentity() { Identity(); }
    void Clear();
    void CalcMeanStress();
    void CalcDeviatoric();
    void CalcPhi();
    void CalcPrincipal();
    void CalcLMN();
    void CalcFI_1();
    void CalcFI_2();
    void CalcFI_3();
    void CalcFJ_2();
    void CalcInvariant();
    void Eigen();

    MkStress &operator=(const MkStress &str);
    double &operator()(int i, int j)
    {
        static double Zero;
        Zero = 0;
        if (i < 3 && i >= 0 && j < 3 && j >= 0)
            return FStress[i][j];
        else
            exit(-11);
    }

    void Out(char *);
    void Out();
};

extern MkStress NullStress;
//---------------------------------------------------------------------------
#endif
