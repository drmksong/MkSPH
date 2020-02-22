//---------------------------------------------------------------------------
#ifndef MkFailureH
#define MkFailureH

#include <stdio.h>
//#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include "MkDouble.h"
#include "MkInt.h"
#include "MkMisc.h"
#include "MkStress.h"

#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif

class MkFailure
{
public:
    double FAlpha;
    double FKay;
    bool FTensionStatus; // true for failure, false for stable
    bool FShearStatus;

public:
    MkFailure();
    MkFailure(double fri, double coh);
    void Set(double fri, double coh);
    void Clear();
    bool CheckTension(MkStress &str);
    bool CheckShear(MkStress &str);
    void TreatTension(MkStress &str);
    void TreatScaling(MkStress &str);

    MkFailure &operator=(const MkFailure &str);

    void Out(char *);
    void Out();
};

extern MkFailure NullFailure;
//---------------------------------------------------------------------------
#endif
