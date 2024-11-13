#ifndef MkDOF_H
#define MkDOF_H
//--------------------------------------------------------------------
#include <stdio.h>
#include <assert.h>
#include <string>
#include "MkInt.h"
#include "MkFloat.h"
#ifdef __BCPLUSPLUS__
#include <stdctrls.hpp>
#endif
//--------------------------------------------------------------------
// MkDOF.h
//
// Last Revised. 2004. July. 13.
//--------------------------------------------------------------------
enum MkDOFType {doftNone=0, doftXDis, doftYDis, doftZDis,
                            doftXAng, doftYAng, doftZAng,
                            doftTemp, doftHead, doftHydroPress};
enum MkBNDType {bndtNone=0, bndtFree,bndtFix};

//--------------------------------------------------------------------
class MkDOF {
public:
  MkDOFType DOFType;
  MkBNDType BNDType;
  int       SteerNum;
public:
  MkDOF(){DOFType = doftNone; BNDType = bndtNone;SteerNum=-1;}
  MkDOF(int){DOFType = doftNone; BNDType = bndtNone;SteerNum=-1;}
  MkDOF(MkDOFType doft, MkBNDType bndt) {DOFType = doft; BNDType = bndt;SteerNum=-1;}
  MkDOF(MkDOFType doft) {DOFType = doft;SteerNum=-1;}
  MkDOF(MkBNDType bndt) {BNDType = bndt;SteerNum=-1;}
  ~MkDOF(){}
  void SetType(MkDOFType doft, MkBNDType bndt){DOFType=doft;BNDType=bndt;}
  void SetDOFType(MkDOFType doft){DOFType=doft;}
  void SetBNDType(MkBNDType bndt){BNDType = bndt;}
  void SetSteer(int str_num){SteerNum = str_num;}
  bool isSameKindWith(MkDOF &dof);
  MkDOFType GetDOFType(){return DOFType;}
  MkBNDType GetBNDType(){return BNDType;}
  int GetSteer(){return SteerNum;}
  MkDOF &operator=(MkDOF dof){DOFType=dof.DOFType;BNDType=dof.BNDType;SteerNum=dof.SteerNum;return *this;}
  bool operator==(MkDOF dof){return DOFType==dof.DOFType && BNDType==dof.BNDType && SteerNum==dof.SteerNum;}
  bool operator!=(MkDOF dof){return !(*this==dof);} 
};

class MkDOFs {
protected:
    MkDOF *FDOF;
    int FSize;
public:
    MkDOFs(int Size);
    MkDOFs(){FSize = 0; FDOF = NULL;}
    ~MkDOFs(){if (FDOF) {delete[] FDOF;FDOF = NULL;}}
    void Initialize(int Size);
    void Clear();
    int GetSize(){return FSize;}
    virtual MkDOF & operator[](int);
    MkDOFs & operator=(MkDOFs &dofs);
    bool operator==(MkDOFs &dofs)
      {
        if (FSize!=dofs.FSize) return false;  
        for(int i=0;i<FSize;i++) 
          if(FDOF[i]!=dofs.FDOF[i]) 
            return false; 
        return true;
      }
    bool operator!=(MkDOFs &dofs) {return !(*this==dofs);}
#ifdef __BCPLUSPLUS__
    void Out(TMemo *);
#endif
    void Out(char *);
    void Out();
};
//--------------------------------------------------------------------
extern MkDOF NullDOF;
#endif
