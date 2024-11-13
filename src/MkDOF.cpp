//--------------------------------------------------------------------
// MkDOF.cpp
// 
// Last Revised. 2004. July. 13
//--------------------------------------------------------------------
#include "MkDOF.h"

MkDOF NullDOF(0);

std::string chDOFType[] = {"doftNone", "doftXDis", "doftYDis", "doftZDis",
                            "doftXAng", "doftYAng", "doftZAng",
                            "doftTemp", "doftHead", "doftHydroPress"};
std::string chBNDType[] = {"bndtNone", "bndtFree", "bndtFix"};


#ifdef __BCPLUSPLUS__
#endif

bool MkDOF::isSameKindWith(MkDOF &dof)
{
  bool b1,b2;
  MkDOFType doftype;
  doftype = dof.DOFType;

  if(DOFType==doftype) return true;
  b1 = DOFType == doftXDis || DOFType == doftYDis || DOFType == doftZDis;
  b2 = doftype == doftXDis || doftype == doftYDis || doftype == doftZDis;
  if(b1&&b2) return true;

  b1 = DOFType == doftXAng || DOFType == doftYAng || DOFType == doftZAng;
  b2 = doftype == doftXAng || doftype == doftYAng || doftype == doftZAng;
  if(b1&&b2) return true;
  
  return false;
}
//--------------------------------------------------------------------
MkDOFs::MkDOFs(int size)
{

    if (size <= 0) {
        MkDebug("::MkDOFs - MkDOFs(int size)");
        return;
    }

    FSize = size;
    FDOF = new MkDOF[FSize];
    assert(FDOF);
}

void MkDOFs::Initialize(int size)
{
    Clear();
    FSize = size;
    if (FSize == 0) {
       FDOF = NULL;
       return;
    }
    FDOF = new MkDOF[FSize];
    assert(FDOF);
}

void MkDOFs::Clear()
{
    FSize = 0;
    if (FDOF) delete[] FDOF;
    FDOF = NULL;
}

MkDOF &  MkDOFs::operator[](int i)
{
    if (i >=0 && i < FSize) return FDOF[i];
    else return NullDOF;
}

MkDOFs &  MkDOFs::operator=(MkDOFs &dofs)
{
    int i;
    FSize = dofs.FSize;
    FDOF = new MkDOF[FSize];

    for (i=0;i<FSize;i++)
      FDOF[i] = dofs.FDOF[i];

    return *this;
}

#ifdef __BCPLUSPLUS__
void MkDOFs::Out(TMemo *memo)
{
  char str[256];
  memo->Lines->Add("DOF Output");
  for (int i=0;i<FSize;i++) {
    sprintf(str,"%d-th DOF is %d ",i,FDOF[i].SteerNum);
    memo->Lines->Add(str);
  }
}
#endif
void MkDOFs::Out(char *fname)
{
  FILE *fp;
  char str[256];

  fp = fopen(fname,"a");
  if(!fp) return;

  fprintf(fp,"DOF Output");
  for (int i=0;i<FSize;i++) {
    sprintf(str,"%d-th DOF is %d ",i,FDOF[i].SteerNum);
    fprintf(fp,str);
  }
  fclose(fp);
}

void MkDOFs::Out()
{
  char str[256];
  MkDebug("DOF Output\n");
  for (int i=0;i<FSize;i++) {
    sprintf(str,"%d-th DOF : (DOFType, %s), (BNDType, %s),(SteerNum, %d)\n",
                  i,chDOFType[FDOF[i].DOFType].c_str(),chBNDType[FDOF[i].BNDType].c_str(),FDOF[i].SteerNum);
    MkDebug(str);
  }
}

//--------------------------------------------------------------------
