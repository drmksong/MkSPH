//---------------------------------------------------------------------------
#include "MkBndCon.hpp"

#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif

MkBndCond NullBndCond(0);
//--------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif
//--------------------------------------------------------------------

MkBndCond::MkBndCond()
{
}

MkBndCond::MkBndCond(int n)
{
}

MkBndCond::~MkBndCond()
{
}

MkBndConds::MkBndConds(int size,MkBndCond *bnds)
{

    if (size < 0) {
      MkDebug("::MkBndConds - MkBndConds(int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0) {
       FBndCond = NULL;
       return;
    }

    FBndCond = new MkBndCond[FSize];
    for (int i=0;i<FSize;i++) (*this)[i] = bnds[i];
}

MkBndConds::MkBndConds(int size)
{
    if (size < 0) {
      MkDebug("::MkBndConds - MkBndConds(int size)");;
      return;
    }

    FSizeOfArray = size;
    FSize = 0;
    if (FSizeOfArray == 0) {
       FBndCond = NULL;
       return;
    }

    FBndCond = new MkBndCond[FSizeOfArray];
}

MkBndConds::~MkBndConds()
{
   FSizeOfArray = FSize = 0;
   if (FBndCond) {
      delete[] (MkBndCond*)FBndCond;
      FBndCond = NULL;
   }
}

void MkBndConds::Initialize(int size)
{
    if (size < 0) {
      MkDebug("::MkBndConds - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;

    FSizeOfArray = size;
    FSize = 0;
    
    if (FSizeOfArray == 0) {
       if (FBndCond!=NULL) delete[] (MkBndCond*)FBndCond;
       FBndCond = NULL;
       return;
    }

    if (FBndCond!=NULL) delete[] (MkBndCond*)FBndCond;
    FBndCond = new MkBndCond[FSizeOfArray];
}

void MkBndConds::Initialize(int size,MkBndCond *bnds)
{

    if (size < 0 || bnds == NULL) {
      MkDebug("::MkBndConds - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0) {
       if (FBndCond!=NULL) delete[] (MkBndCond*)FBndCond;
       FBndCond = NULL;
       return;
    }

    if (FBndCond!=NULL) delete[] (MkBndCond*)FBndCond;
    FBndCond = new MkBndCond[FSizeOfArray];
    for (int i=0;i<FSizeOfArray;i++) FBndCond[i] = bnds[i];
}

int MkBndConds::Grow(int delta)
{
    int i;
    MkBndCond *bndcond=NULL;

    if (!(bndcond = new MkBndCond[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        bndcond[i] = FBndCond[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        bndcond[i] = NullBndCond;
    if (FBndCond) {
       delete[] (MkBndCond*)FBndCond;
       FBndCond = NULL;
    }
    FBndCond = bndcond;
    FSizeOfArray = FSizeOfArray+delta;
    return FSizeOfArray;
}

int MkBndConds::Shrink(int delta)
{
    int i;
    MkBndCond *bndcond=NULL;

    if (!(bndcond = new MkBndCond[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        bndcond[i] = FBndCond[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        bndcond[i] = NullBndCond;
    if (FBndCond) {
       delete[] (MkBndCond*)FBndCond;
       FBndCond = NULL;
    }
    FBndCond = bndcond;
    FSizeOfArray = FSizeOfArray-delta;
    return FSizeOfArray;
}

bool MkBndConds::Add(MkBndCond &bndcond)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) Grow(FSize-FSizeOfArray+1);
    if(tmp==FSizeOfArray) return false;

    FSize++;
    FBndCond[FSize-1] = bndcond;
    return true;
}

bool MkBndConds::Add(int index, MkBndCond &bndcond)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) Grow(FSize-FSizeOfArray+1);
    if(tmp==FSizeOfArray) return false;

    for (int i=FSize-1;i>=index;i--)
      FBndCond[i+1] = FBndCond[i];
    FSize++;
    FBndCond[index] = bndcond;
    return true;
}

bool MkBndConds::Delete(MkBndCond &bndcond)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FBndCond[i] == bndcond) break;
    }
    if(i==FSize) return false;
    if(FBndCond[i] == bndcond) {
      for (int j=i;j<FSize-1;j++)
        FBndCond[j] = FBndCond[j+1];
    }
    FSize--;
    FBndCond[FSize] = NullBndCond;
    return true;
}

bool MkBndConds::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FBndCond[j] = FBndCond[j+1];

    FSize--;
    FBndCond[FSize] = NullBndCond;
    return true;
}

bool MkBndConds::Clear()
{
   FSizeOfArray = FSize = 0;
   if (FBndCond) {
      delete[] FBndCond;
      FBndCond = NULL;
   }
   return true;
}

MkBndCond & MkBndConds::operator[](int i)
{
    if (FSizeOfArray == 0) return NullBndCond;
    if (0<=i && i<FSizeOfArray) FSize = i+1;

    if (i >=0 && i < FSize) return FBndCond[i];
    else return NullBndCond;
}

MkBndConds & MkBndConds::operator=(MkBndConds &bnds)
{
    int i;

    Clear();
    FSize = bnds.FSize;
    FSizeOfArray = bnds.FSizeOfArray;
    if (FSize == 0) {
       FBndCond = NULL;
       return *this;
    }
   
    FBndCond = new MkBndCond[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FBndCond[i] = bnds.FBndCond[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FBndCond[i] = NullBndCond;

    return *this;
}

bool MkBndConds::operator==(MkBndConds &bnds)
{
    int i;

    if (FSize != bnds.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FBndCond[i] != bnds.FBndCond[i]) return false;

    return true;
}

#ifdef __BCPLUSPLUS__
void MkBndConds::Draw(TObject *Sender)
{
    TColor C;
    float Offset;
    if (FSize == 0) return;
    if (String(Sender->ClassName()) == String("MkPaintBox")) {
       MkPaintBox *pb=(MkPaintBox*)Sender;
       C = pb->Canvas->Pen->Color;
       pb->Canvas->Pen->Color = Color;
       for (int i = 0; i < FSize;i++) {
           Offset = pb->Offset(3);
           FBndCond[i].Draw(pb);
       }
       pb->Canvas->Pen->Color = C;
    }
}
#endif

//---------------------------------------------------------------------------

