//---------------------------------------------------------------------------
#include "MkEntity.h"

//---------------------------------------------------------------------------
MkEntity NullEntity(0);

void MkEntity::Clear()
{
#ifdef __BCPLUSPLUS__
  Grid=NULL;
  Sheet=NULL;
  Name="";
#else
  Name.clear();
#endif

  Division=0;
  Number=0;
  Depth=0;
  Length=0;
  className = "MkEntity";
}

MkEntities::MkEntities(int size,MkEntity *ents)
{
    if (size < 0) {
      MkDebug("::MkEntities - MkEntities(int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0) {
       FEntity = NULL;
       return;
    }

    FEntity = new MkEntity[FSizeOfArray];
    for (int i=0;i<FSizeOfArray;i++) (*this)[i] = ents[i];
}

MkEntities::MkEntities(int size)
{
    if (size < 0) {
      MkDebug("::MkEntities - MkEntities(int size)");;
      return;
    }

    FSizeOfArray = size;
    FSize = 0;
    if (FSizeOfArray == 0) {
       FEntity = NULL;

       return;
    }

    FEntity = new MkEntity[FSizeOfArray];
}

MkEntities::~MkEntities()
{
   FSizeOfArray = FSize = 0;
   if (FEntity) {
      delete[] FEntity;
      FEntity = NULL;
   }
}

void MkEntities::Initialize(int size)
{
    if (size < 0) {
      MkDebug("::MkEntities - Initialize(int size)");
      return;
    }
    if (FSizeOfArray == size) return;

    FSizeOfArray = size;
    FSize = 0;

    if (FSizeOfArray == 0) {
       if (FEntity!=NULL) delete[] (MkEntity*)FEntity;
       FEntity = NULL;
       return;
    }

    if (FEntity!=NULL) delete[] (MkEntity*)FEntity;
    FEntity = new MkEntity[FSizeOfArray];
}

void MkEntities::Initialize(int size,MkEntity *ents)
{

    if (size < 0 || ents == NULL) {
      MkDebug("::MkEntities - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0) {
       if (FEntity!=NULL) delete[] (MkEntity*)FEntity;
       FEntity = NULL;
       return;
    }

    if (FEntity!=NULL) delete[] (MkEntity*)FEntity;
    FEntity = new MkEntity[FSizeOfArray];
    for (int i=0;i<FSizeOfArray;i++) FEntity[i] = ents[i];
}

int MkEntities::Grow(int delta)
{
    int i;
    MkEntity *ent=NULL;

    if (!(ent = new MkEntity[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        ent[i] = FEntity[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        ent[i] = NullEntity;
    if (FEntity) {
       delete[] (MkEntity*)FEntity;
       FEntity = NULL;
    }
    FEntity = ent;
    FSizeOfArray = FSizeOfArray+delta;
    return FSizeOfArray;
}

int MkEntities::Shrink(int delta)
{
    int i;
    MkEntity *ent=NULL;

    if (!(ent = new MkEntity[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        ent[i] = FEntity[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        ent[i] = NullEntity;
    if (FEntity) {
       delete[] (MkEntity*)FEntity;
       FEntity = NULL;
    }
    FEntity = ent;
    FSizeOfArray = FSizeOfArray-delta;
    return FSizeOfArray;
}

bool MkEntities::Add(MkEntity &ent)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }

    FSize++;
    FEntity[FSize-1] = ent;
    return true;
}

bool MkEntities::Add(int index, MkEntity &ent)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }

    for (int i=FSize-1;i>=index;i--)
      FEntity[i+1] = FEntity[i];
    FSize++;
    FEntity[index] = ent;
    return true;
}

bool MkEntities::Delete(MkEntity &ent)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FEntity[i] == ent) break;
    }
    if(i==FSize) return false;
    if(FEntity[i] == ent) {
      for (int j=i;j<FSize-1;j++)
        FEntity[j] = FEntity[j+1];
    }
    FSize--;
    FEntity[FSize] = NullEntity;
    return true;
}

bool MkEntities::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FEntity[j] = FEntity[j+1];

    FSize--;
    FEntity[FSize] = NullEntity;
    return true;
}

bool MkEntities::Clear()
{
   FSizeOfArray = FSize = 0;
   if (FEntity) {
      delete[] FEntity;
      FEntity = NULL;
   }
   return true;
}

MkEntity & MkEntities::operator[](int i)
{
    if (0<=i && i<FSize) return FEntity[i];
    else if(FSize<=i && i<FSizeOfArray) {
      FSize=i+1;
      return FEntity[i];
    }
    else if (FSizeOfArray <= i) {
      Grow(i-FSizeOfArray+5);
      return NullEntity;
    }
    else return NullEntity;
}

MkEntities & MkEntities::operator=(MkEntities &ents)
{
    int i;

    Clear();
    FSize = ents.FSize;
    FSizeOfArray = ents.FSizeOfArray;
    if (FSizeOfArray == 0) {
       FEntity = NULL;
       return *this;
    }
    this->FEntity = new MkEntity[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FEntity[i] = ents.FEntity[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FEntity[i] = NullEntity;

    return *this;
}

bool MkEntities::operator==(MkEntities &ents)
{
    int i;

    if (FSize != ents.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FEntity[i] != ents.FEntity[i]) return false;

    return true;
}
//---------------------------------------------------------------------------

