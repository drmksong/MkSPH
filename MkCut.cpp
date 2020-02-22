//---------------------------------------------------------------------------
#include "MkCut.h"
#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif
#pragma hdrstop
//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif

MkCut  NullCut(0);
MkCut::MkCut()
{
  Wall[0]=NULL;
  Wall[1]=NULL;
  className = "MkCut";
}

MkCut::MkCut(int n)
{
  Wall[0]=NULL;
  Wall[1]=NULL;
  className = "MkCut";
}

float MkCut::GetDepth(float x)
{
  MkPolygon &poly = Profile.GetProfile();
  MkPoints pnts;
  float depth;
  int i;

  poly.getCrossWithX(x,pnts);
  if(pnts.GetSize()<=0) return 0;

  depth = pnts[0].Y;
  for(i=1;i<pnts.GetSize();i++) depth = (depth<pnts[i].Y)?depth:pnts[i].Y;
  return depth;
}

#ifdef __BCPLUSPLUS__
bool MkCut::UpdateFrom()
{
  if(!Grid) return false;

  Number =Grid->Cells[1][0].ToInt();

  return true;
}
bool MkCut::UpdateTo()
{
  if(!Grid) return false;

  Grid->Cells[0][0] = "Number";

  Grid->Cells[1][0] = Number;

  return true;
}

void MkCut::Out(TObject *Sender)
{

}
#endif

void MkCut::Out(char *fname)
{
  FILE *fp;
  MkPolygon &cutline = Profile.GetProfile();
  fp = fopen(fname,"a");
  if(!fp) {
    MkDebug(fname); MkDebug(" is not found, so fp is null and return false\n");
    return;
  }
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
//fprintf(fp,"  Cut    G.L.\n");
//fprintf(fp,"  No.    (m) \n");

  fprintf(fp,"  %3d   \n",Number);
  for(int i=0;i<cutline.GetSize();i++) {
    fprintf(fp,"        %8.4f %8.4f\n ",cutline[i].X, cutline[i].Y);
  }

  fclose(fp);
}

bool MkCut::operator==(MkCut &cut)
{
  return Depth==cut.Depth && *Wall[0] == *cut.Wall[0] && *Wall[1] == *cut.Wall[1];

}
bool MkCut::operator!=(MkCut &cut)
{
  return !operator==(cut);
}

MkCut& MkCut::operator=(MkCut &cut)
{
  Depth = cut.Depth;
  Wall[0] = cut.Wall[0];
  Wall[1] = cut.Wall[1];
  Profile = cut.Profile;
  return *this;
}

#ifdef __BCPLUSPLUS__
void  MkCut::Draw(TObject *Sender)
{

}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void  MkCut::Draw(MkPaint *pb)
{

}
#endif
//---------------------------------------------------------------------------
MkCuts::MkCuts(int size,MkCut *cuts)
{
    if (size < 0) {
      MkDebug("::MkCuts - MkCuts(int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0) {
       FCut = NULL;
       return;
    }

    FCut = new MkCut[FSize];
    for (int i=0;i<FSize;i++) (*this)[i] = cuts[i];
}

MkCuts::MkCuts(int size)
{
    if (size < 0) {
      MkDebug("::MkCuts - MkCuts(int size)");;
      return;
    }

    FSizeOfArray = size;
    FSize = 0;
    if (FSizeOfArray == 0) {
       FCut = NULL;
       return;
    }

    FCut = new MkCut[FSizeOfArray];
}

MkCuts::~MkCuts()
{
   FSizeOfArray = FSize = 0;
   if (FCut) {
      delete[] FCut;
      FCut = NULL;
   }
}

void MkCuts::Initialize(int size)
{
    int i;
    if (size < 0) {
      MkDebug("::MkCuts - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;

    MkCut *cut = new MkCut[size];
    if(cut) {
      for(i=0;i<min(size,FSizeOfArray);i++)
        cut[i] = FCut[i];
    }

    FSizeOfArray = FSize =size;

    if (FSizeOfArray == 0) {
       if (FCut!=NULL) delete[] (MkCut*)FCut;
       FCut = NULL;
       return;
    }

    if (FCut!=NULL) delete[] (MkCut*)FCut;
    FCut = (cut!=NULL)?cut:new MkCut[FSizeOfArray];
    for (i=0;i<FSize;i++) FCut[i].Number = i;
}

void MkCuts::Initialize(int size,MkCut *cuts)
{
    int i;
    if (size < 0 || cuts == NULL) {
      MkDebug("::MkCuts - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0) {
       if (FCut!=NULL) delete[] (MkCut*)FCut;
       FCut = NULL;
       return;
    }

    if (FCut!=NULL) delete[] (MkCut*)FCut;
    FCut = new MkCut[FSizeOfArray];
    for (i=0;i<FSizeOfArray;i++) FCut[i] = cuts[i];
    for (i=0;i<FSize;i++) FCut[i].Number = i;
}

int MkCuts::Grow(int delta)
{
    int i;
    MkCut *cut=NULL;

	cut = new MkCut[FSizeOfArray+delta];
    if (!(cut)) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        cut[i] = FCut[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        cut[i] = NullCut;
    if (FCut) {
       delete[] (MkCut*)FCut;
       FCut = NULL;
    }
    FCut = cut;
    FSizeOfArray = FSizeOfArray+delta;
    return FSizeOfArray;
}

int MkCuts::Shrink(int delta)
{
    int i;
    MkCut *cut=NULL;

    if (!(cut = new MkCut[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        cut[i] = FCut[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        cut[i] = NullCut;
    if (FCut) {
       delete[] (MkCut*)FCut;
       FCut = NULL;
    }
    FCut = cut;
    FSizeOfArray = FSizeOfArray-delta;
    return FSizeOfArray;
}

bool MkCuts::Add(MkCut &cut)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }

    FSize++;
    FCut[FSize-1] = cut;
    return true;
}

bool MkCuts::Add(int index, MkCut &cut)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }

    for (int i=FSize-1;i>=index;i--)
      FCut[i+1] = FCut[i];
    FSize++;
    FCut[index] = cut;
    return true;
}

bool MkCuts::Delete(MkCut &cut)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FCut[i] == cut) break;
    }
    if(i==FSize) return false;
    if(FCut[i] == cut) {
      for (int j=i;j<FSize-1;j++)
        FCut[j] = FCut[j+1];
    }
    FSize--;
    FCut[FSize] = NullCut;
    return true;
}

bool MkCuts::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FCut[j] = FCut[j+1];

    FSize--;
    FCut[FSize] = NullCut;
    return true;
}

bool MkCuts::Clear()
{
   FSizeOfArray = FSize = 0;
   if (FCut) {
      delete[] FCut;
      FCut = NULL;
   }
   return true;
}

void MkCuts::Out(char *fname)
{
  FILE *fp;
  fp = fopen(fname,"a");
  if(!fp) {
    return;
  }

  if(FSize==0) return;
    
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
  fprintf(fp,"\n");
  fprintf(fp,"                           <Information of Cut Stage>\n");
  fprintf(fp,"\n");
  fprintf(fp,"  Cut     Dist       G.L.\n");
  fprintf(fp,"  No.     (m)        (m) \n");
  fclose(fp);

  for(int i=0;i<FSize;i++) 
    FCut[i].Out(fname);
}

MkCut & MkCuts::operator[](int i)
{
    if (0<=i && i<FSize) return FCut[i];
    else if(FSize<=i && i<FSizeOfArray) {
      FSize=i+1;
      return FCut[i];
    }
    else if (FSizeOfArray <= i) {
      int size = FSizeOfArray;
      if(Grow(i-FSizeOfArray+5)==size) return NullCut;
      return FCut[i];
    }
    else return NullCut;
}

MkCuts & MkCuts::operator=(MkCuts &cuts)
{
    int i;

    Clear();
    FSize = cuts.FSize;
    FSizeOfArray = cuts.FSizeOfArray;
    if (FSizeOfArray == 0) {
       FCut = NULL;
       return *this;
    }
    this->FCut = new MkCut[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FCut[i] = cuts.FCut[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FCut[i] = NullCut;

    return *this;
}

bool MkCuts::operator==(MkCuts &cuts)
{
    int i;

    if (FSize != cuts.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FCut[i] != cuts.FCut[i]) return false;

    return true;
}

#ifdef __BCPLUSPLUS__
void MkCuts::Draw(TObject *Sender)
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
           FCut[i].Draw(pb);
       }
       pb->Canvas->Pen->Color = C;
    }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void  MkCuts::Draw(MkPaint *pb)
{

}
#endif

//---------------------------------------------------------------------------

