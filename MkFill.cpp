//---------------------------------------------------------------------------
#include "MkFill.hpp"
#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif
#pragma hdrstop
//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif
//---------------------------------------------------------------------------
MkFill NullFill(0);
MkFill::MkFill()
{
  Depth=0;
  Wall[0]=NULL;
  Wall[1]=NULL;
  className = "MkFill";
}

MkFill::MkFill(int n)
{
  Depth=0;
  Wall[0]=NULL;
  Wall[1]=NULL;
  className = "MkFill";
}

float MkFill::GetDepth(float x)
{
  MkPolygon &poly = Profile.GetProfile();
  MkPoints pnts;
  float depth;
  int i;

  poly.getCrossWithX(x,pnts);
  if(pnts.GetSize()<=0) return 0;

  depth = pnts[0].Y;
  for(i=1;i<pnts.GetSize();i++) depth = (depth>pnts[i].Y)?depth:pnts[i].Y;
  return depth;
}

#ifdef __BCPLUSPLUS__
bool MkFill::UpdateFrom()
{
  if(!Grid) return false;

  Number = Grid->Cells[1][0].ToInt();

  return true;
}
bool MkFill::UpdateTo()
{
  if(!Grid) return false;

  Grid->Cells[0][0] = "Number";
  Grid->Cells[0][1] = "Ground Level";

  Grid->Cells[1][0] = Number;

  return true;
}


void MkFill::Out(TObject *Sender)
{
}
#else

#endif

void MkFill::Out(char *fname)
{
  FILE *fp;
  MkPolygon &fillline = Profile.GetProfile();
  fp = fopen(fname,"a");
  if(!fp) {
    MkDebug(fname); MkDebug(" is not found, so fp is null and return false\n");
    return ;
  }
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
//fprintf(fp,"  Fill    G.L.\n");
//fprintf(fp,"   No.   (m) \n");
  fprintf(fp,"  %3d   \n",Number);
  for(int i=0;i<fillline.GetSize();i++) {
    fprintf(fp,"        %8.4f %8.4f\n ",fillline[i].X, fillline[i].Y);
  }

  fclose(fp);
}

bool MkFill::operator==(MkFill &fill)
{
  return Depth==fill.Depth && *Wall[0] == *fill.Wall[0] && *Wall[1] == *fill.Wall[1];
}
bool MkFill::operator!=(MkFill &fill)
{
  return !operator==(fill);
}

MkFill& MkFill::operator=(MkFill& fill)
{
  Depth = fill.Depth;
  Wall[0] = fill.Wall[0];
  Wall[1] = fill.Wall[1];
  Profile = fill.Profile;
  return *this;
}

#ifdef __BCPLUSPLUS__
void  MkFill::Draw(TObject *Sender)
{

}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void  MkFill::Draw(MkPaint *pb)
{

}
#endif
//---------------------------------------------------------------------------
MkFills::MkFills(int size,MkFill *fills)
{
    if (size < 0) {
      MkDebug("::MkFills - MkFills(int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0) {
       FFill = NULL;
       return;
    }

    FFill = new MkFill[FSize];
    for (int i=0;i<FSize;i++) (*this)[i] = fills[i];
}

MkFills::MkFills(int size)
{
    if (size < 0) {
      MkDebug("::MkFills - MkFills(int size)");;
      return;
    }

    FSizeOfArray = size;
    FSize = 0;
    if (FSizeOfArray == 0) {
       FFill = NULL;
       return;
    }

    FFill = new MkFill[FSizeOfArray];
}

MkFills::~MkFills()
{
   FSizeOfArray = FSize = 0;
   if (FFill) {
      delete[] FFill;
      FFill = NULL;
   }
}

void MkFills::Initialize(int size)
{
    int i;
    if (size < 0) {
      MkDebug("::MkFills - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;

    FSizeOfArray = FSize = size;

    if (FSizeOfArray == 0) {
       if (FFill!=NULL) delete[] (MkFill*)FFill;
       FFill = NULL;
       return;
    }

    if (FFill!=NULL) delete[] (MkFill*)FFill;
    FFill = new MkFill[FSizeOfArray];
    for (i=0;i<FSize;i++) FFill[i].Number = i;
}

void MkFills::Initialize(int size,MkFill *fills)
{
    int i;
    if (size < 0 || fills == NULL) {
      MkDebug("::MkFills - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0) {
       if (FFill!=NULL) delete[] (MkFill*)FFill;
       FFill = NULL;
       return;
    }

    if (FFill!=NULL) delete[] (MkFill*)FFill;
    FFill = new MkFill[FSizeOfArray];
    for (i=0;i<FSizeOfArray;i++) FFill[i] = fills[i];
    for (i=0;i<FSize;i++) FFill[i].Number = i;
}

int MkFills::Grow(int delta)
{
    int i;
    MkFill *fill=NULL;

    if (!(fill = new MkFill[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        fill[i] = FFill[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        fill[i] = NullFill;
    if (FFill) {
       delete[] (MkFill*)FFill;
       FFill = NULL;
    }
    FFill = fill;
    FSizeOfArray = FSizeOfArray+delta;
    return FSizeOfArray;
}

int MkFills::Shrink(int delta)
{
    int i;
    MkFill *fill=NULL;

    if (!(fill = new MkFill[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        fill[i] = FFill[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        fill[i] = NullFill;
    if (FFill) {
       delete[] (MkFill*)FFill;
       FFill = NULL;
    }
    FFill = fill;
    FSizeOfArray = FSizeOfArray-delta;
    return FSizeOfArray;
}

bool MkFills::Add(MkFill &fill)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }

    FSize++;
    FFill[FSize-1] = fill;
    return true;
}

bool MkFills::Add(int index, MkFill &fill)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }

    for (int i=FSize-1;i>=index;i--)
      FFill[i+1] = FFill[i];
    FSize++;
    FFill[index] = fill;
    return true;
}

bool MkFills::Delete(MkFill &fill)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FFill[i] == fill) break;
    }
    if(i==FSize) return false;
    if(FFill[i] == fill) {
      for (int j=i;j<FSize-1;j++)
        FFill[j] = FFill[j+1];
    }
    FSize--;
    FFill[FSize] = NullFill;
    return true;
}

bool MkFills::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FFill[j] = FFill[j+1];

    FSize--;
    FFill[FSize] = NullFill;
    return true;
}

bool MkFills::Clear()
{
   FSizeOfArray = FSize = 0;
   if (FFill) {
      delete[] FFill;
      FFill = NULL;
   }
   return true;
}

void MkFills::Out(char *fname)
{
  FILE *fp;
  fp = fopen(fname,"a");
  if(!fp) {
    return;
  }

  if(FSize==0) return;
    
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
  fprintf(fp,"\n");
  fprintf(fp,"                           <Information of Fill Stage>\n");
  fprintf(fp,"\n");
  fprintf(fp," Fill     Dist       G.L.\n");
  fprintf(fp,"  No.     (m)        (m) \n");
  fclose(fp);

  for(int i=0;i<FSize;i++)
    FFill[i].Out(fname);
}

MkFill & MkFills::operator[](int i)
{
    if (0<=i && i<FSize) return FFill[i];
    else if(FSize<=i && i<FSizeOfArray) {
      FSize=i+1;
      return FFill[i];
    }
    else if (FSizeOfArray <= i) {
      int size = FSizeOfArray;
      if(Grow(i-FSizeOfArray+5)==size) return NullFill;
      return FFill[i];
    }
    else return NullFill;
}

MkFills & MkFills::operator=(MkFills &fills)
{
    int i;

    Clear();
    FSize = fills.FSize;
    FSizeOfArray = fills.FSizeOfArray;
    if (FSizeOfArray == 0) {
       FFill = NULL;
       return *this;
    }
    this->FFill = new MkFill[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FFill[i] = fills.FFill[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FFill[i] = NullFill;

    return *this;
}

bool MkFills::operator==(MkFills &fills)
{
    int i;

    if (FSize != fills.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FFill[i] != fills.FFill[i]) return false;

    return true;
}

#ifdef __BCPLUSPLUS__
void MkFills::Draw(TObject *Sender)
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
           FFill[i].Draw(pb);
       }
       pb->Canvas->Pen->Color = C;
    }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void  MkFills::Draw(MkPaint *pb)
{

}
#endif
//---------------------------------------------------------------------------
