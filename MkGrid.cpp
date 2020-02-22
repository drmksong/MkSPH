//---------------------------------------------------------------------------
#include "MkGrid.h"
//---------------------------------------------------------------------------
MkSPHGrid NullSPHGrid;
MkSPHGrids NullSPHGrids(0);

MkTopoGrid::MkTopoGrid()
{
    XMin = 0;
    XMax = 0;
    YMin = 0;
    YMax = 0;
    ElevMin = 0;
    ElevMax = 0;
    Assigned = false;
}

MkTopoGrid::MkTopoGrid(MkFloat &x,MkFloat &y,MkFloat &elev)
{
    X.CopyFrom(x);
    NX = X.getSzX();
    Y.CopyFrom(y);
    NY = Y.getSzX();
    Elev.CopyFrom(elev);
    XMin = X(0);XMax = X(X.getSzX());
    YMin = Y(0);YMax = Y(Y.getSzX());
    ElevMin = Elev(0);ElevMax = Elev(Elev.getSzX());
    Assigned = true;
}

void MkTopoGrid::SetTopoGrid(MkFloat &x,MkFloat &y,MkFloat &elev)
{
    X.CopyFrom(x);
    NX = X.getSzX();
    Y.CopyFrom(y);
    NY = Y.getSzX();
    Elev.CopyFrom(elev);
    XMin = X(0);XMax = X(X.getSzX()-1);
    YMin = Y(0);YMax = Y(Y.getSzX()-1);
    ElevMin = Elev(0);ElevMax = Elev(Elev.getSzX()-1);
    Assigned = true;
}

float MkTopoGrid::operator()(float x,float y)
{
  MkPoint rp[4];
  MkTriangle rt[2];
  MkLine rl[2];
  float elev[4];
  int I,J;
  int i,j;

  if (x < XMin || x > XMax || y < YMin || y > YMax) return 0;

  for (i = 0 ; i < NX-1 ; i++) {
    if (x >= X(i) && x <= X(i+1)) {
       I = i;
       break;
    }
  }
  for (j = 0; j < NY-1 ; j++ ) {
    if (y >= Y(j) && y <= Y(j+1)) {
       J = j;
       break;
    }
  }
  if (i == NX-1 || j == NY-1) return 0;

  rp[0].SetPoint(X(I  ),Y(J  ),Elev(I,J));
  rp[1].SetPoint(X(I+1),Y(J  ),Elev((I+1),J));
  rp[2].SetPoint(X(I+1),Y(J+1),Elev((I+1),J+1));
  rp[3].SetPoint(X(I  ),Y(J+1),Elev(I,J+1));

  rt[0].Reset(rp[0],rp[1],rp[3]);
  rt[1].Reset(rp[1],rp[2],rp[3]);

  if (rt[0].isIn(x,y)) {
     return rt[0](x,y);
  }

  else if (rt[1].isIn(x,y)) {
     return rt[1](x,y);
  }
  else return 0;
}

MkTopoGrid & MkTopoGrid::operator=(MkTopoGrid &tg)
{
     X.CopyFrom(tg.X);
     Y.CopyFrom(tg.Y);
     Elev.CopyFrom(tg.Elev);
     XMin     =  tg.XMin;
     XMax     =  tg.XMax;
     YMin     =  tg.YMin;
     YMax     =  tg.YMax;
     ElevMin  =  tg.ElevMin;
     ElevMax  =  tg.ElevMax;
     NX       =  tg.NX;
     NY       =  tg.NY;
     Assigned =  true;
}
//-----------------------------------------------------------
MkSPHGrid::MkSPHGrid()
{
  StartX = StartY = StartZ = EndX = EndY = EndZ = 0;
  ParticleRef.Clear();
}

MkSPHGrid::~MkSPHGrid()
{
  ParticleRef.Clear();
}

void MkSPHGrid::Clear(void)
{
  StartX = StartY = StartZ = EndX = EndY = EndZ = 0;
  ParticleRef.Clear();
}

void MkSPHGrid::Setup(float sx, float sy, float sz, float ex, float ey, float ez)
{
  StartX = sx;
  StartY = sy;
  StartZ = sz;
  EndX = ex;
  EndY = ey;
  EndZ = ez;
}

void MkSPHGrid::Setup(int num)
{
  ParticleRef.Initialize(num);
  for (int i=0;i<ParticleRef.getSzX();i++) {
    ParticleRef[i]=-1;
  }
}

bool MkSPHGrid::Check(int parnum)
{
  bool flag=false; int i, I;
  for (i=0;i<ParticleRef.getSzX();i++) {
    if (ParticleRef[i] == parnum) {
      flag = true;
      I = i;
    }
  }
  return flag;
}

bool MkSPHGrid::Register(int parnum)
{
  bool flag=false;
  bool isFull=true; 
  int i, I;
  for (i=0;i<ParticleRef.getSzX();i++) {
    if (ParticleRef[i] == parnum) {
      flag = true;
      I = i;
    }
  }
  if (!flag) {
    for (i=0;i<ParticleRef.getSzX();i++) {
      if(ParticleRef[i] == -1) {
	ParticleRef[i] = parnum;
	isFull = false;
	break;
      }
    }
    if(isFull) {
      MkInt parref(ParticleRef.getSzX()+1);
      for (i=0;i<ParticleRef.getSzX();i++) parref[i]= ParticleRef[i];
      parref[i] = parnum;
      ParticleRef.Initialize(parref.getSzX());
      for (i=0;i<ParticleRef.getSzX();i++) ParticleRef[i] = parref[i];
    }
  }  

  return flag;
}

bool MkSPHGrid::Unregister(int parnum)
{
  bool flag=false;
  int i, I;
  for (i=0;i<ParticleRef.getSzX();i++) {
    if (ParticleRef[i] == parnum) {
      flag = true;
      I = i;
      ParticleRef[i] = -1;
    }
  }
  
  return flag;
}

int MkSPHGrid::Out()
{
  for (int i=0;i<ParticleRef.getSzX();i++) {
    if (ParticleRef[i]!=-1) printf("%3d, ",ParticleRef[i]);
  }
  printf("\n");
  return -1;
}


//-----------------------------------------------------------

 MkSPHGrids::MkSPHGrids(int size,MkSPHGrid *rl)
{
   if (size < 0) {
#ifdef __BCPLUSPLUS__
     ShowMessage("::MkSPHGrids - MkSPHGrids(int size)");
#else 
     MkDebug("::MkSPHGrids - MkSPHGrids(int size)");
#endif
     return;
   }

   FSize = size;
   if (FSize == 0) {
      FSPHGrid = NULL;
      return;
   }

   FSPHGrid = new MkSPHGrid[FSize];
   for (int i=0;i<FSize;i++) FSPHGrid[i] = rl[i];
}

 MkSPHGrids::MkSPHGrids(int size)
{
   if (size < 0) {
#ifdef __BCPLUSPLUS__
     ShowMessage("::MkSPHGrids - MkSPHGrids(int size)");
#else 
     MkDebug("::MkSPHGrids - MkSPHGrids(int size)");
#endif

     return;
   }

   FSize = size;
   if (FSize == 0) {
      FSPHGrid = NULL;
      return;
   }

   FSPHGrid = new MkSPHGrid[FSize];
}

 MkSPHGrids::~MkSPHGrids()
{
   FSize = 0;
   if (FSPHGrid) {
      delete[] (MkSPHGrid*)FSPHGrid;
      FSPHGrid = NULL;
   }
}

void MkSPHGrids::Initialize(int size)
{
    if (size < 0) {
#ifdef __BCPLUSPLUS__
     ShowMessage("::MkSPHGrids - Initialize(int size)");
#else 
     MkDebug("::MkSPHGrids - Initialize(int size)");
#endif

      return;
    }
    if (FSize == size) return;

    FSize = size;
    if (FSize == 0) {
       if (FSPHGrid!=NULL) delete[] (MkSPHGrid*)FSPHGrid;
       FSPHGrid = NULL;
       return;
    }

    if (FSPHGrid!=NULL) delete[] (MkSPHGrid*)FSPHGrid;
    FSPHGrid = new MkSPHGrid[FSize];
}

void MkSPHGrids::Initialize(int size,MkSPHGrid *rl)
{
    if (size < 0) {
#ifdef __BCPLUSPLUS__
     ShowMessage("::MkSPHGrids - Initialize(int size)");
#else 
     MkDebug("::MkSPHGrids - Initialize(int size)");
#endif

      return;
    }

    FSize = size;
    if (FSize == 0) {
       if (FSPHGrid!=NULL) delete[] (MkSPHGrid*)FSPHGrid;
       FSPHGrid = NULL;
       return;
    }

    if (FSPHGrid!=NULL) delete[] (MkSPHGrid*)FSPHGrid;
    FSPHGrid = new MkSPHGrid[FSize];
    for (int i=0;i<FSize;i++) FSPHGrid[i] = rl[i];
}

void MkSPHGrids::Grow(int size) // Grow allocate extra memory
{
    if (size <=0) {
      return;
    }

    MkSPHGrid *pLine = new MkSPHGrid[FSize+size];

    for (int i = 0 ; i < FSize ; i++)
        pLine[i] = FSPHGrid[i];

    if (FSPHGrid!=NULL) {
       delete[] (MkSPHGrid*)FSPHGrid;
       FSPHGrid = NULL;
    }

    FSPHGrid = pLine;
    FSize = FSize + size;
}

bool MkSPHGrids::Clear()
{
   FSize = 0;
   if (FSPHGrid) {
      delete[] (MkSPHGrid*)FSPHGrid;
      FSPHGrid = NULL;
   }
   return true;
}

MkSPHGrid & MkSPHGrids::operator[](int i)
{
    if (FSize == 0) return NullSPHGrid;
    else if (i >=0 && i < FSize) return FSPHGrid[i];
    else return NullSPHGrid;
}


MkSPHGrids & MkSPHGrids::operator=(MkSPHGrids &rls)
{
    int i;

    Clear();
    FSize = rls.FSize;
    if (FSize == 0) {
       FSPHGrid = NULL;
       return *this;
    }
    FSPHGrid = new MkSPHGrid[FSize];

    for (i=0;i<FSize;i++)
      FSPHGrid[i] = rls.FSPHGrid[i];

    return *this;
}
//---------------------------------------------------------------------------
