#include "MkLiuGrid.hpp"

MkLiuGrid NullLiuGrid(0);

MkLiuGrid::MkLiuGrid()
{
  Initialize();
}

MkLiuGrid::~MkLiuGrid()
{
  Clear();
}

void MkLiuGrid::Initialize()
{
  StartX = StartY = StartZ = EndX = EndY = EndZ = 0;
  ParticleRef.Clear();
}

void MkLiuGrid::Clear(void)
{
  StartX = StartY = StartZ = EndX = EndY = EndZ = 0;
  ParticleRef.Clear();
}

void MkLiuGrid::Setup(float sx, float sy, float sz, float ex, float ey, float ez)
{
  StartX = sx;
  StartY = sy;
  StartZ = sz;
  EndX = ex;
  EndY = ey;
  EndZ = ez;
}

void MkLiuGrid::Setup(int num)
{
  ParticleRef.Initialize(num);
  for (int i = 0; i < ParticleRef.getSzX(); i++)
  {
    ParticleRef(i) = -1;
  }
}

int MkLiuGrid::CountNParticles()
{
  int i, cnt = 0;
  for (i = 0; i < ParticleRef.getSzX(); i++)
  {
    if (ParticleRef(i) != -1)
      cnt++;
  }
  NParticles = cnt;
  return NParticles;
}

bool MkLiuGrid::Check(int parnum)
{
  bool flag = false;
  int i;
  for (i = 0; i < ParticleRef.getSzX(); i++)
  {
    if (ParticleRef(i) == parnum)
      flag = true;
  }
  return flag;
}

bool MkLiuGrid::Register(int parnum)
{
  bool flag = false;
  bool isFull = true;
  int i, I;

  if (ParticleRef.getSzX() <= 0)
  {
    Setup(1);
    //MkDebug("MkLiuGrid::Register ParticleRef setup");
  }

  for (i = 0; i < ParticleRef.getSzX(); i++)
  {
    if (ParticleRef(i) == parnum)
    {
      flag = true;
      I = i;
    }
  }
  if (!flag)
  {
    for (i = 0; i < ParticleRef.getSzX(); i++)
    {
      if (ParticleRef(i) == -1)
      {
        ParticleRef(i) = parnum;
        isFull = false;
        break;
      }
    }
    if (isFull)
    {
      MkInt parref(ParticleRef.getSzX() + 1);
      for (i = 0; i < ParticleRef.getSzX(); i++)
        parref(i) = ParticleRef(i);
      parref(i) = parnum;
      ParticleRef.Initialize(parref.getSzX());
      for (i = 0; i < ParticleRef.getSzX(); i++)
        ParticleRef(i) = parref(i);
    }
  }
  //MkDebug("Register %d\n",parnum);
  return flag;
}

bool MkLiuGrid::Unregister(int parnum)
{
  bool flag = false;
  int i, I;
  for (i = 0; i < ParticleRef.getSzX(); i++)
  {
    if (ParticleRef(i) == parnum)
    {
      flag = true;
      I = i;
      ParticleRef(i) = -1;
    }
  }
  //MkDebug("Unregister %d\n",parnum);
  return flag;
}

int MkLiuGrid::Out()
{
  for (int i = 0; i < ParticleRef.getSzX(); i++)
  {
    if (ParticleRef(i) != -1)
      printf("%3d, ", ParticleRef(i));
  }
  printf("\n");
  return -1;
}

//-----------------------------------------------------------
MkLiuGrids::MkLiuGrids()
{
  FLiuGrid.reset();
  FSize = NX = NY = NZ = 0;
  FSmoothLen = 0;
  XMin = YMin = ZMin = XMax = YMax = ZMax = 0;
}

MkLiuGrids::MkLiuGrids(int size, MkLiuGrid *rl)
{
  if (size < 0)
  {
    MkDebug("::MkLiuGrids - MkLiuGrids(int size)");
    throw Size(std::string("MkLiuGrids::Constructor negative size"), size);
  }
  if (rl == NULL)
  {
    MkDebug("::MkLiuGrids - MkLiuGrids(int size)");
    throw Size(std::string("MkLiuGrids::Constructor null parameter"), size);
  }

  FSize = size;
  if (FSize == 0)
  {
    FLiuGrid.reset();
    return;
  }

  try
  {
    FLiuGrid.reset(new MkLiuGrid[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuGrids::Setup memory allocation error\n");
    throw Alloc(a.what());
  }

  for (int i = 0; i < FSize; i++)
    FLiuGrid[i] = rl[i];
}

MkLiuGrids::MkLiuGrids(int size)
{
  if (size < 0)
  {
    MkDebug("::MkLiuGrids - MkLiuGrids(int size)");
    throw Size(std::string("MkLiuGrids::Constructor negative size"), size);
  }

  FSize = size;
  if (FSize == 0)
  {
    FLiuGrid.reset();
    return;
  }

  try
  {
    FLiuGrid.reset(new MkLiuGrid[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuGrids::Setup memory allocation error\n");
    throw Alloc(a.what());
  }
}

MkLiuGrids::~MkLiuGrids()
{
  FSize = 0;
  FLiuGrid.reset();
}

void MkLiuGrids::Initialize(int size)
{
  if (size < 0)
  {
    MkDebug("::MkLiuGrids - Initialize(int size)");
    throw Size(std::string("MkLiuGrids::Initialize negative size"), size);
  }
  if (FSize == size)
    return;

  FSize = size;
  if (FSize == 0)
  {
    FLiuGrid.reset();
    return;
  }

  try
  {
    FLiuGrid.reset(new MkLiuGrid[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuPairs::Constructor mem allocation error");
    throw Alloc(a.what());
  }
}

void MkLiuGrids::Initialize(int size, MkLiuGrid *rl)
{
  if (size < 0)
  {
    MkDebug("::MkLiuGrids - Initialize(int size)");
    throw Size(std::string("MkLiuGrids::Initialize negative size"), size);
  }
  if (rl == NULL)
  {
    MkDebug("::MkLiuGrids - Initialize(int size)");
    throw Size(std::string("MkLiuGrids::Initialize null parameter"), size);
  }

  FSize = size;
  if (FSize == 0)
  {
    FLiuGrid.reset();
    return;
  }

  try
  {
    FLiuGrid.reset(new MkLiuGrid[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuPairs::Initialize mem allocation error");
    throw Alloc(a.what());
  }
  for (int i = 0; i < FSize; i++)
    FLiuGrid[i] = rl[i];
}

void MkLiuGrids::Grow(int size) // Grow allocate extra memory
{
  if (size <= 0)
  {
    MkDebug("::MkLiuGrids - Grow(int size)");
    throw Size(std::string("MkLiuGrids::Grow negative size"), size);
  }

  boost::shared_array<MkLiuGrid> pLine;

  try
  {
    pLine.reset(new MkLiuGrid[FSize + size]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuPairs::Grow mem allocation error");
    throw Alloc(a.what());
  }

  for (int i = 0; i < FSize; i++)
    pLine[i] = FLiuGrid[i];

  FLiuGrid = pLine;
  FSize = FSize + size;
}

void MkLiuGrids::Clear()
{
  FSize = 0;
  FLiuGrid.reset();
}

void MkLiuGrids::Setup()
{
  double sx, sy, sz, ex, ey, ez;
  GridLen = FSmoothLen * 2;
  Clear();

  if ((XMax - XMin) < 0.001 && (YMax - YMin) < 0.001 && (ZMax - ZMin) < 0.001)
  {
    MkDebug("MkLiuGrids::Setup boundary strange\n");
    exit(-1);
  }
  NX = int((XMax - XMin) / GridLen) + 3;
  NY = int((YMax - YMin) / GridLen) + 3;
  NZ = int((ZMax - ZMin) / GridLen) + 3;
  FSize = NX * NY * NZ;
  MkDebug("GridLen %f NX %d,NY %d,NZ %d,FSize %d\n", GridLen, NX, NY, NZ, FSize);
  MkDebug("XMin %f XMax %f YMin %f YMax %f ZMin %f ZMax %f\n", XMin, XMax, YMin, YMax, ZMin, ZMax);
  try
  {
    FLiuGrid.reset(new MkLiuGrid[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuGrids::Setup memory allocation error\n");
    throw Alloc(a.what());
  }

  for (int i = 0; i < NX; i++)
  {
    sx = XMin + GridLen * (i - 1);
    ex = XMin + GridLen * (i + 0);
    for (int j = 0; j < NY; j++)
    {
      sy = YMin + GridLen * (j - 1);
      ey = YMin + GridLen * (j + 0);
      for (int k = 0; k < NZ; k++)
      {
        sz = ZMin + GridLen * (k - 1);
        ez = ZMin + GridLen * (k + 0);
        (*this)(i, j, k).Setup(sx, sy, sz, ex, ey, ez);
        //(*this)(i,j,k).Setup(0); // need to be fixed...
      }
    }
  }
}

void MkLiuGrids::Setup(MkLiuParticles &par)
{
  MkDebug("MkLiuGrids::Setup( with par) \n");
  Setup();
  MkDebug("MLG::Setup 01\n");
  Update(par);
  MkDebug("~MkLiuGrids::Setup (with par) \n");
}

void MkLiuGrids::Update(MkLiuParticles &par)
{
  int i, j, k, NP = par.GetSize();
  int I, J, K;
  int GI, GJ, GK;
  double X, Y, Z, SX, SY, SZ, EX, EY, EZ;

  for (i = 0; i < NP; i++)
  {
    MkLiuParticle &par_i = par[i];
    X = par_i.X;
    Y = par_i.Y;
    Z = par_i.Z;
    I = int((X - XMin) / GridLen) + 1;
    J = int((Y - YMin) / GridLen) + 1;
    K = int((Z - ZMin) / GridLen) + 1;
    GI = par_i.GI;
    GJ = par_i.GJ;
    GK = par_i.GK;

    if ((GI == I) && (GJ == J) && (GK == K))
      continue;

    FLiuGrid[GI + GJ * NX + GK * NX * NY].Unregister(i);
    FLiuGrid[I + J * NX + K * NX * NY].Register(i);
    par_i.SetGrid(I, J, K);
  }
}

void MkLiuGrids::Update_backup(MkLiuParticles &par)
{
  int i, j, k, NP = par.GetSize();
  int I, J, K;
  int GI, GJ, GK;
  double X, Y, Z, SX, SY, SZ, EX, EY, EZ;

  for (i = 0; i < NP; i++)
  {
    for (I = 0; I < NX; I++)
    {
      for (J = 0; J < NY; J++)
      {
        for (K = 0; K < NZ; K++)
        {
          MkLiuGrid &grid = FLiuGrid[I + J * NX + K * NX * NY]; //(*this)(I,J,K);
          MkLiuParticle &par_i = par[i];
          X = par_i.X;
          Y = par_i.Y;
          Z = par_i.Z;
          SX = grid.StartX;
          SY = grid.StartY;
          SZ = grid.StartZ;
          EX = grid.EndX;
          EY = grid.EndY;
          EZ = grid.EndZ;
          //MkDebug("SX %f,X %f,EX %f,SY %f,Y %f,EY %f,SZ %f,Z %f,EZ %f \n",SX,X,EX,SY,Y,EY,SZ,Z,EZ);
          if ((SX < X && X < EX) && (SY < Y && Y < EY) && (SZ < Z && Z < EZ))
          {
            GI = par_i.GI;
            GJ = par_i.GJ;
            GK = par_i.GK;

            if ((I != GI) || (J != GJ) || (K != GK))
            {
              (*this)(GI, GJ, GK).Unregister(i);
              grid.Register(i);
              par_i.SetGrid(I, J, K);
              //MkDebug("In the grid cell (I %d,J %d ,K %d),particle i registered %d \n",I,J,K, i);
            }
          }
          //grid.CountNParticles();
        }
      }
    }
  }
  //MkDebug("NX %d,NY %d,NZ %d\n",NX,NY,NZ);
  /*
  for (i=0;i<FSize;i++) {

    //    i = I+J*NX+K*NX*NY;
    float sx,sy,ex,ey;
    K = i/(NX*NY);
    J = (i - K*NX*NY)/NX;
    I = i - J*NX - K*NX*NY;
    sx = FLiuGrid[i].StartX;
    sy = FLiuGrid[i].StartY;
    ex = FLiuGrid[i].EndX;
    ey = FLiuGrid[i].EndY;

    
    MkDebug("grid no %3d, (I %3d J %3d K %3d) Par nos %3d sx %f sy %f ex %f ey %f ",i,I,J,K,FLiuGrid[i].NParticles,sx,sy,ex,ey);
    for (j=0;j<FLiuGrid[i].ParticleRef.getSzX();j++) {
      MkDebug(" %3d ",FLiuGrid[i].ParticleRef(j));
    } 
    MkDebug("\n");
  }
  */
  //usleep(999999);usleep(999999);usleep(999999);
}

int MkLiuGrids::GetNumOfParticles(int I, int J, int K)
{
  int II, JJ, KK;
  int np = 0;
  int np_ = 0;

  //MkDebug("GNP2 np %d\n",np);

  II = I - 1;
  JJ = J - 1;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J - 1;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J - 1;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP3 np %d\n",np);

  II = I - 1;
  JJ = J;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP4 np %d\n",np);

  II = I - 1;
  JJ = J + 1;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J + 1;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J + 1;
  KK = K - 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP5 np %d\n",np);

  II = I - 1;
  JJ = J - 1;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J - 1;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J - 1;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP6 np %d\n",np);

  II = I - 1;
  JJ = J;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  np_ += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP7 np %d\n",np);

  II = I - 1;
  JJ = J + 1;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J + 1;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J + 1;
  KK = K;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP8 np %d\n",np);

  II = I - 1;
  JJ = J - 1;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J - 1;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J - 1;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP9 np %d\n",np);

  II = I - 1;
  JJ = J;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("GNP10 np %d\n",np);

  II = I - 1;
  JJ = J + 1;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I;
  JJ = J + 1;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();
  II = I + 1;
  JJ = J + 1;
  KK = K + 1;
  np += (*this)(II, JJ, KK).GetNumOfParticle();

  //MkDebug("Get (I %d, J %d, K %d) Num of Particles = %d, No of Neighbor %d \n",I,J,K,np_,np);

  return np;
}

MkInt &MkLiuGrids::GetParticles(int I, int J, int K)
{
  MkInt *par;
  //MkDebug("MkLiuGrids::GetParticles()\n");
  int NP = GetNumOfParticles(I, J, K);
  //MkDebug("GP1 ");
  int np = 0, i, II, JJ, KK;

  static MkInt ParRef;
  ParRef.Clear();

  if (0 == NP)
    return ParRef;

  ParRef.Initialize(NP);

  //MkDebug("NP = %d, size of ParRef %d \n",NP,ParRef.getSzX());

  //MkDebug("GP2 np %d\n",np);

  II = I - 1;
  JJ = J - 1;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J - 1;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J - 1;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP3 np %d\n",np);

  II = I - 1;
  JJ = J;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP4 np %d\n",np);

  II = I - 1;
  JJ = J + 1;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J + 1;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J + 1;
  KK = K - 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP5 np %d\n",np);

  II = I - 1;
  JJ = J - 1;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J - 1;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J - 1;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP6 np %d\n",np);

  II = I - 1;
  JJ = J;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef(); //MkDebug("II %d JJ %d KK %d \n",II,JJ,KK);
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef(); //MkDebug("II %d JJ %d KK %d \n",II,JJ,KK);
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef(); //MkDebug("II %d JJ %d KK %d \n",II,JJ,KK);
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP7 np %d\n",np);

  II = I - 1;
  JJ = J + 1;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef(); //MkDebug("II %d JJ %d KK %d \n",II,JJ,KK);
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J + 1;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef(); //MkDebug("II %d JJ %d KK %d \n",II,JJ,KK);
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J + 1;
  KK = K;
  par = &(*this)(II, JJ, KK).GetParticleRef(); //MkDebug("II %d JJ %d KK %d \n",II,JJ,KK);
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP8 np %d\n",np);

  II = I - 1;
  JJ = J - 1;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J - 1;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J - 1;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP9 np %d\n",np);

  II = I - 1;
  JJ = J;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }

  //MkDebug("GP10 np %d\n",np);

  II = I - 1;
  JJ = J + 1;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I;
  JJ = J + 1;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  II = I + 1;
  JJ = J + 1;
  KK = K + 1;
  par = &(*this)(II, JJ, KK).GetParticleRef();
  for (i = 0; i < par->getSzX(); i++)
    if ((*par)[i] != -1)
    {
      ParRef[np] = (*par)[i];
      np++;
    }
  /*
  for (i=0;i<ParRef.getSzX();i++) {
    MkDebug(" ParRef[%d] =  %d,\n",i,ParRef[i]);
  }
  MkDebug("\n");
  */
  return ParRef;
}

MkLiuGrid &MkLiuGrids::operator()(int i, int j, int k)
{
  if ((0 <= i && i < NX) && (0 <= j && j < NY) && (0 <= k && k < NZ))
  {
    return FLiuGrid[i + j * NX + k * NX * NY];
  }
  else
    return NullLiuGrid;
}

MkLiuGrid &MkLiuGrids::operator[](int i)
{
  if (FSize == 0)
    return NullLiuGrid;
  else if (i >= 0 && i < FSize)
    return FLiuGrid[i];
  else
    return NullLiuGrid;
}

MkLiuGrids &MkLiuGrids::operator=(MkLiuGrids &rls)
{
  int i;

  Clear();
  FSize = rls.FSize;
  if (FSize == 0)
  {
    FLiuGrid.reset();
    return *this;
  }

  try
  {
    FLiuGrid.reset(new MkLiuGrid[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuGrids::Setup memory allocation error\n");
    throw Alloc(a.what());
  }

  for (i = 0; i < FSize; i++)
    FLiuGrid[i] = rls.FLiuGrid[i];

  return *this;
}
//---------------------------------------------------------------------------
