//234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
//------------------------------------------------------------------------------
#include "MkLiuParticle.hpp"

MkLiuParticle NullParticle(0, 0, 0);
MkLiuParticles NullParticles(0);
//double Zero = 0;

void MkLiuParticle::Initialize()
{
  ColorType = ctNone;
  R = G = B = 0;
  Radius = 0;
  Rho = 0; //??
  Mass = 0;
  Volume = 0;
  Eta = 0;
  Temp = 0; //??
  Press = 0;
  Energy = 0;
  SoundSpeed = 0;
  SmoothLen = 0;
  CountIac = 0;
  GI = GJ = GK = -1;

  X = 0;
  Y = 0;
  Z = 0;

  XVel = 0;
  YVel = 0;
  ZVel = 0;

  XAVel = 0;
  YAVel = 0;
  ZAVel = 0;

  ParticleType = 0; // IType Particle Type
  Virtual = false;

  DUDt = 0; //**
  INDUDt = 0;
  AVDUDt = 0;
  AHDUDt = 0;
  ADDUDt = 0; // artificial drag internal energy
  APDUDt = 0; // artificial repel internal energy
  DRhoDt = 0;
  DVXDt = 0; //**
  INDVXDt = 0;
  ARDVXDt = 0;
  ADDVXDt = 0;
  APDVXDt = 0;
  EXDVXDt = 0;
  DVYDt = 0; //**
  INDVYDt = 0;
  ARDVYDt = 0;
  ADDVYDt = 0;
  APDVYDt = 0;
  EXDVYDt = 0;
  DVZDt = 0; //**
  INDVZDt = 0;
  ARDVZDt = 0;
  ADDVZDt = 0;
  APDVZDt = 0;
  EXDVZDt = 0;
  TDSDt = 0;
}

void MkLiuParticle::SetupColor()
{
  double rho;
  R = G = B = 0;
  //rho = Rho - 600;
  switch (ColorType)
  {
  case ctRho:
    /*    R = (Rho>1000) ? (Rho-1000)/1000.0 : 0;
    B = (Rho<1000) ? 1-Rho/1000.0 : 0;
    G = (Rho<1000) ? Rho/1000.0:1-(Rho-1000.0)/1000.0;
    //   MkDebug("X : %f, Y : %f, Rho : %f, Mass : %f, R : %f, G : %f, B : %f \n",X,Y,Rho,Mass,R,G,B);
    break;
    R = (rho<200) ? rho/200.0 : (rho-200)/200.0;
    B = (200<rho&&rho<400) ? (rho-200)/200.0 : 0;
    G = (rho>200) ? (rho-200)/100.0:1-(Rho-200.0)/100.0;
*/
    if (Rho < 800)
    {
      B = (800 - Rho) / 400.0;
      G = (Rho - 400) / 400.0;
      R = 0;
    }
    else if (800 < Rho && Rho < 1200)
    {
      B = 0;
      G = 1;
      R = (Rho - 800) / 400.0;
    }
    else if (1200 < Rho && Rho < 1600)
    {
      B = 0;
      G = (1600 - Rho) / 400.0;
      R = 1;
    }
    else if (400 > Rho)
    {
      B = 1;
      G = 0;
      R = 0;
    }
    else if (1600 < Rho)
    {
      B = 0;
      G = 0;
      R = 1;
    }
    else
    {
      R = G = B = 1;
    }

    break;

  case ctPress:
    R = (Press > 50) ? (Press - 50) / 50.0 : 0;
    B = (Press < 50) ? 1 - Press / 50.0 : 0;
    G = (Press < 50) ? Press / 50.0 : 1 - (Press - 50.0) / 50.0;
    break;

  default:
    R = 1;
    G = 0;
    B = 0;
  }
}

#if defined(__GL_H__)
void MkLiuParticle::Draw()
{
  SetupColor();
  glPushMatrix();
  //    glLoadIdentity();
  //    glScalef(2,2,2);
  //    glTranslatef(1,1,0);
  glColor3f(R, G, B);
  glTranslatef(-X, -Y, -Z);
  glutSolidSphere(SmoothLen / 2.0, 16, 16);
  glPopMatrix();
}
#endif

MkLiuParticle &MkLiuParticle::operator=(const MkLiuParticle &lp)
{

  MkPoint::operator=((MkPoint &)lp);

  R = lp.R;
  G = lp.G;
  B = lp.B;
  Radius = lp.Radius;
  ColorType = lp.ColorType;

  XVel = lp.XVel;
  YVel = lp.YVel;
  ZVel = lp.ZVel;
  XAVel = lp.XAVel;
  YAVel = lp.YAVel;
  ZAVel = lp.ZAVel;
  ParticleType = lp.ParticleType;
  Virtual = lp.Virtual;

  Mass = lp.Mass;             // real mass checked
  Rho = lp.Rho;               // density checked,
  Rho_Norm = lp.Rho_Norm;     // norm density for Norm_Density()
  Volume = lp.Volume;         // volume to calculate the selfdens (density when no neighbor particles)  mk
  Eta = lp.Eta;               // dynamic viscosity
  Press = lp.Press;           // P  checked
  Temp = lp.Temp;             // T  checked
  Energy = lp.Energy;         //U  checked
  SoundSpeed = lp.SoundSpeed; //C
  SmoothLen = lp.SmoothLen;   //Hsml checked

  CountIac = lp.CountIac;
  GI = lp.GI;
  GJ = lp.GJ;
  GK = lp.GK;

  DUDt = lp.DUDt; //**
  INDUDt = lp.INDUDt;
  AVDUDt = lp.AVDUDt;
  AHDUDt = lp.AHDUDt;
  ADDUDt = lp.ADDUDt; // artificial drag internal energy
  APDUDt = lp.APDUDt; // artificial repel internal energy
  DRhoDt = lp.DRhoDt;
  DVXDt = lp.DVXDt; //**
  INDVXDt = lp.INDVXDt;
  ARDVXDt = lp.ARDVXDt;
  ADDVXDt = lp.ADDVXDt; // artificial drag
  APDVXDt = lp.APDVXDt; // artificial repel
  EXDVXDt = lp.EXDVXDt;
  DVYDt = lp.DVYDt; //**
  INDVYDt = lp.INDVYDt;
  ARDVYDt = lp.ARDVYDt;
  ADDVYDt = lp.ADDVYDt; // artificial drag
  APDVYDt = lp.APDVYDt; // artificail repel
  EXDVYDt = lp.EXDVYDt;
  DVZDt = lp.DVZDt; //**
  INDVZDt = lp.INDVZDt;
  ARDVZDt = lp.ARDVZDt;
  ADDVZDt = lp.ADDVZDt; // artificial drag
  APDVZDt = lp.APDVZDt; // artificial repel
  EXDVZDt = lp.EXDVZDt;
  TDSDt = lp.TDSDt;
  /*
   W=lp.W;  // kernal of all interaction pairs
   DWDX=lp.DWDX;//derivative of kernal with respect to x
   DWDY=lp.DWDY;//derivative of kernal with respect to x
   DWDZ=lp.DWDZ;//derivative of kernal with respect to x
   */
  return (*this);
}

MkLiuParticle &operator*(MkLiuParticle &rp, MkMatrix4<double> &rm)
{
  static MkLiuParticle rp_t;
  rp_t = rp;

  rp_t.X = rp.X * rm(0, 0) + rp.Y * rm(0, 1) + rp.Z * rm(0, 2) + 1 * rm(0, 3);
  rp_t.Y = rp.X * rm(1, 0) + rp.Y * rm(1, 1) + rp.Z * rm(1, 2) + 1 * rm(1, 3);
  rp_t.Z = rp.X * rm(2, 0) + rp.Y * rm(2, 1) + rp.Z * rm(2, 2) + 1 * rm(2, 3);

  return rp_t;
}

MkLiuParticle &operator*(MkLiuParticle &rp, double f)
{
  static MkLiuParticle rp_t;
  rp_t = rp;

  rp_t.X *= f;
  rp_t.Y *= f;
  rp_t.Z *= f;
  return rp_t;
}

//------------------------------------------------------------------------------

MkLiuParticles::MkLiuParticles(int size, MkLiuParticle *rps)
{

  if (size < 0)
  {
    MkDebug("::MkLiuParticles - MkLiuParticles(int size) negative size");
    throw Size(std::string("MkLiuParticles Constructor negative size"), size);
  }

  if (rps == NULL)
  {
    MkDebug("::MkLiuParticles - MkLiuParticles(int size) NULL MkLiuParticle");
    throw Size(std::string("MkLiuParticles Constructor null pointer argument"), 0);
  }

  FSizeOfArray = FSize = size;
  if (FSize == 0)
  {
    FParticle.reset();
    return;
  }
  try
  {
    FParticle.reset(new MkLiuParticle[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuParticles::Constructor throws bad_alloc during allocating the memory");
    throw Alloc(a.what());
  }

  for (int i = 0; i < FSize; i++)
    (*this)[i] = rps[i];
}

MkLiuParticles::MkLiuParticles(int size)
{
  if (size < 0)
  {
    MkDebug("::MkLiuParticles - MkLiuParticles(int size)");
    throw Size(std::string("MkLiuParticles Constructor negative size"), size);
  }

  FSize = FSizeOfArray = size;

  if (FSizeOfArray == 0)
  {
    FParticle.reset();
    return;
  }

  try
  {
    FParticle.reset(new MkLiuParticle[FSizeOfArray]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuParticles::Constructor throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }
}

MkLiuParticles::~MkLiuParticles()
{
  FSizeOfArray = FSize = 0;
  if (FParticle)
  {
    FParticle.reset();
  }
}

void MkLiuParticles::Initialize(int size)
{
  if (size < 0)
  {
    MkDebug("::MkLiuParticles - Initialize(int size)");
    throw Size(std::string("MkLiuParticles Constructor negative size"), size);
  }
  if (FSizeOfArray == size)
    return;

  FSize = FSizeOfArray = size;

  if (FSizeOfArray == 0)
  {
    FParticle.reset();
    return;
  }

  try
  {
    FParticle.reset(new MkLiuParticle[FSizeOfArray]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuParticles::Initialize throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }
}

void MkLiuParticles::Initialize(int size, MkLiuParticle *rps)
{

  if (size < 0)
  {
    MkDebug("::MkLiuParticles - Initialize(int size)");
    throw Size(std::string("MkLiuParticles Constructor negative size"), size);
  }
  if (rps == NULL)
  {
    MkDebug("::MkLiuParticles - Initialize(int size)");
    throw Size(std::string("MkLiuParticles Constructor null pointer argument"), 0);
  }

  if (FSizeOfArray == size)
  {
    FSize = size;
    for (int i = 0; i < FSizeOfArray; i++)
      FParticle[i] = rps[i];
    return;
  }

  FSize = FSizeOfArray = size;
  if (FSizeOfArray == 0)
  {
    FParticle.reset();
    return;
  }

  try
  {
    FParticle.reset(new MkLiuParticle[FSizeOfArray]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuParticles::Initialize throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }

  for (int i = 0; i < FSizeOfArray; i++)
    FParticle[i] = rps[i];
}

void MkLiuParticles::FindCenter() // find maximum and minimum for each axis and average it, the center
{
  if (FSize == 0)
  {
    FCenter = NullParticle;
    return;
  }

  double x1, y1, z1, x2, y2, z2;
  x1 = (*this)[0].X;
  y1 = (*this)[0].Y;
  z1 = (*this)[0].Z;
  x2 = (*this)[0].X;
  y2 = (*this)[0].Y;
  z2 = (*this)[0].Z;

  for (int i = 1; i < GetSize(); i++)
  {
    x1 = x1 < (*this)[i].X ? x1 : (*this)[i].X;
    y1 = y1 < (*this)[i].Y ? y1 : (*this)[i].Y;
    z1 = z1 < (*this)[i].Z ? z1 : (*this)[i].Z;
    x2 = x2 > (*this)[i].X ? x2 : (*this)[i].X;
    y2 = y2 > (*this)[i].Y ? y2 : (*this)[i].Y;
    z2 = z2 > (*this)[i].Z ? z2 : (*this)[i].Z;
  }
  FCenter.X = (x1 + x2) / 2.0;
  FCenter.Y = (y1 + y2) / 2.0;
  FCenter.Z = (z1 + z2) / 2.0;
}

int MkLiuParticles::Grow(int delta)
{
  int i;
  boost::shared_array<MkLiuParticle> rp;

  try
  {
    rp.reset(new MkLiuParticle[FSizeOfArray + delta]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuParticles::Grow throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }

  for (i = 0; i < FSize; i++)
    rp[i] = FParticle[i];
  for (i = FSize; i < FSizeOfArray + delta; i++)
    rp[i] = NullParticle;

  FParticle = rp; // it does not leak memory as it is boost smart pointers, shared_array
  FSizeOfArray = FSizeOfArray + delta;
  return FSizeOfArray;
}

int MkLiuParticles::Shrink(int delta)
{
  int i;
  boost::shared_array<MkLiuParticle> rp;

  try
  {
    rp.reset(new MkLiuParticle[FSizeOfArray - delta]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuParticles::Shrink throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }

  for (i = 0; i < FSize; i++)
    rp[i] = FParticle[i];
  for (i = FSize; i < FSizeOfArray - delta; i++)
    rp[i] = NullParticle;

  FParticle = rp; // it does not leak memory as it is boost smart pointers, shared_array
  FSizeOfArray = FSizeOfArray - delta;
  return FSizeOfArray;
}

bool MkLiuParticles::Add(MkLiuParticle point)
{
  if (FSize == FSizeOfArray)
  {
    try
    {
      Grow(FSize - FSizeOfArray + 10);
    }
    catch (Alloc &a)
    {
      MkDebug("MkLiuParticles::Add throws Alloc while trying to Grow memory");
      throw Alloc(a.what());
    }
  }
  if (FSize == FSizeOfArray)
    return false;

  FSize++;
  FParticle[FSize - 1] = point;
  return true;
}

bool MkLiuParticles::Add(int index, MkLiuParticle point)
{
  if (FSize == FSizeOfArray)
  {
    try
    {
      Grow(FSize - FSizeOfArray + 10);
    }
    catch (Alloc &a)
    {
      MkDebug("MkLiuParticles::Add throws Alloc while trying to Grow memory");
      throw Alloc(a.what());
    }
  }
  if (FSize == FSizeOfArray)
    return false;

  for (int i = FSize - 1; i >= index; i--)
    FParticle[i + 1] = FParticle[i];
  FSize++;
  FParticle[index] = point;
  return true;
}

bool MkLiuParticles::Delete(MkLiuParticle point)
{
  int i;
  for (i = 0; i < FSize; i++)
  {
    if (FParticle[i] == point)
      break;
  }
  if (i == FSize)
    return false;
  if (FParticle[i] == point)
  {
    for (int j = i; j < FSize - 1; j++)
      FParticle[j] = FParticle[j + 1];
  }
  FSize--;
  FParticle[FSize] = NullParticle;
  return true;
}

bool MkLiuParticles::Delete(int index)
{
  for (int j = index; j < FSize - 1; j++)
    FParticle[j] = FParticle[j + 1];

  FSize--;
  FParticle[FSize] = NullParticle;
  return true;
}

bool MkLiuParticles::Swap(int i, int j)
{
  MkLiuParticle p;
  if (i >= FSize || j >= FSize)
    return false;
  p = FParticle[i];
  FParticle[i] = FParticle[j];
  FParticle[j] = p;
  return true;
}

bool MkLiuParticles::Clear()
{
  FSizeOfArray = FSize = 0;
  if (FParticle)
  {
    FParticle.reset();
  }
  return true;
}

MkLiuParticle &MkLiuParticles::operator[](int i)
{
  //    if (FSizeOfArray == 0) return NullParticle;
  //    if (i >= FSize && i < FSizeOfArray) FSize = i+1;

  if (0 <= i && i < FSize)
    return FParticle[i];
  else
  {
    MkDebug("MkLiuParticles index out of range\n");
    exit(-1); //return NullParticle;
  }
}

MkLiuParticles &MkLiuParticles::operator=(MkLiuParticles &points)
{
  int i;

  Clear();
  FSize = points.FSize;
  FSizeOfArray = points.FSizeOfArray;
  if (FSize == 0)
  {
    FParticle.reset();
    return *this;
  }

  try
  {
    FParticle.reset(new MkLiuParticle[FSizeOfArray]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuParticles::operator= throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }

  for (i = 0; i < FSize; i++)
    FParticle[i] = points.FParticle[i];
  for (i = FSize; i < FSizeOfArray; i++)
    FParticle[i] = NullParticle;

  return *this;
}

bool MkLiuParticles::operator==(MkLiuParticles &points)
{
  int i;

  if (FSize != points.FSize)
    return false;
  for (i = 0; i < FSize; i++)
    if (this->FParticle[i] != points.FParticle[i])
      return false;

  return true;
}

MkLiuParticles &MkLiuParticles::operator*=(MkMatrix4<double> &rm)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i] = this->FParticle[i] * rm;
  return *this;
}

MkLiuParticles &operator*(MkLiuParticles &rps, MkMatrix4<double> &rm)
{
  static MkLiuParticles rps_t;
  rps_t = rps;
  for (int i = 0; i < rps.FSize; i++)
    rps_t.FParticle[i] = rps.FParticle[i] * rm;
  return rps_t;
}

MkLiuParticles &MkLiuParticles::Translate(MkLiuParticle rp)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].Translate(rp);
  return *this;
}

MkLiuParticles &MkLiuParticles::Translate(double x, double y, double z)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].Translate(x, y, z);
  return *this;
}

MkLiuParticles &MkLiuParticles::Rotate(double alpha, double beta, double gamma)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].Rotate(alpha, beta, gamma);
  return *this;
}

MkLiuParticles &MkLiuParticles::RotateInX(double ang)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].RotateInX(ang);
  return *this;
}

MkLiuParticles &MkLiuParticles::RotateInY(double ang)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].RotateInY(ang);
  return *this;
}

MkLiuParticles &MkLiuParticles::RotateInZ(double ang)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].RotateInZ(ang);
  return *this;
}

MkLiuParticles &MkLiuParticles::RotateInA(double ang, double l, double m, double n)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].RotateInA(ang, l, m, n);
  return *this;
}

MkLiuParticles &MkLiuParticles::Scale(double sx, double sy, double sz)
{
  for (int i = 0; i < FSize; i++)
    this->FParticle[i].Scale(sx, sy, sz);
  return *this;
}

#ifdef __BCPLUSPLUS__
void MkLiuParticles::Draw(TObject *Sender)
{
  for (int i = 0; i < FSize; i++)
    (*this)[i].Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkLiuParticles::Draw(MkPaint *pb)
{
  for (int i = 0; i < FSize; i++)
    (*this)[i].Draw(pb);
}
#endif

#if defined(__GL_H__)
void MkLiuParticles::Draw()
{
  for (int i = 0; i < FSize; i++)
    (*this)[i].Draw();
}
#endif

//------------------------------------------------------------------------------
#pragma package(smart_init)
