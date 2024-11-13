//234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#include "MkLiuBound.hpp"

MkLiuBoundary NullBoundary(0);
MkLiuBoundarys NullBoundarys(0);
//double Zero = 0;

void MkLiuBoundary::Initialize()
{
  ColorType = ctNone;
  R = G = B = 0;
  Transparancy = 0;
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

  XVel = 0;
  YVel = 0;
  ZVel = 0;

  XAVel = 0;
  YAVel = 0;
  ZAVel = 0;

  BoundaryType = 0; // IType Boundary Type

  DUDt = 0; //**
  INDUDt = 0;
  AVDUDt = 0;
  AHDUDt = 0;
  DRhoDt = 0;
  DVXDt = 0; //**
  INDVXDt = 0;
  ARDVXDt = 0;
  ADDVXDt = 0;
  EXDVXDt = 0;
  DVYDt = 0; //**
  INDVYDt = 0;
  ARDVYDt = 0;
  ADDVYDt = 0;
  EXDVYDt = 0;
  DVZDt = 0; //**
  INDVZDt = 0;
  ARDVZDt = 0;
  ADDVZDt = 0;
  EXDVZDt = 0;
  TDSDt = 0;
}

void MkLiuBoundary::SetupColor()
{
}

#if defined(__GL_H__)
void MkLiuBoundary::Draw()
{
  /*  SetupColor();
  glPushMatrix();
    glTranslatef(FCenter.X*2-1,FCenter.Y*2-1,FCenter.Z*1);
    glColor3f (R, G, B);
    glPopMatrix();*/

  //  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // ApplyMatrix();
  glPushMatrix();
  //  glTranslatef(-(*this)[0].X,-(*this)[0].Y,-(*this)[0].Z);
  glColor3f(1, 1, 1);
  glLineWidth(10.0);

  glBegin(GL_LINES);
  glVertex3f(LD.X, LD.Y, LD.Z);
  glVertex3f(RD.X, RD.Y, RD.Z);
  glVertex3f(RU.X, RU.Y, RU.Z);
  glVertex3f(LU.X, LU.Y, LU.Z);
  glEnd();
  //  glFlush();
  glPopMatrix();
}
#endif

MkLiuBoundary &MkLiuBoundary::operator=(const MkLiuBoundary &lb)
{

  MkPlane::operator=((MkPlane &)lb);

#if defined(__GL_H__)
  R = lb.R;
  G = lb.G;
  B = lb.B;
  Transparancy = lb.Transparancy;
#endif

  Radius = lb.Radius;
  ColorType = lb.ColorType;

  XVel = lb.XVel;
  YVel = lb.YVel;
  ZVel = lb.ZVel;
  XAVel = lb.XAVel;
  YAVel = lb.YAVel;
  ZAVel = lb.ZAVel;
  BoundaryType = lb.BoundaryType;

  Mass = lb.Mass;             // real mass checked
  Rho = lb.Rho;               // density checked,
  Volume = lb.Volume;         // volume to calculate the selfdens (density when no neighbor particles)  mk
  Eta = lb.Eta;               // dynamic viscosity
  Press = lb.Press;           // P  checked
  Temp = lb.Temp;             // T  checked
  Energy = lb.Energy;         //U  checked
  SoundSpeed = lb.SoundSpeed; //C
  SmoothLen = lb.SmoothLen;   //Hsml checked

  CountIac = lb.CountIac;

  DUDt = lb.DUDt; //**
  INDUDt = lb.INDUDt;
  AVDUDt = lb.AVDUDt;
  AHDUDt = lb.AHDUDt;
  ADDUDt = lb.ADDUDt; // artificial drag internal energy
  APDUDt = lb.APDUDt; // artificial repel internal energy
  DRhoDt = lb.DRhoDt;
  DVXDt = lb.DVXDt; //**
  INDVXDt = lb.INDVXDt;
  ARDVXDt = lb.ARDVXDt;
  ADDVXDt = lb.ADDVXDt; // artificial drag
  APDVXDt = lb.APDVXDt; // artificial repel
  EXDVXDt = lb.EXDVXDt;
  DVYDt = lb.DVYDt; //**
  INDVYDt = lb.INDVYDt;
  ARDVYDt = lb.ARDVYDt;
  ADDVYDt = lb.ADDVYDt; // artificial drag
  APDVYDt = lb.APDVYDt; // artificail repel
  EXDVYDt = lb.EXDVYDt;
  DVZDt = lb.DVZDt; //**
  INDVZDt = lb.INDVZDt;
  ARDVZDt = lb.ARDVZDt;
  ADDVZDt = lb.ADDVZDt; // artificial drag
  APDVZDt = lb.APDVZDt; // artificial repel
  EXDVZDt = lb.EXDVZDt;
  TDSDt = lb.TDSDt;

  return (*this);
}
//-----------------------------------------------------------------------------
MkLiuBoundarys::MkLiuBoundarys(int size, MkLiuBoundary *jp)
{
  if (size < 0)
  {
    MkDebug("::MkLiuBoundarys - MkLiuBoundarys(int size)");
    throw Size(std::string("MkLiuBoundarys::Constructor throws mem allocation negative size exception"), size);
  }
  if (jp == NULL)
  {
    MkDebug("::MkLiuBoundarys - MkLiuBoundarys(int size)");
    throw Size(std::string("MkLiuBoundarys::Constructor throws Size exception due to NULL argument "), size);
  }

  FSize = size;
  if (FSize == 0)
  {
    FBoundary.reset();
    return;
  }

  try
  {
    FBoundary.reset(new MkLiuBoundary[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuBoundarys::Constructor throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }

  for (int i = 0; i < FSize; i++)
    FBoundary[i] = jp[i];
}

MkLiuBoundarys::MkLiuBoundarys(int size)
{
  if (size < 0)
  {
    MkDebug("::MkLiuBoundarys - MkLiuBoundarys(int size)");
    throw Size(std::string("MkLiuBoundarys::Constructor throws mem allocation negative size exception"), size);
  }

  FSize = size;
  if (FSize == 0)
  {
    FBoundary.reset();
    return;
  }

  try
  {
    FBoundary.reset(new MkLiuBoundary[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuBoundarys::Constructor throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }
}

MkLiuBoundarys::~MkLiuBoundarys()
{
  FSize = 0;
  FBoundary.reset();
}

void MkLiuBoundarys::Initialize(int size)
{
  if (size < 0)
  {
    MkDebug("::MkLiuBoundarys - Initialize(int size)");
    throw Size(std::string("MkLiuBoundarys::Constructor throws mem allocation negative size exception"), size);
  }
  if (FSize == size)
    return;

  FSize = size;
  if (FSize == 0)
  {
    FBoundary.reset();
    return;
  }

  try
  {
    FBoundary.reset(new MkLiuBoundary[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuBoundarys::Constructor throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }
}

void MkLiuBoundarys::Initialize(int size, MkLiuBoundary *jp)
{
  if (size < 0)
  {
    MkDebug("::MkLiuBoundarys - Initialize(int size)");
    throw Size(std::string("MkLiuBoundarys::Constructor throws mem allocation negative size exception"), size);
  }

  FSize = size;
  if (FSize == 0)
  {
    FBoundary.reset();
    return;
  }

  try
  {
    FBoundary.reset(new MkLiuBoundary[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuBoundarys::Constructor throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }

  for (int i = 0; i < FSize; i++)
    FBoundary[i] = jp[i];
}

bool MkLiuBoundarys::Clear()
{
  FSize = 0;
  FBoundary.reset();

  return true;
}

MkLiuBoundary &MkLiuBoundarys::operator[](int i)
{
  if (FSize == 0)
    return NullBoundary;
  else if (i >= 0 && i < FSize)
    return FBoundary[i];
  else
  {
    MkDebug("MkLiuBoundarys::operator[] range error");
    throw Range(std::string("MkLiuBoundarys::operator[] throws range exception"), i);
  }
}

MkLiuBoundarys &MkLiuBoundarys::operator=(MkLiuBoundarys &jps)
{
  int i;

  Clear();
  FSize = jps.FSize;
  if (FSize == 0)
  {
    FBoundary.reset();
    return *this;
  }

  try
  {
    FBoundary.reset(new MkLiuBoundary[FSize]);
  }
  catch (std::bad_alloc &a)
  {
    MkDebug("MkLiuBoundarys::Constructor throws bad_alloc!!!\n");
    throw Alloc(a.what());
  }

  for (i = 0; i < FSize; i++)
    this->FBoundary[i] = jps.FBoundary[i];

  return *this;
}

bool MkLiuBoundarys::operator==(MkLiuBoundarys &Reals)
{
  int i;

  if (FSize != Reals.FSize)
    return false;
  for (i = 0; i < FSize; i++)
    if (this->FBoundary[i] != Reals.FBoundary[i])
      return false;

  return true;
}

void MkLiuBoundarys::Translate(double x, double y, double z)
{
  for (int i = 0; i < GetSize(); i++)
    FBoundary[i].Translate(x, y, z);
}

void MkLiuBoundarys::Translate(MkPoint rp)
{
  for (int i = 0; i < GetSize(); i++)
    FBoundary[i].Translate(rp);
}

#ifdef __BCPLUSPLUS__
void MkLiuBoundarys::Draw(TObject *Sender)
{
  for (int i = 0; i < GetSize(); i++)
    FBoundary[i].Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkLiuBoundarys::Draw(MkPaint *pb)
{
  for (int i = 0; i < GetSize(); i++)
    FBoundary[i].Draw(pb);
}
#endif

#if defined(__GL_H__)
void MkLiuBoundarys::Draw()
{
  for (int i = 0; i < GetSize(); i++)
    FBoundary[i].Draw();
}
#endif

//---------------------------------------------------------------------------
#pragma package(smart_init)
