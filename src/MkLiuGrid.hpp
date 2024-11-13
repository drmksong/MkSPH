#ifndef LiuGRID_HPP
#define LiuGRID_HPP

#include <time.h>
#include <unistd.h>
#include "MkArray.hpp"
#include "MkLiuParticle.hpp"

class MkLiuGrid
{
public:
  float StartX, StartY, StartZ;
  float EndX, EndY, EndZ;
  MkInt ParticleRef;
  int NParticles;

  MkLiuGrid();
  MkLiuGrid(int i) { Initialize(); }
  ~MkLiuGrid();
  void Initialize();
  void Clear();
  void Setup(float sx, float sy, float sz, float ex, float ey, float ez);
  void Setup(int num);
  int CountNParticles();
  int GetNumOfParticle(void) { return NParticles; }
  MkInt &GetParticleRef(void) { return ParticleRef; }
  bool Check(int parnum);
  bool Register(int parnum);
  bool Unregister(int parnum);
  int Out();
  int operator[](int n) { return ParticleRef[n]; }
};
//---------------------------------------------------------------------------
class MkLiuGrids
{
protected:
  boost::shared_array<MkLiuGrid> FLiuGrid;
  int FSize;
  int NX, NY, NZ;
  double FSmoothLen;
  double GridLen;
  double XMin, YMin, ZMin, XMax, YMax, ZMax;

public:
  MkLiuGrids();
  MkLiuGrids(int size, MkLiuGrid *rl);
  MkLiuGrids(int FSize);
  ~MkLiuGrids();
  void Initialize(int size);
  void Initialize(int size, MkLiuGrid *rl);
  void Clear();
  void Setup();
  void Setup(MkLiuParticles &par);
  void Update(MkLiuParticles &par);
  void Update_backup(MkLiuParticles &par);
  void Grow(int sz);
  void Add(MkLiuGrid &l)
  {
    Grow(1);
    FLiuGrid[FSize - 1] = l;
  }
  void SetSmoothLen(double sm) { FSmoothLen = sm; }
  void SetXMin(double xm)
  {
    XMin = xm;
    MkDebug("SetXMin \n");
  }
  void SetYMin(double ym)
  {
    YMin = ym;
    MkDebug("SetYMin \n");
  }
  void SetZMin(double zm)
  {
    ZMin = zm;
    MkDebug("SetZMin \n");
  }
  void SetXMax(double xm)
  {
    XMax = xm;
    MkDebug("SetXMax \n");
  }
  void SetYMax(double ym)
  {
    YMax = ym;
    MkDebug("SetYMax \n");
  }
  void SetZMax(double zm)
  {
    ZMax = zm;
    MkDebug("SetZMax \n");
  }
  int GetSize() { return FSize; }
  int GetNumber() { return FSize; }
  boost::shared_array<MkLiuGrid> GetLine() { return FLiuGrid; }
  double GetSmoothLen() { return FSmoothLen; }
  double GetXMin() { return XMin; }
  double GetYMin() { return YMin; }
  double GetZMin() { return ZMin; }
  double GetXMax() { return XMax; }
  double GetYMax() { return YMax; }
  double GetZMax() { return ZMax; }

  int GetNX() { return NX; }
  int GetNY() { return NY; }
  int GetNZ() { return NZ; }

  int GetNumOfParticles(int I, int J, int K);
  MkInt &GetParticles(int I, int J, int K);

  virtual MkLiuGrid &operator[](int);
  MkLiuGrid &operator()(int I, int J, int K);

  MkLiuGrids &operator=(MkLiuGrids &);

  class Alloc
  {
  public:
    std::string What;
    Alloc(std::string what) : What(what) {}
    std::string what() { return What; }
  };
  class Size
  {
  public:
    std::string What;
    int N;
    Size(std::string what, int n) : What(what), N(n) {}
    std::string what() { return What; }
  };
  class Range
  {
  public:
    std::string What;
    int N;
    Range(std::string what, int n) : What(what), N(n) {}
    std::string what() { return What; }
  };
};

extern MkLiuGrid NullLiuGrid;
#endif

/*

class MkLiuGrid {
protected:
  int MaxNGridX;
  int MaxNGridY;
  int MaxNGridZ;
  double SmoothLen;

  int NGridX[3]; // number of sorting grid cells in x,y, z direction
  double MaxGridX[3]; // maximum x,y,z coordinate
  double MinGridX[3]; // maximum x,y,z coordinate
  double DGeomX[3];
  double GHsmlX[3];
  MkInt Grid;

public:
  MkLiuGrid(){Initialize();}
  ~MkLiuGrid(){Clear();}

  void Initialize(int nx, int ny, int nz);
  void Initialize(){Clear();}
  void Clear();
  void Dump(){}

  void SetMaxNGridX(int m){MaxNGridX = m;}
  void SetMaxNGridY(int m){MaxNGridY = m;}
  void SetMaxNGridZ(int m){MaxNGridZ = m;}
  void SetSmoothLen(double sl){SmoothLen = sl;}

  void SetNGridX(int n, int g){(0<=n&&n<2)?NGridX[n] = g: MkDebug("in MkLiuGrid::SetNGridX Invalid Index \n");}
  void SetMaxGridX(int n,int m){(0<=n&&n<2)?MaxGridX[n] = m: MkDebug("in MkLiuGrid::SetMaxGridX Invalid Index \n");} 
  void SetMinGridX(int n,int m){(0<=n&&n<2)?MinGridX[n] = m: MkDebug("in MkLiuGrid::SetMinGridX Invalid Index \n");} 
  void SetDGeomX(int n,int m){(0<=n&&n<2)?DGeomX[n] = m: MkDebug("in MkLiuGrid::SetDGeomX Invalid Index \n");} 
  void SetGHsmlX(int n,int m){(0<=n&&n<2)?GHsmlX[n] = m: MkDebug("in MkLiuGrid::SetGHsmlX Invalid Index \n");} 

  void SetGrid(int i, int j, int k, int g){Grid(i,j,j) = g;}

  int GetMaxNGridX(){return MaxNGridX;}
  int GetMaxNGridY(){return MaxNGridY;}
  int GetMaxNGridZ(){return MaxNGridZ;}
  double GetSmoothLen(){return SmoothLen;}

  int GetNGridX(int n, int g){(0<=n&&n<2)?NGridX[n] = g: MkDebug("in MkLiuGrid::GetNGridX Invalid Index \n");}
  int GetMaxGridX(int n,int m){(0<=n&&n<2)?MaxGridX[n] = m: MkDebug("in MkLiuGrid::GetMaxGridX Invalid Index \n");} 
  int GetMinGridX(int n,int m){(0<=n&&n<2)?MinGridX[n] = m: MkDebug("in MkLiuGrid::GetMinGridX Invalid Index \n");} 
  int GetDGeomX(int n,int m){(0<=n&&n<2)?DGeomX[n] = m: MkDebug("in MkLiuGrid::GetDGeomX Invalid Index \n");} 
  int GetGHsmlX(int n,int m){(0<=n&&n<2)?GHsmlX[n] = m: MkDebug("in MkLiuGrid::GetGHsmlX Invalid Index \n");} 

  int GetGrid(int i, int j, int k, int g){Grid(i,j,j) = g;}


};

 */
