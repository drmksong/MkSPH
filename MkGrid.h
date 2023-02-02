//---------------------------------------------------------------------------
#ifndef MkGridH
#define MkGridH

#include "MkMisc.hpp"
#include "MkPoint.hpp"
#include "MkFloat.h"
#include "MkCube.hpp"
#include "MkFault.hpp"
#include "MkPlane.hpp"
#include "MkKrig.hpp"

class MkTopoGrid {
private:
    MkFloat X,Y,Elev;
    float XMin,XMax,YMin,YMax,ElevMin,ElevMax;//Z mean value.
    int NX,NY;
    bool Assigned;
public:
    MkTopoGrid();
    MkTopoGrid(MkFloat &x,MkFloat &y,MkFloat &elev);
    void SetX(MkFloat &x){X.CopyFrom(x);}
    void SetY(MkFloat &y){Y.CopyFrom(y);}
    void SetElev(MkFloat &elev){Elev.CopyFrom(elev);}
    void SetTopoGrid(MkFloat &x,MkFloat &y,MkFloat &elev);
    float operator()(float x,float y);
    float operator()(int i,int j){return Elev(i,j);}
    float GetX(int i){return X(i);}
    float GetY(int i){return Y(i);}
    MkTopoGrid & operator=(MkTopoGrid &tg);
    bool isAssigned(){return Assigned;}
};

class MkSPHGrid {
 public:
  float StartX, StartY, StartZ;
  float EndX, EndY, EndZ;
  MkInt ParticleRef;
  
  MkSPHGrid();
  ~MkSPHGrid();
  void Clear();
  void Setup(float sx,float sy, float sz, float ex, float ey, float ez);
  void Setup(int num);
  int GetNumOfParticle(void){return ParticleRef.getSzX();}
  bool Check(int parnum);
  bool Register(int parnum);
  bool Unregister(int parnum);
  int Out();
  int operator[](int n){return ParticleRef[n];}
};
//---------------------------------------------------------------------------
class MkSPHGrids {
protected:
    MkSPHGrid *FSPHGrid;
    int FSize;
public:
    MkSPHGrids(int size,MkSPHGrid *rl);
    MkSPHGrids(int FSize);
    MkSPHGrids(){FSize = 0;FSPHGrid = NULL;}
    ~MkSPHGrids();
    void Initialize(int size);
    void Initialize(int size,MkSPHGrid *rl);
    void Grow(int sz);
    void Add(MkSPHGrid &l){Grow(1);FSPHGrid[FSize-1] = l;}
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    MkSPHGrid * GetLine(){return FSPHGrid;}
    bool Clear();

    virtual MkSPHGrid & operator[](int);

    MkSPHGrids & operator=(MkSPHGrids &);

};
extern MkSPHGrid NullSPHGrid;
extern MkSPHGrids NullSPHGrids;

//---------------------------------------------------------------------------
#endif 
