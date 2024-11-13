//---------------------------------------------------------------------------
#ifndef MkWallH
#define MkWallH
#include <string>
#include "MkEntity.h" 
//---------------------------------------------------------------------------
class MkWall : public MkEntity {
protected:
  MkLine Wall;
public:
  MkWall(){}
  MkWall(int){}
  ~MkWall(){}
public:
  void SetLine(MkLine &w)
    {
      Wall=w;
      Wall.SetFiniteness(true);
      Depth = Wall[0].Y<Wall[1].Y?-Wall[0].Y:-Wall[1].Y;
      Length = Wall.GetLength();
    }
public:
  MkLine &GetLine(){return Wall;}
  bool operator==(MkWall &);
  bool operator!=(MkWall &);
  MkWall & operator=(MkWall &);

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkWall");};
#else
  std::string ClassName(){return std::string("MkWall");}  
#endif

#ifdef __BCPLUSPLUS__
  void  Draw(TObject *Sender){Wall.Draw(Sender);}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *pb){Wall.Draw(pb);}
#endif
};

class MkWalls {
protected:
    MkWall *FWall;
    int FSize;//Actual size of walls
    int FSizeOfArray;
#ifdef __BCPLUSPLUS__
    TColor Color;
#endif
public:
    MkWalls(int size,MkWall *wall);
    MkWalls(int size);
    MkWalls(){FSizeOfArray = FSize = 0;FWall = NULL;}
     ~MkWalls();
    virtual void Initialize(int size);
    virtual void Initialize(int size,MkWall *);
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    bool Add(MkWall &wall);  // change of size of wall
    bool Add(int index,MkWall &wall);
    bool Delete(MkWall &wall);  // change of size of wall
    bool Delete(int index);
    int Grow(int Delta);            // change of size of array
    int Shrink(int Delta);          // change of size of array
    bool Clear();
#ifdef __BCPLUSPLUS__
    TColor GetColor(){return Color;};
    void SetColor(TColor c){Color = c;}
#endif

#ifdef __BCPLUSPLUS__
  void  Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif

    virtual MkWall & operator[](int);
    MkWalls & operator=(MkWalls &walls);
    bool operator==(MkWalls &walls);
};
//---------------------------------------------------------------------------
#endif



