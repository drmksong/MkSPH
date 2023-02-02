//---------------------------------------------------------------------------
#ifndef MkEntityH
#define MkEntityH

#ifdef __BCPLUSPLUS__
#include "GlobalVarUnit.h"
#endif

#include <string>
#include "MkShape.hpp"
#include "MkPolygon.hpp"
#include "MkRect.hpp"
#include "MkLine.hpp"
#include "assert.h"
//---------------------------------------------------------------------------

enum MkSide {mkLeft, mkRight, mkNone};

class MkEntity : public MkShape {
protected:
#ifdef __BCPLUSPLUS__
  TStringGrid *Grid;
  TTabSheet *Sheet;
  AnsiString Name;
#else
  std::string Name;
#endif
  MkDouble Prop;
  int Division;
  std::string className;
public:
  int Number;
  double Depth;
  double Length;
public:
  MkEntity()
    {
      Division = 1;
#ifdef __BCPLUSPLUS__
      Grid=NULL;
#endif
      className = "MkEntity";
    };
  MkEntity(int n)
    {
      Division = 1;
#ifdef __BCPLUSPLUS__
      Grid=NULL;
#endif
      className = "MkEntity";
    };
  ~MkEntity(){};
public: //setting function
#ifdef __BCPLUSPLUS__
  void SetGrid(TStringGrid* grid){Grid = grid;}
  void SetTabSheet(TTabSheet *sheet){Sheet = sheet;}
  void SetName(AnsiString name){Name=name;}
  void SetName(char *name){Name=name;}
#else
  void SetName(std::string name){Name=name;}
  void SetName(char * name){Name=name;}
#endif
  void SetEntNum(int n){Number=n;}
  void SetDivision(int n){Division=n;}
  void SetDepth(double d){Depth=d;}
  void SetLength(double l){Length=l;}
public: //getting function
  int  GetEntNum(){return Number;}
#ifdef __BCPLUSPLUS__
  AnsiString GetName(){return Name;}
#else
  std::string GetName(){return Name;}
#endif
  int  GetDivision(){return Division;}
  MkDouble &GetProp(){return Prop;}
  double GetDepth(void){return Depth;}
  double GetDepth(double x);
  double GetLength(){return Length;}
#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkEntity");}
  virtual bool UpdateFrom(TStringGrid*){return false;}
  virtual bool UpdateTo(TStringGrid*){return false;}
  virtual bool UpdateFrom(){return false;}
  virtual bool UpdateTo(){return false;}
  virtual void Out(TObject *){};
#else
  std::string ClassName(){return className;}
#endif
  virtual void Clear();
  virtual void Out(char *fname){};

#ifdef __BCPLUSPLUS__
  virtual void Draw(TObject *){MkDebug("pure virtual function MkEntity::Draw() is called\n");}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  virtual void Draw(MkPaint *){MkDebug("pure virtual function MkEntity::Draw() is called\n");}
#endif
  virtual bool operator==(MkEntity &){return true;}
  virtual bool operator!=(MkEntity &){return false;}
};

class MkEntities {
protected:
    MkEntity *FEntity;
    int FSize;//Actual size of entities
    int FSizeOfArray;
#ifdef __BCPLUSPLUS__
    TColor Color;
#endif
public:
    MkEntities(int size,MkEntity *ent);
    MkEntities(int size);
    MkEntities(){FSizeOfArray = FSize = 0;FEntity = NULL;}
     ~MkEntities();
    virtual void Initialize(int size);
    virtual void Initialize(int size,MkEntity *);
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    bool Add(MkEntity &entity);  // change of size of entity
    bool Add(int index,MkEntity &entity);
    bool Delete(MkEntity &entity);  // change of size of entity
    bool Delete(int index);
    int Grow(int Delta);            // change of size of array
    int Shrink(int Delta);          // change of size of array
    bool Clear();
#ifdef __BCPLUSPLUS__
    TColor GetColor(){return Color;};
    void SetColor(TColor c){Color = c;}
    virtual void Out(TObject *){};
#else

#endif
    virtual void Out(char *fname){};

#ifdef __BCPLUSPLUS__
  void Draw(TObject *){MkDebug("pure virtual function MkEntities::Draw() is called\n");}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *){MkDebug("pure virtual function MkEntities::Draw() is called\n");}
#endif

    virtual MkEntity & operator[](int);
    MkEntities & operator=(MkEntities &entities);
    bool operator==(MkEntities &entities);
};
//---------------------------------------------------------------------------
extern MkEntity NullEntity;
//---------------------------------------------------------------------------
#endif
