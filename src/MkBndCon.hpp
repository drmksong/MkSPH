//---------------------------------------------------------------------------
#ifndef MkBndConH
#define MkBndConH
#include "MkDOF.h"
#include "MkRange.hpp"

//-------------------------------------------------------------
enum MkBndCondType {bcNone, bcDisX, bcDisY, bcDisZ, bcTemp, bcPress};

class MkBndCond {
protected:
  MkDOFs DOFs;
  MkRangeTree RangeTree;
public: //constructors
  MkBndCond();
  MkBndCond(int n);
  MkBndCond(MkDOFs &dofs, MkRangeTree &range);  
  ~MkBndCond();
public: //setting function
  void SetDOFs(MkDOFs &dofs){DOFs = dofs;}
  void SetRangeTree(MkRangeTree &range){RangeTree = range;}
public: //getting function
  MkDOFs & GetDOFs(){return DOFs;}
  MkRangeTree & GetRangeTree(){return RangeTree;}
  bool operator==(MkBndCond &bc){return DOFs==bc.DOFs && RangeTree == bc.RangeTree;}
  bool operator!=(MkBndCond &bc){return DOFs!=bc.DOFs || RangeTree != bc.RangeTree;}
#ifdef __BCPLUSPLUS__
  void Draw(MkPaintBox *){}
#endif
};

class MkBndConds {
protected:
  MkBndCond *FBndCond;
  int FSize;//Actual size of boundary conditions
  int FSizeOfArray;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif

public:
  MkBndConds(int size,MkBndCond *bc);
  MkBndConds(int size);
  MkBndConds(){FSize = 0;FBndCond = NULL;}
  ~MkBndConds();
  virtual void Initialize(int size);
  virtual void Initialize(int size,MkBndCond *);
  int GetSize(){return FSize;};
  int GetNumber(){return FSize;};
  bool Add(MkBndCond &bndcond);  // change of size of bndcond
  bool Add(int index,MkBndCond &bndcond);
  bool Delete(MkBndCond &bndcond);  // change of size of bndcond
  bool Delete(int index);
  int Grow(int Delta);            // change of size of array
  int Shrink(int Delta);          // change of size of array
  bool Clear();
#ifdef __BCPLUSPLUS__
  TColor GetColor(){return Color;};
  void SetColor(TColor c){Color = c;}
  virtual void Draw(TObject *);
#endif

  virtual MkBndCond & operator[](int);
  MkBndConds & operator=(MkBndConds &bndconds);
  bool operator==(MkBndConds &bndconds);
};
//---------------------------------------------------------------------------
extern MkBndCond NullBndCond;

#endif
