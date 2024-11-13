//---------------------------------------------------------------------------
#ifndef MkRangeH
#define MkRangeH
#include "MkPoint.hpp"
#include "MkShape.hpp"
#include "MkCube.hpp"
#include "MkSphere.hpp"
#include "MkPlane.hpp"
#include "MkCylinder.hpp"
#include "MkRect.hpp"
#include "MkCircle.hpp"

//---------------------------------------------------------------------------
enum MkRangeOperatorType { rotNone, rotUnion, rotIntersect, rotComplement, rotDifference }; // ������(+),������(*),������(!),������(-)
class MkRangeOperator;
class MkRange {
  friend class MkRangeTree;
protected:
  MkRange *Parent;
  MkRange *Left;
  MkRange *Right;
  MkRangeOperator *op;
  std::string className;
public:
  MkRange(){MkDebug("  Allocation of base is performed.\n");op=NULL;Parent=Left=Right=NULL;className = "MkRange";}
  MkRange(MkRange *left, MkRange *right);
  ~MkRange(){
    Clear();
  }
  void Clear();
  virtual bool Operate(MkPoint &pnt){return false;}
  virtual bool IsIn(MkPoint &pnt){return false;}
  bool SetParent(MkRange *parent);
  bool SetLeft(MkRange *left);
  bool SetRight(MkRange *right);
  MkRange *GetLeft(){return Left;}
  MkRange *GetRight(){return Right;}
  MkRange *GetParent(){return Parent;}
  MkRange *GetChild();
  MkRange *GetChild(int n);
  MkRange *operator[](int n);
  virtual bool operator==(MkRange *);
  virtual bool operator!=(MkRange *);

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkRange");}
#else
  std::string ClassName(){return className;}
#endif

  virtual bool isShape(){return false;}
  virtual bool isOperator(){return false;}
  virtual MkRangeOperatorType GetOperator(){return rotNone;}
  MkRangeOperator* operator+(MkRange &);
  MkRangeOperator* operator*(MkRange &);
  MkRangeOperator* operator-(MkRange &);
  MkRangeOperator* operator!();
};

class MkRangeShape : public MkRange {
protected:
  MkShape *Shape;
  bool MustBeInSpace; // in case of plane, it is true to be in space for normal to plane
  bool MustBeInSurface;
  float Thick;
public:
  MkRangeShape();
  MkRangeShape(MkRange *left, MkRange *right);
  MkRangeShape(MkRange *left, MkRange *right, MkShape *shape, bool inspace=true, bool insurface=false,float thick=0);
  ~MkRangeShape();
  void SetShape(MkShape *shape);
  bool Operate(MkPoint &pnt);
  bool IsIn(MkPoint &pnt);
  virtual bool isShape(){return true;}

  MkRangeShape &operator=(MkRangeShape &);
  bool operator==(MkRangeShape &);
  bool operator!=(MkRangeShape &);

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkRangeShape");}
#else
  std::string ClassName(){return className;}
#endif
};

class MkRangeOperator : public MkRange {
protected:
  MkRangeOperatorType RangeOperatorType;
public:
  MkRangeOperator();
  MkRangeOperator(MkRange *left, MkRange *right=NULL, MkRangeOperatorType op=rotNone);
  ~MkRangeOperator(){};
  void SetOperator(MkRangeOperatorType type){RangeOperatorType = type;}
  MkRangeOperatorType GetOperator(){return RangeOperatorType;}
  bool Operate(MkPoint &pnt);
  bool IsIn(MkPoint &pnt);
  virtual bool isOperator(){return true;}

  MkRangeOperator &operator=(MkRangeOperator &ro);
  bool operator==(MkRangeOperator &);
  bool operator!=(MkRangeOperator &);

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkOperator");};
#else
  std::string ClassName(){return className;}
#endif
};

class MkRangeTree {
friend class MkRange;
friend class MkRangeShape;
friend class MkRangeOperator;
protected:
  MkRange *Root;
public:
  MkRangeTree();
  ~MkRangeTree();
  bool Verify();
  bool Verify(MkRange *parent);
  bool Clear();
  void SetRoot(MkRange *root){Clear(); CopyChild(Root,root);} // memory allocation occured
  void DeleteChild(MkRange*&parent);
  bool CompChild(MkRange *ra, MkRange *rb);
  bool CopyChild(MkRange *&ra, MkRange *rb);
  bool Operate(MkPoint &pnt){return Root->Operate(pnt);}
  MkRangeTree &operator=(MkRangeTree &rt); // memory allocated, make a complete clone
  bool operator==(MkRangeTree &rt);
  bool operator!=(MkRangeTree &rt);
};
bool Xor(MkRange *a, MkRange *b);
#endif
