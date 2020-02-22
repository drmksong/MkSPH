//---------------------------------------------------------------------------
#pragma hdrstop
#include "stdafx.h"
#include "MkRange.h"

#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif

//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif
//---------------------------------------------------------------------------
MkRange::MkRange(MkRange *left, MkRange *right)
{
  //MkDebug("  Allocation of base is performed.\n");
  Parent=Left=Right=NULL;op=NULL;
  SetLeft(left);
  SetRight(right);
  className = "MkRange";
}

void MkRange::Clear()
{
  if(op) {
    delete (MkRangeOperator*)op;
    op=NULL;
  }
}

bool MkRange::SetParent(MkRange *parent)
{
  if(!parent) return false;

  if(parent->isShape()){
    Parent = (MkRangeShape*)parent;
    return true;
  }
  else if(parent->isOperator()) {
    Parent = (MkRangeOperator*)parent;
    return true;
  }
  return false;
}

bool MkRange::SetLeft(MkRange *left) // allocation is performed
{
  if(!left) return false;

  if(left->isShape()) {
    Left = (MkRangeShape*)left;
    return true;
  }
  else if(left->isOperator()) {
    Left = (MkRangeOperator*)left;
    return true;
  }
  return false;
}

bool MkRange::SetRight(MkRange *right)
{
  if(!right) return false;

  if(right->isShape()) {
    Right = (MkRangeShape*)right;
    return true;
  }
  else if(right->isOperator()) {
    Right = (MkRangeOperator*)right;
    return true;
  }
  return false;
}

MkRange *MkRange::GetChild()
{
  if(Left) return Left;
  else if(Right) return Right;
  else return NULL;
}

MkRange *MkRange::GetChild(int n)
{
  if(n==0) return Left; 
  else if(n==1) return Right;
  else return NULL;
}

MkRange *MkRange::operator[](int n)
{
  if(n==0) return Left; 
  else if(n==1) return Right;
  else return NULL;
}

bool MkRange::operator==(MkRange *ra)
{
  //MkDebug("Comparison between MkRange object occured.");
  return false;
}

bool MkRange::operator!=(MkRange *ra)
{
  //MkDebug("Comparison between MkRange object occured.");
  return false;
}

MkRangeOperator* MkRange::operator+(MkRange &right)
{
  if(op) {
    delete op;
    op = NULL;
  }
  op = new MkRangeOperator(this,&right,rotUnion);
  SetParent(op);
  right.SetParent(op);
  return op;
}

MkRangeOperator* MkRange::operator*(MkRange &right)
{
  if(op) {
    delete op;
    op = NULL;
  }
  op = new MkRangeOperator(this,&right,rotIntersect);
  SetParent(op);
  right.SetParent(op);
  return op;
}

MkRangeOperator* MkRange::operator-(MkRange &right)
{
  if(op) {
    delete op;
    op = NULL;
  }
  op = new MkRangeOperator(this,&right,rotDifference);
  SetParent(op);
  right.SetParent(op);
  return op;
}

MkRangeOperator* MkRange::operator!()
{
  if(op) {
    delete op;
    op = NULL;
  }
  op = new MkRangeOperator(this,NULL,rotComplement);
  SetParent(op);
  return op;
}
//---------------------------------------------------------------------------
MkRangeShape::MkRangeShape() : MkRange()
{
  //MkDebug("  Allocation of shape is performed.\n");
  Shape = NULL;
  Thick=0;
  MustBeInSpace=true;
  MustBeInSurface=false;
  className = "MkRangeShape";
}

MkRangeShape::MkRangeShape(MkRange *left, MkRange *right) : MkRange(left,right)
{
  //MkDebug("  Allocation of shape is performed.\n");
  Shape = NULL;
  Thick=0;
  MustBeInSpace=true;
  MustBeInSurface=false;
  className = "MkRangeShape";
}

MkRangeShape::MkRangeShape(MkRange *left, MkRange *right, MkShape *shape, bool inspace, bool insurface,float thick)
    : MkRange(left,right)
{
  SetShape(shape);
  Thick=thick;
  MustBeInSpace=inspace;
  MustBeInSurface=insurface;
  className = "MkRangeShape";
}

MkRangeShape::~MkRangeShape()
{
  if(Shape) {
    delete Shape;
    Shape = NULL;
  }
}

void MkRangeShape::SetShape(MkShape *shape)
{
  if(!shape) return;
  if(shape->isTriangle()) {
    if(Shape) delete Shape;
    MkTriangle *t = new MkTriangle();
    *t = *(MkTriangle*)shape;
    Shape = t;
  }
  else if(shape->isRect()) {
    if(Shape) delete Shape;
    MkRect *r= new MkRect();
    *r = *(MkRect*)shape;
    Shape = r;
  }
  else if(shape->isCube()) {
    if(Shape) delete Shape;
    MkCube *c = new MkCube();
    *c = *(MkCube*)shape;
    Shape = c;
  }
  else if(shape->isCircle()) {
    if(Shape) delete Shape;
    MkCircle *c = new MkCircle();
    *c = *(MkCircle*)shape;
    Shape = c;
  }
  else if(shape->isSphere()) {
    if(Shape) delete Shape;
    MkSphere *s = new MkSphere();
    *s = *(MkSphere*)shape;
    Shape = s;
  }
  else if(shape->isCylinder()) {
    if(Shape) delete Shape;
    MkCylinder *c = new MkCylinder();
    *c = *(MkCylinder*)shape;
    Shape = c;
  }
}

bool MkRangeShape::Operate(MkPoint &pnt)
{
  bool flag;
  //MkDebug("MkRangeShape::Operate is entering\n");

  if(!Shape) {
    //MkDebug("Shape is NULL\n");
    return false;
  }
  
//  //MkDebug("Shape is ");//MkDebug((char*)Shape->ClassName());
  //MkDebug("\n");
  flag = IsIn(pnt);
  //MkDebug("MkRangeShape::Operate is leaving\n");
  return flag;
}

bool MkRangeShape::IsIn(MkPoint &pnt)
{
  if (MustBeInSpace)
    return Shape->IsInSpace(pnt);
  else if (MustBeInSurface)
    return Shape->IsInSurface(pnt,Thick);
  return false;
}

MkRangeShape & MkRangeShape::operator=(MkRangeShape &rs)
{
  SetShape(rs.Shape);
  MustBeInSpace=rs.MustBeInSpace; // in case of plane, it is true to be in space for normal to plane
  MustBeInSurface=rs.MustBeInSurface;
  Thick=rs.Thick;
  return *this;
}

bool MkRangeShape::operator==(MkRangeShape &rs)
{
  return fabs(Thick-rs.Thick) < EPS && 
         MustBeInSpace == rs.MustBeInSpace && 
         MustBeInSurface == rs.MustBeInSurface &&
         *Shape==*rs.Shape;
}

bool MkRangeShape::operator!=(MkRangeShape &rs)
{
  return fabs(Thick-rs.Thick) > EPS || 
         MustBeInSpace != rs.MustBeInSpace || 
         MustBeInSurface != rs.MustBeInSurface || 
         *Shape!=*rs.Shape;
}

//---------------------------------------------------------------------------
MkRangeOperator::MkRangeOperator() : MkRange()
{
  //MkDebug("  Allocation of operator is performed.\n");
  RangeOperatorType=rotNone;
  className = "MkRangeOperator";
}

MkRangeOperator::MkRangeOperator(MkRange *left, MkRange *right, MkRangeOperatorType op) : MkRange(left,right)
{
  //MkDebug("  Allocation of operator is performed.\n");
  RangeOperatorType=op;
  className = "MkRangeOperator";
}

bool MkRangeOperator::Operate(MkPoint &pnt)
{
  bool lf=false,rf=false;
  //MkDebug("MkRangeOperator::Operate entering\n");

  if (Left&&Left->isShape()) {
     //MkDebug("MkRangeOperator::Operate 1\n");
     MkRangeShape *left = (MkRangeShape*)Left;
     lf = left->Operate(pnt);
  }
  else if (Left&&Left->isOperator()) {
     //MkDebug("MkRangeOperator::Operate 2\n");
     MkRangeOperator *left = (MkRangeOperator*)Left;
     lf = left->Operate(pnt);
  }

  if (Right&&Right->isShape()) {
     //MkDebug("MkRangeOperator::Operate 3\n");
     MkRangeShape *right = (MkRangeShape*)Right;
     rf = right->Operate(pnt);
     //MkDebug("MkRangeOperator::Operate 3-1\n");     
  }
  else if (Right&&Right->isOperator()) {
     //MkDebug("MkRangeOperator::Operate 4\n");
     MkRangeOperator *right = (MkRangeOperator*)Right;
     rf = right->Operate(pnt);
  }

  //MkDebug("MkRangeOperator::Operate ^^\n");

  switch(RangeOperatorType) {
    case rotUnion: //MkDebug("MkRangeOperator::Operate 5\n");
                   return lf || rf;
    case rotIntersect: //MkDebug("MkRangeOperator::Operate 6\n");
                   return lf && rf;
    case rotDifference: //MkDebug("MkRangeOperator::Operate 7\n");
                   return lf && !rf;
    case rotComplement: //MkDebug("MkRangeOperator::Operate 8\n");
                   return Left?!lf:(Right?!rf :false);
  }
  //MkDebug("MkRangeOperator::Operate leaving\n");
  return false;
}

bool MkRangeOperator::IsIn(MkPoint &pnt)
{
  return Operate(pnt);
}

MkRangeOperator &MkRangeOperator::operator=(MkRangeOperator &ro)
{
  RangeOperatorType = ro.RangeOperatorType;
  return *this;
}

bool MkRangeOperator::operator==(MkRangeOperator &ro)
{
  return RangeOperatorType == ro.RangeOperatorType;
}

bool MkRangeOperator::operator!=(MkRangeOperator &ro)
{
  return RangeOperatorType != ro.RangeOperatorType;
}

//---------------------------------------------------------------------------
MkRangeTree::MkRangeTree()
{
  Root = NULL;
}

MkRangeTree::~MkRangeTree()
{
  Clear();
}

bool MkRangeTree::Verify()
{
  if(!Root) return false;
  return Verify(Root);
}

bool MkRangeTree::Verify(MkRange *parent)
{
  bool flag;
  MkRange *left,*right;
  
  assert(parent);
  left = parent->GetLeft();
  right = parent->GetRight();
  if (parent->GetOperator()==rotUnion
    ||parent->GetOperator()==rotIntersect
    ||parent->GetOperator()==rotDifference)
    if (!left || !right) return false;
  else if(parent->GetOperator()==rotComplement)       
    if ((left || !right) && (!left || right)) return false;
  else if(parent->isShape())
    if (left || right) return false;
    else return true;

  if(left) flag = Verify(left);
  if(right) flag = flag && Verify(right);
  return flag;
}

bool MkRangeTree::Clear()
{
  DeleteChild(Root);
  return true;
}

void MkRangeTree::DeleteChild(MkRange*&parent)
{
  if (!parent) return;
  if(parent->Left)
    DeleteChild(parent->Left);
  if(parent->Right)
    DeleteChild(parent->Right);
  if(!parent->GetChild()) {
    //MkDebug("Delete Range\n");
    if(parent) {
      if(parent->isShape()) {
	MkRangeShape *s=(MkRangeShape*)parent;
	delete s;
      }
      else if(parent->isOperator()) {
	MkRangeOperator *o=(MkRangeOperator*)parent;
	delete o;
      }
      else delete parent;
    }
    parent = NULL; 
    return;
  }
}

MkRangeTree &MkRangeTree::operator=(MkRangeTree &rt)
{
  if(!CopyChild(Root,rt.Root)) Clear();
  return *this;
}

bool MkRangeTree::operator==(MkRangeTree &rt)
{
  return CompChild(Root,rt.Root);
}

bool MkRangeTree::operator!=(MkRangeTree &rt)
{
  return !CompChild(Root,rt.Root);
}

bool MkRangeTree::CompChild(MkRange *ra,MkRange *rb)
{
  bool flag=true;
  if(!ra || !rb) return false;
  if (*ra!=rb) return false;
  if (Xor(ra->Left,rb->Left)) return false;
  if (Xor(ra->Right,rb->Right)) return false;
  for(int i=0;i<2;i++) if(!CompChild(ra->GetChild(i),rb->GetChild(i))) return false;
  return true;
}

bool MkRangeTree::CopyChild(MkRange *&ra, MkRange *rb)
{
  bool flag=true;
  MkRange *parent;
  MkRangeShape *rs;
  MkRangeOperator *ro;
  //MkDebug(" Enter the copy child.  \n");

  if(!rb) return false;

  if(ra) {
    if(ra->Parent) {
      //MkDebug(" ra->Parent is not null  \n");
      if(ra->Parent->isShape()) {
	//MkDebug("   ra->Parent is shape  \n");
	parent = (MkRangeShape*)ra->Parent;
      }
      else if(ra->Parent->isOperator()) {
	//MkDebug("   ra->Parent is operator  \n");
	parent = (MkRangeOperator*)ra->Parent;
      }
    }
  }
  else {
    //MkDebug(" ra is null  \n");
    parent = NULL;
  }

  if(ra) DeleteChild(ra);

  if(rb->isShape()) {
    //MkDebug("   Shape is copied.\n");
    rs = new MkRangeShape();
    if(!rs) return false;
    *rs = *(MkRangeShape*)rb;
    rs->Parent = parent;
    rs->Left = NULL;
    rs->Right = NULL;
    ra = rs;
  }
  else if(rb->isOperator()) {
    //MkDebug("   Operator is copied.\n");
    ro = new MkRangeOperator();
    if(!ro) return false;
    *ro = *(MkRangeOperator*)rb;
    ro->Parent = parent;
    ro->Left = NULL;
    ro->Right = NULL;
    ra = ro;
  }

  if(rb->Left) flag = CopyChild(ra->Left,rb->Left) && flag;
  if(rb->Right) flag = CopyChild(ra->Right,rb->Right) && flag;
  //MkDebug(" Leaving the copy child.  \n");
  return flag;
}

bool Xor(MkRange *a, MkRange *b)
{
  return (a&&b)!=(a||b);
}

//---------------------------------------------------------------------------

