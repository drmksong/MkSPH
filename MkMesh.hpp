#ifndef MkMesh_H
#define MkMesh_H
//--------------------------------------------------------------------
#include <stdio.h>
#include <string>
#include "MkLayer.h"
#include "MkEntity.h"
#include "MkInt.h"
#include "MkFloat.h"
#include "MkPoint.h"
#include "MkLine.h"
#include "MkPlane.h"
#include "MkDOF.h"
#include "MkBndCon.h"

#ifndef dyn_cp
#define dyn_cp(_dest_,_src_) {\
  MkTriElement *_tri_,*_tr_;\
  MkBeamElement *_beam_,*_b_;\
  MkTrussElement *_truss_,*_t_;\
  MkLineElement *_line_,*_l_;\
  if(_tri_ = dynamic_cast<MkTriElement*>((_src_))) {\
    assert(_tri_);\
    MkTriElement *_tr_ = new MkTriElement();\
    *_tr_ = *_tri_;\
    (_dest_) = _tr_;\
  }\
  if(_beam_ = dynamic_cast<MkBeamElement*>((_src_))) {\
    assert(_beam_);\
    MkBeamElement *_b_ = new MkBeamElement();\
    *_b_ = *_beam_;\
    (_dest_) = _b_;\
  }\
  else if(_truss_ = dynamic_cast<MkTrussElement*>((_src_))) {\
    assert(_truss_);\
    MkTrussElement *_t_ = new MkTrussElement();\
    *_t_ = *_truss_;\
    (_dest_) = _t_;\
  }\
  else if(_line_ = dynamic_cast<MkLineElement*>((_src_))) {\
    assert(_line_);\
    MkLineElement *_l_ = new MkLineElement();\
    *_l_ = *_line_;\
    (_dest_) = _l_;\
  }\
}
#endif

//--------------------------------------------------------------------
// MkMesh.h
//
// Last Revised. 2004. Mar. 19.
//--------------------------------------------------------------------
enum MkNodeType {ntNone=0,  ntMech1D, ntMech2D, ntMech3D,
                            ntTher1D, ntTher2D, ntTher3D,
                            ntHydr1D, ntHydr2D, ntHydr3D,
                            ntTherMech1D, ntTherMech2D, ntTherMech3D,
                            ntHydrMech1D, ntHydrMech2D, ntHydrMech3D,
                            ntTherHydr1D, ntTherHydr2D, ntTherHydr3D,
                            ntTherHydrMech1D, ntTherHydrMech2D, ntTherHydrMech3D,
                            ntStrt1D, ntStrt2D, ntStrt3D,
                            ntTherStrt1D, ntTherStrt2D, ntTherStrt3D,ntGeneral};

enum MkResultType {rtNone=0, rtMoment, rtCompFrc, rtShearFrc, rtCompStr, rtShearStr, rtTensStr};
//---------------------------------------------------------------------------
class MkNode {
protected:
  bool Movable;
  MkPoint FPoint;
  MkNodeType FNodeType;
  MkDOFs FDOFs;
  MkFloat FPrimaryVar; // displacement, head, etc.

  MkFloat & GetPrimaryVar(){return FPrimaryVar;}
  float &PrimaryVar(int i){return FPrimaryVar(i);}
  float & operator[](int i){return FPrimaryVar(i);}
public:
  friend class MkNodes;
  friend class MkElement;
  friend class MkLineElement;
  friend class MkTrussElement;
  friend class MkBeamElement;
  MkNode();
  MkNode(int n){};
  friend float CalDist(MkNode &node1, MkNode &node2);
  MkPoint & GetPoint(){return FPoint;}
  friend float GetVoronoiRadius(MkNode &node1, MkNode &node2, MkNode &node3);
  // setting functing
  void SetNodeType(MkNodeType nt); // warning::it deletes previous data
  void SetDOF(int ndof);           // warning::it also deletes previous data
  void SetDOF(MkDOFs &dofs);       // warning::it also deletes previous data
  void SetPoint(MkPoint &pnt){FPoint = pnt;}
  void SetMovable(bool flag){Movable = flag;}
  void Reset();  // reset FPrimaryVar for the setting of FDOFs
  //getting function
  bool isMovable(){return Movable;}
  MkDOFs & GetDOFs(){return FDOFs;}
  bool Apply(MkBndConds &bc);
  void Out();

  void SetXDis(float v);
  void SetYDis(float v);
  void SetZDis(float v);
  void SetXAng(float v);
  void SetYAng(float v);
  void SetZAng(float v);

  float GetXDis();
  float GetYDis();
  float GetZDis();
  float GetXAng();
  float GetYAng();
  float GetZAng();

  bool isEq(MkNode &);
  bool isNE(MkNode &);
  bool operator==(MkNode &);
  bool operator!=(MkNode &);
  MkNode &operator=(MkNode &node)
  {
    Movable=node.Movable;
    FPoint=node.FPoint;
    FNodeType=node.FNodeType;
    FDOFs=node.FDOFs;
    FPrimaryVar.CopyFrom(node.FPrimaryVar); // displacement, head, etc.
    return *this;
  }

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

void Swap(MkNode &,MkNode &);

class MkNodes {
protected:
  MkNode *FNode;
  int FSize;//Actual size of nodes
  int FSizeOfArray;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif

public:
  MkNodes(int size,MkNode *node);
  MkNodes(int size);
  MkNodes(){FSizeOfArray = FSize = 0;FNode = NULL;}
  ~MkNodes();
  virtual void Initialize(int size);
  virtual void Initialize(int size,MkNode *);
  int GetSize(){return FSize;};
  int GetNumber(){return FSize;};
  bool Add(MkNode &node);  // change of size of node
  bool Add(MkNodes &node) {
       bool flag=true;
       MkDebug("MkNodes::Add(Nodes) is called.\n");
       for (int i=0;i<node.GetSize();i++)
         flag = Add(node[i]) && flag;
       return flag;
  }  // change of size of node
  bool Add(int index,MkNode &node);
  bool Delete(MkNode &node);  // change of size of node
  bool Delete(int index);
  bool Apply(MkBndConds &bc){bool flag=true; /*MkDebug("MkNodes::Apply is entered\n");*/for (int i=0;i<FSize;i++) flag = FNode[i].Apply(bc) && flag;/*MkDebug("MkNodes::Apply leaving\n");*/return flag;}
  int Grow(int Delta);            // change of size of array
  int Shrink(int Delta);          // change of size of array
  bool Clear();
  void Out(){for (int i=0;i<FSize;i++) FNode[i].Out();}
#ifdef __BCPLUSPLUS__
  TColor GetColor(){return Color;};
  void SetColor(TColor c){Color = c;}

#endif
  virtual MkNode & operator[](int);
  MkNodes & operator=(MkNodes &nodes);
  bool operator==(MkNodes &nodes);
  void IniPrimaryVar(int n){for (int i=0;i<FSize;i++) FNode[i].GetPrimaryVar().Initialize(n);}

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif

};
//--------------------------------------------------------------------
enum MkCrossType {ctNot, ctNodeToNode, ctEdgeToNode, ctNodeToEdge, ctIntersect}; 

class MkElement {
protected:
  MkInt ElemNode; // Global Node Number of Elemental Node;
  MkFloat Prop;
  MkNodes * NodeRef; // use for reference, do not allocate, deallocate it
  MkMatrix Stiff;    // local stifness matrix
  MkMatrix TranMat;  // local to global
  MkPoint Center, Centroid;
  bool haveResult;
  bool isStiffed;
  std::string className;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif

public:
#ifdef __BCPLUSPLUS__
  MkElement(){Center=NullPoint;Centroid=NullPoint;NodeRef=NULL;Color=clBlack;isStiffed=false;className = "MkElement";}
  MkElement(MkNodes &nodes){Center=NullPoint;Centroid=NullPoint;NodeRef=NULL;Color=clBlack;SetNodes(nodes);isStiffed=false;className = "MkElement";}
#else 
  MkElement(){Center=NullPoint;Centroid=NullPoint;NodeRef=NULL;isStiffed=false;}
  MkElement(MkNodes &nodes){Center=NullPoint;Centroid=NullPoint;NodeRef=NULL;SetNodes(nodes);isStiffed=false;}
#endif

  MkElement(int n){NodeRef=NULL;isStiffed=false;className = "MkElement";}
  ~MkElement(){}

public:
  virtual bool IsInside(int node){return false;};  // is node located inside the element
  virtual bool IsIn(MkNode &node){return false;}
  virtual bool IsIncluded(int node){return false;}; // is node included in ElemNode
  virtual bool IsCross(MkElement *elem){return false;}
  virtual void FindCenter();
  virtual void FindCentroid();
  virtual bool Post(){MkDebug("Pure virtual function MkElement::Post called()\n");return false;}
  virtual bool HaveResult(){return haveResult;}
public:
  virtual bool SetupStiff(){MkDebug("SetupStiff called on MkElement!");return false;}
  virtual bool SetupTrans(){return false;}
  virtual bool SetupProp(MkFloat &prop){Prop = prop;return false;}
  virtual bool SetAxialLoad(float axil){return false;}
  void SetNodes(MkNodes & nodes){NodeRef = &nodes;}
  virtual void SetElemNode(MkInt &n){ElemNode = n;}
  virtual bool SetJackingForce(float jacking){return false;}
  virtual bool SetIniForceCorrection(float iniforce){return false;}
  void SetStiffed(bool flag){isStiffed=flag;}
public:
  MkNodes & GetNodes(){return *NodeRef;}
  int GetNodeNumber(int i){if(i<ElemNode.getSzX()&&i>=0) return ElemNode(i);else return -1;}
  MkInt &GetElemNode(){return ElemNode;}
  MkNode & GetElemNode(int i){return (*NodeRef)[ElemNode(i)];}
  MkNode & GetGlobNode(int i){return (*NodeRef)[i];}
  MkPoint & GetCenter(){FindCenter(); return Center;}
  MkPoint & GetCentroid(){FindCentroid(); return Centroid;}
  MkMatrix & GetStiff(){return Stiff;}
  MkMatrix & GetTrans(){return TranMat;}
  int GetNumberOfNode(){return ElemNode.getSzX();}
  bool GetStiffed(){return isStiffed;}

  virtual float GetAxialForce(float xi){return 0;}
  virtual float GetAxialStress(float xi){return 0;}
  virtual float GetAxialStrain(float xi){return 0;}
  virtual float GetMoment(float xi){return 0;} // 0 < xi < 1
  virtual float GetShearForce(float xi){return 0;}
  virtual float GetCompStress(float xi){return 0;}
  virtual float GetShearStress(float xi){return 0;}
  virtual float GetTensStress(float xi){return 0;}
  virtual float GetDisp(int s){return 0;}
  virtual float GetAngDisp(int s){return 0;}
  virtual float GetMoment(int s){return 0;} //s:0 first, s:1 last
  virtual float GetShearForce(int s){return 0;} //s:0 first, s:1 last

  virtual MkCrossType GetCrossType(MkElement *){MkDebug("Pure virtual function MkElement::GetCrossType()\n");return ctNot;}
  virtual bool isElement(){return true;}
  virtual bool isLineElement(){return false;}
  virtual bool isTrussElement(){return false;}
  virtual bool isBeamElement(){return false;}
  virtual bool isTriElement(){return false;}
  virtual bool isQuadElement(){return false;}
  virtual bool isTetraElement(){return false;}
  virtual bool isPrismElement(){return false;}
  virtual bool isCubicElement(){return false;}
  virtual bool isGeneralEmenent(){return false;}
public:
  virtual bool operator==(MkElement &);
  virtual bool operator!=(MkElement &);
  virtual MkElement & operator=(MkElement &);
#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkElement");}
  bool isClass(AnsiString as){return as==ClassName();}
  virtual void Out(TMemo *);
#else
  std::string ClassName(){return className;}
  bool isClass(char* as){return className.compare(as);}
#endif

  virtual void Out(char *);// it always append data, so you must initialize the file.
  virtual void Out();

#ifdef __BCPLUSPLUS__
  virtual void Draw(TObject *){MkDebug("Pure virtual function MkElement::Draw() is called\n");}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  virtual void Draw(MkPaint *){MkDebug("Pure virtual function MkElement::Draw() is called\n");}
#endif
};

class MkLineElement : public MkElement {
protected:
  float Length;
  MkLine Line;
public:
  MkLineElement():MkElement(){Length=0;className="MkLineElement";}
  MkLineElement(MkNodes &nodes):MkElement(nodes){Length=0;className="MkLineElement";}
  MkLineElement(MkInt & elemnode, MkNodes &nodes);
  ~MkLineElement(){};
  bool SetupTrans();
  bool SetupStiff(){MkDebug("SetupStiff called on MkLineElement!");return false;}
  void SetElemNode(int n1, int n2){ElemNode.Initialize(2);ElemNode(0) = n1;ElemNode(1)=1;}
  void SetElemNode(MkInt &n){ElemNode = n;}
  bool IsIn(MkNode &node);
  bool IsInside(int node);  // is node located in this element
  bool IsIncluded(int node);
  bool IsCross(MkElement *elem);
  float GetLength();
  MkLine &GetLine()
    {
      MkPoint sp,ep;
      assert(ElemNode.getSzX()==2);
      assert(NodeRef);
      sp = GetElemNode(0).GetPoint();
      ep = GetElemNode(1).GetPoint();
      Line.SetLine(sp,ep);
      return Line;
    }
  MkCrossType GetCrossType(MkElement *);
  bool isLineElement(){return true;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkLineElement");}
#else
  std::string ClassName(){return className;}
#endif

  MkLineElement & operator=(MkLineElement &);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

class MkTrussElement : public MkLineElement {
 protected:
  float JackingForce,IniForceCorrection; //
 public:
  MkTrussElement():MkLineElement(){JackingForce = 0;className = "MkTrussElement";}
  MkTrussElement(MkNodes &nodes):MkLineElement(nodes){JackingForce = 0;className="MkTrussElement";}
  MkTrussElement(MkInt & elemnode, MkNodes &nodes):MkLineElement(elemnode,nodes){JackingForce=0;className="MkTrussElement";}
 public:
  bool Post();
 public:
  bool SetupStiff();
  bool SetJackingForce(float jacking){JackingForce = jacking;return true;}
  bool SetIniForceCorrection(float iniforce){IniForceCorrection = iniforce;return true;}
  bool SetupProp(float ym, float area, float jacking, float iniforce)
    {
      Prop.Initialize(5);
      Prop(0) = ym;
      Prop(1) = area;
      Prop(2) = jacking;
      JackingForce = jacking;
      IniForceCorrection = iniforce;
      return true;
    }
 public:
  float GetAxialForce(float xi);
  float GetAxialStress(float xi);
  float GetAxialStrain(float xi);
  float GetDisp(int s);
  float N1(float xi);
  float N2(float xi);
  float dN1(float xi);
  float dN2(float xi);

  bool isTrussElement(){return true;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkTrussElement");}
  void Out(TMemo *);
#else
  std::string ClassName(){return className;}
#endif

  void Out(char *);
  void Out();
  MkTrussElement & operator=(MkTrussElement &);

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
  virtual void DrawResult(TObject *,MkResultType rt);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
  virtual void DrawResult(MkPaint *,MkResultType rt);
#endif
};

class MkBeamElement : public MkLineElement {
protected:
  float AxialLoad, Height, Width, Thickness;
  MkVector FixedEnd,Resultant;
public:
  MkBeamElement() : MkLineElement() {AxialLoad = 0;FixedEnd.Initialize(12);Resultant.Initialize(12);className="MkBeamElement";}
  MkBeamElement(MkNodes &nodes) : MkLineElement(nodes){AxialLoad=0;FixedEnd.Initialize(12);Resultant.Initialize(12);className="MkBeamElement";}
  MkBeamElement(MkInt & elemnode, MkNodes & nodes):MkLineElement(elemnode,nodes)
    {AxialLoad = Height = Width = Thickness = 0;FixedEnd.Initialize(12);Resultant.Initialize(12);className="MkBeamElement";}
  ~MkBeamElement(){};
public:
  bool Post();
public:
  bool SetAxialLoad(float axil){AxialLoad = axil;return true;}
  bool SetupStiff();
  bool SetupProp(float ym, float area, float gj, float iy, float iz)
    {
      Prop.Initialize(5);
      Prop(0) = ym;
      Prop(1) = area;
      Prop(2) = gj;
      Prop(3) = iy;
      Prop(4) = iz;
      return true;
    }
public:
  float N1(float xi);
  float N2(float xi);
  float N3(float xi);
  float N4(float xi);
  float dN1(float xi);
  float dN2(float xi);
  float dN3(float xi);
  float dN4(float xi);
  float d2N1(float xi);
  float d2N2(float xi);
  float d2N3(float xi);
  float d2N4(float xi);
  float d3N1(float xi);
  float d3N2(float xi);
  float d3N3(float xi);
  float d3N4(float xi);
  MkVector & GetFixedEnd(){return FixedEnd;}
  MkVector & GetResultant(){return Resultant;}
  float GetDisp(float xi);
  float GetAngDis(float xi);
  float GetMoment(float xi); // 0 < xi < 1
  float GetShearForce(float xi);
  float GetCompStress(float xi);
  float GetShearStress(float xi);
  float GetTensStress(float xi);
  float GetDisp(int s);
  float GetAngDisp(int s);  
  float GetMoment(int s);//s:0 first, s:1 last
  float GetShearForce(int s); //s:0 first, s:1 last
  bool isBeamElement(){return true;}

public:
  MkBeamElement & operator=(MkBeamElement &);

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkBeamElement");}
  void Out(TMemo *);
#else
  std::string ClassName(){return className;}
#endif

  void Out(char *);
  void Out();
  
#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
  void DrawResult(TObject *,MkResultType rt);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
  void DrawResult(MkPaint *,MkResultType rt);
#endif
};

class MkTriElement : public MkElement { // Tri - 2Dim
protected:
  MkPoint FCenter;
  float Radius;
public:
  MkTriElement();
  ~MkTriElement(){};
  float Area(){return 0;};
  bool  IsIn(MkNode &node){return false;};
  void FindCenter(){};
#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkTriElement");}
#else
  std::string ClassName(){return className;}
#endif

  bool isTriElement(){return false;}

#ifdef __BCPLUSPLUS__
  void Draw(TObject *){};
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *){};
#endif
};

class MkQuadElement : public MkElement { //Quad - 2Dim
protected:  // NofNode = 4

public:
  MkQuadElement(){className = "MkQuadElement";};
  ~MkQuadElement(){};
  float Area(){return 0;};
  bool  IsIn(MkNode &node){return false;};
  void FindCenter(){};

  bool isQuadElement(){return true;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkQuadElement");}
#else
  std::string ClassName(){return className;}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *){};
#endif

};

class MkTetraElement : public MkElement { //Tetra - 3Dim
protected:  // NofNode = 4

public:
  MkTetraElement(){className = "MkTetraElement";};
  ~MkTetraElement(){};
  float Volume(){return 0;}
  bool  IsIn(MkNode &node){return false;};
  void FindCenter(){};

  bool isTetraElement(){return true;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkTetraElement");}
#else
  std::string ClassName(){return className;}
#endif

#ifdef __BCPLUSPLUS__
  void Draw(TObject *){};
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *){};
#endif

};

class MkPrismElement : public MkElement { //Prism - 3Dim
protected:  // NofNode = 6

public:
  MkPrismElement(){className = "MkPrismElement";};
  ~MkPrismElement(){};
  float Volume(){return 0;}
  bool  IsIn(MkNode &node){return false;};
  void FindCenter(){};
  bool isPrismElement(){return true;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkPrismElement");}
#else
  std::string ClassName(){return className;}
#endif

#ifdef __BCPLUSPLUS__
  void Draw(TObject *){};
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *){};
#endif

};

class MkCubicElement : public MkElement { //Cubic - 3Dim
protected:  // NofNode = 8

public:
  MkCubicElement(){className="MkCubicElement";};
  ~MkCubicElement(){};
  float Volume(){return 0;}
  bool IsIn(MkNode &node){return false;}
  void FindCenter(){};

  bool isCubicElement(){return true;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkCubicElement");}
#else
  std::string ClassName(){return className;}
#endif

#ifdef __BCPLUSPLUS__
  void Draw(TObject *){};
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *){};
#endif

};

class MkGenericElement : public MkElement {
protected:

public:
  MkGenericElement(){className = "MkGenericElement";};
  ~MkGenericElement(){};
  bool IsIn(MkNode &node){return false;};
  virtual bool IsIntersect(MkLine &rl){return false;};
  virtual MkPoint GetIntersection(MkLine &){return NullPoint;};
  virtual MkLine GetIntersection(MkPlane &){return NullLine;};

  bool isGeneralEmenent(){return true;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkGenericElement");}  
#else
  std::string ClassName(){return className;}  
#endif

#ifdef __BCPLUSPLUS__
  void Draw(TObject *){};
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *){};
#endif
};

class MkElements { // Element Container
protected:
  MkElement **FElement;
  int FSize,FSizeOfArray;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif
private:
  bool dynamic_copy(MkElement *dest, MkElement *src);
public:
  MkElements(int size,MkElement **elem);
  MkElements(){FSizeOfArray = FSize = 0;FElement = NULL;}
  ~MkElements();
  virtual void Initialize(int size,MkElement **);
  bool Add(MkElement *);
  bool AddNew(MkElement *);
  bool Delete(MkElement *);
  bool DeleteAElement(MkElement *elem);
  void SetNodes(MkNodes &nodes)
    {
      for(int i=0;i<FSize;i++) FElement[i]->SetNodes(nodes);
    }
  void SetupStiff()
    {                                       
      for(int i=0;i<FSize;i++)
       FElement[i]->SetupStiff();
    }
  int GetSize(){return FSize;};
  int GetNumber(){return FSize;};
  int Grow(int delta);
  int Shrink(int delta);
  bool Clear();
  MkElements & operator=(MkElements &);
  MkElement & operator[](int);
  bool operator==(MkElements &);
  bool operator!=(MkElements &);
#ifdef __BCPLUSPLUS__
  void Out(TObject *);
#else
  void Out(char *fname);
#endif
  void Out();

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};
/*-------------------------------------------------------------
class of MkBound is in charge of describing the model boundary,
and generate the nodes at the boundary. The other resposibilities
of class MkBound is to investigate the boundaries wether these
cover the entire domain or not. The container MkBounds of MkBound
grant the sequence of each boundaries.
--------------------------------------------------------------*/
enum MkBoundType {btNone=0, btOneDimension, btTwoDimension};

class MkMesh;

class MkBound {
protected:
  MkMesh *FMesh;
  MkBoundType FBoundType;
  int FBoundGroup;//Exterior boundary is group 0, Interior boundary starts 1 to...
  MkBound *NextBound;
  MkBound *PrevBound;
  char FName[256];
  std::string className;
public:
  MkBound();
  MkBound(int n);
  MkBound(MkBoundType bt);
  ~MkBound();
  std::string ClassName(){return className;}
  void SetBoundType(MkBoundType bt){FBoundType = bt;}
  MkBound *Next(){return NextBound;}
  MkBound *Prev(){return PrevBound;}
  void SetNext(MkBound *next){NextBound = next;}
  void SetPrev(MkBound *prev){PrevBound = prev;}
  void SetName(char *name){strcpy(FName,name);}
  char*GetName(){return FName;}
  virtual bool IsNeighborWith(MkBound &bnd);
  virtual bool GenerateMesh(){return false;};
  bool IsMeshGenerated(){return FMesh!=NULL;}
  MkMesh *GetMesh(){return FMesh;}
  virtual bool operator==(MkBound &bnd){return false;}
  virtual bool operator!=(MkBound &bnd){return true;}

#ifdef __BCPLUSPLUS__
  virtual void Draw(MkPaintBox *pb){};
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

class MkLineBound : public MkBound {
protected:
  MkNodes FNodes;
  MkLine FLine;
public:
  MkLineBound(){className = "MkLineBound";};
  ~MkLineBound();
  std::string ClassName(){return className;}
};

class MkCurveBound : public MkBound {
protected:
//  MkNurbsCurve FNurbsCurve;
public:
  MkCurveBound(){className = "MkCurveBound";};
  ~MkCurveBound();
  std::string ClassName(){return className;}
};

class MkPlaneBound : public MkBound {
protected:
  MkPoints FPlane;
  MkElements FElems;
public:
  MkPlaneBound(){className = "MkPlaneBound";};
  ~MkPlaneBound();
  std::string ClassName(){return className;}
};

class MkSurfBound : public MkBound {
protected:
//  MkNurbsSurf FNurbsSurf;
  MkElements FElems;
public:
  MkSurfBound(){className = "MkSurfBound";};
  ~MkSurfBound();
  std::string ClassName(){return className;}
};

class MkBounds {
protected:
    MkBound *FBound;
    int FSize;//Actual size of bounds
    int FSizeOfArray;
    MkNodes FNodes;
#ifdef __BCPLUSPLUS__
    TColor Color;
#endif

public:
    MkBounds(int size,MkBound *bound);
    MkBounds(int size);
    MkBounds(){FSize = 0;FBound = NULL;}
    ~MkBounds();
    virtual void Initialize(int size);
    virtual void Initialize(int size,MkBound *);
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    bool Add(MkBound &bound);  // change of size of bound
    bool Add(int index,MkBound &bound);
    bool Delete(MkBound &bound);  // change of size of bound
    bool Delete(int index);
    int Grow(int Delta);            // change of size of array
    int Shrink(int Delta);          // change of size of array
    bool Clear();
    virtual MkBound & operator[](int);
    MkBounds & operator=(MkBounds &bounds);
    bool operator==(MkBounds &bounds);

    bool IsClosed();    // check if it covers the entire boundary
    bool FindSequence();     // classify boundaries into several group and find the order of each group of boundaries
    bool GenBoundNodes();    // if aboves are ok, generate nodes on the boundaries with respect to FNumOfSeg.
#ifdef __BCPLUSPLUS__
    TColor GetColor(){return Color;};
    void SetColor(TColor c){Color = c;}
    virtual void Draw(TObject *);
#endif

};
//--------------------------------------------------------------------
class MkFront {
protected:
  MkNodes FFrontNodes;
  MkNodes FInternalNodes;
public:
  MkFront(){};
  ~MkFront(){};
  void SetFrontNodes(MkNodes &nodes);
  void SetInternalNodes(MkNodes &nodes);
  void GoFoward();
};
//---------------------------------------------------------------------------
enum MkMeshType {mtNone, mtOneDimension, mtTwoDimension, mtThreeDimension};

class MkMesh {
protected:
  MkElements Elements;
  MkNodes Nodes;
  MkLayers Layers;
  MkBounds Bounds;
  MkNodes BoundNodes;
  MkNodes InternalNodes;
  MkFront Front;
  bool    isGenerated;
  char    FName[256];
  float   Division;
  std::string className;
public:
  MkMesh();
  MkMesh(int n);
  ~MkMesh();
#ifdef __BCPLUSPLUS__
  virtual AnsiString ClassName(){return AnsiString("MkMesh");}
#else
  std::string ClassName(){return className;}
#endif
  virtual bool GenerateInternalNodes(){return false;};
  virtual bool Generate();
  virtual bool GenerateWithPolygon(MkPolygon &poly);

  virtual bool SetBoundary(MkBounds &bounds);
  void SetName(char *name){strcpy(FName,name);}
  void SetDivision(float div){Division = div;}
  void SetNodes(MkNodes &nodes){Nodes = nodes;}
  void SetElements(MkElements &elems){Elements = elems;Elements.SetNodes(Nodes);}
  
  char*GetName(){return FName;}
  MkElements & GetElements(){return Elements;}
  MkNodes & GetNodes(){return Nodes;}
  std::string GetClassName(){return className;}

public:
  bool operator==(MkMesh &){return false;};
  bool operator!=(MkMesh &){return true;};
  MkMesh &operator=(MkMesh&);
public:
#ifdef __BCPLUSPLUS__
  void Draw(MkPaintBox *pb){};
#endif

};

class MkMeshes {
protected:
    MkMesh *FMesh;
    int FSize;//Actual size of meshes
    MkNodes FNodes;
#ifdef __BCPLUSPLUS__
    TColor Color;
#endif
public:
    MkMeshes(int size,MkMesh *mesh);
    MkMeshes(int size);
    MkMeshes(){FSize = 0;FMesh = NULL;}
    ~MkMeshes();
    virtual void Initialize(int size);
    virtual void Initialize(int size,MkMesh *);
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    bool Clear();
    virtual MkMesh & operator[](int);
    MkMeshes & operator=(MkMeshes &meshes);
    bool operator==(MkMeshes &meshes);
#ifdef __BCPLUSPLUS__
    TColor GetColor(){return Color;};
    void SetColor(TColor c){Color = c;}
    virtual void Draw(TObject *);
#endif
};
extern MkNode NullNode;
extern MkElement NullElement;
extern MkBound NullBound;
extern MkMesh NullMesh;


#endif
