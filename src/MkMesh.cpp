//--------------------------------------------------------------------
// MkMesh.cpp
// 
// Last Revised. 2004. Feb. 11
//--------------------------------------------------------------------
#include "MkMesh.h"

MkNode NullNode(0);
MkElement NullElement(0);
MkBound NullBound(0);
MkMesh NullMesh(0);

#ifdef __BCPLUSPLUS__
#endif

//--------------------------------------------------------------------
MkNode::MkNode()
{
  FNodeType=ntNone;
  Movable = true;
}

float CalDist(MkNode &node1, MkNode &node2)
{
  return CalDist(node1.GetPoint(),node2.GetPoint());
}

float GetVoronoiRadius(MkNode &node1, MkNode &node2, MkNode &node3)
{
  return -1;
}

void MkNode::SetNodeType(MkNodeType nt) // warning::it deletes previous data
{
  int FNumDOF;
  if(FNodeType==nt) return;
  FNodeType = nt;

  if (FNodeType==ntNone ||FNodeType == ntGeneral) {
     FDOFs.Clear();
     FPrimaryVar.Clear();
     return;
  }

  switch(FNodeType) {
    case ntMech1D: FNumDOF=1;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   break;
    case ntMech2D: FNumDOF=2;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   break;
    case ntMech3D: FNumDOF=3;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftZDis,bndtFree);
                   break;
    case ntTher1D:
    case ntTher2D:
    case ntTher3D: FNumDOF=1;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftTemp,bndtFree);
                   break;
    case ntHydr1D:
    case ntHydr2D:
    case ntHydr3D: FNumDOF=1;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftHydroPress,bndtFree);
                   break;
    case ntTherMech1D: FNumDOF=2;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftTemp,bndtFree);
                   break;
    case ntTherMech2D: FNumDOF=3;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftTemp,bndtFree);
                   break;
    case ntTherMech3D: FNumDOF=4;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftZDis,bndtFree);
                   FDOFs[3].SetType(doftTemp,bndtFree);
                   break;
    case ntHydrMech1D: FNumDOF=2;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftHydroPress,bndtFree);
                   break;
    case ntHydrMech2D: FNumDOF=3;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftHydroPress,bndtFree);
                   break;
    case ntHydrMech3D: FNumDOF=4;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftZDis,bndtFree);
                   FDOFs[3].SetType(doftHydroPress,bndtFree);
                   break;
    case ntTherHydr1D:
    case ntTherHydr2D:
    case ntTherHydr3D: FNumDOF=2;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftTemp,bndtFree);
                   FDOFs[1].SetType(doftHydroPress,bndtFree);
                   break;
    case ntTherHydrMech1D: FNumDOF=3;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftTemp,bndtFree);
                   FDOFs[2].SetType(doftHydroPress,bndtFree);
                   break;
    case ntTherHydrMech2D: FNumDOF=4;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftTemp,bndtFree);
                   FDOFs[3].SetType(doftHydroPress,bndtFree);
                   break;
    case ntTherHydrMech3D: FNumDOF=5;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftZDis,bndtFree);
                   FDOFs[3].SetType(doftTemp,bndtFree);
                   FDOFs[4].SetType(doftHydroPress,bndtFree);
                   break;
    case ntStrt1D: FNumDOF=2;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftXAng,bndtFree);
                   break;
    case ntStrt2D: FNumDOF=4;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftXAng,bndtFree);
                   FDOFs[3].SetType(doftYAng,bndtFree);
                   break;
    case ntStrt3D: FNumDOF=6;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftZDis,bndtFree);
                   FDOFs[3].SetType(doftXAng,bndtFree);
                   FDOFs[4].SetType(doftYAng,bndtFree);
                   FDOFs[5].SetType(doftZAng,bndtFree);
                   break;
    case ntTherStrt1D: FNumDOF=3;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftXAng,bndtFree);
                   FDOFs[2].SetType(doftTemp,bndtFree);
                   break;
    case ntTherStrt2D: FNumDOF=5;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftXAng,bndtFree);
                   FDOFs[3].SetType(doftYAng,bndtFree);
                   FDOFs[4].SetType(doftTemp,bndtFree);
                   break;
    case ntTherStrt3D: FNumDOF=7;
                   FDOFs.Initialize(FNumDOF);
                   FPrimaryVar.Initialize(FNumDOF);
                   FDOFs[0].SetType(doftXDis,bndtFree);
                   FDOFs[1].SetType(doftYDis,bndtFree);
                   FDOFs[2].SetType(doftZDis,bndtFree);
                   FDOFs[3].SetType(doftXAng,bndtFree);
                   FDOFs[4].SetType(doftYAng,bndtFree);
                   FDOFs[5].SetType(doftZAng,bndtFree);
                   FDOFs[6].SetType(doftTemp,bndtFree);
                   break;
    default:       break;
  };
}

void MkNode::SetDOF(int ndof) // warning::it also deletes previous data
{
  FNodeType = ntGeneral;
  FDOFs.Initialize(ndof);
  FPrimaryVar.Initialize(ndof);
}

void MkNode::SetDOF(MkDOFs &dofs) // warning::it also deletes previous data
{
  FNodeType = ntGeneral;
  FDOFs = dofs;
  FPrimaryVar.Initialize(FDOFs.GetSize());
}

bool MkNode::Apply(MkBndConds &bc)
{
  int i,j,k;
  bool flag=false;
  MkDebug("MkNode::Apply is entered\n");
  for(k=0;k<bc.GetSize();k++) {
    MkDOFs &dofs=bc[k].GetDOFs();
    if (!bc[k].GetRangeTree().Operate(FPoint)) continue;
    MkDebug("MkNode::Apply before dof\n");
    for(i=0;i<dofs.GetSize();i++)
      for (j=0;j<FDOFs.GetSize();j++)
	if(dofs[i].GetDOFType() == FDOFs[j].GetDOFType()) {
	  FDOFs[j].SetBNDType(dofs[i].GetBNDType());
	  flag=true;
	}
    MkDebug("MkNode::Apply after dof\n");
  }
  MkDebug("MkNode::Apply leaving\n");
  return flag;
}

void MkNode::Out()
{
  char str[256];
  sprintf(str,"node coordinate (%10.5f, %10.5f, %10.5f)\n",FPoint.X, FPoint.Y, FPoint.Z);
  MkDebug(str);
  //  FDOFs.Out();
}

void MkNode::Reset()  // reset FPrimaryVar for the setting of FDOFs
{                     // clear all the value of FPrimaryVar and resize
  FPrimaryVar.Clear();
  FPrimaryVar.Initialize(FDOFs.GetSize());
}

void MkNode::SetXDis(float v)
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftXDis)
      FPrimaryVar(i) = v;
}

void MkNode::SetYDis(float v)
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftYDis)
      FPrimaryVar(i) = v;
}

void MkNode::SetZDis(float v)
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftZDis)
      FPrimaryVar(i) = v;
}

void MkNode::SetXAng(float v)
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftXAng)
      FPrimaryVar(i) = v;
}

void MkNode::SetYAng(float v)
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftYAng)
      FPrimaryVar(i) = v;
}

void MkNode::SetZAng(float v)
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftZAng)
      FPrimaryVar(i) = v;
}

float MkNode::GetXDis()
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftXDis)
      return FPrimaryVar(i);
  return 0;
}

float MkNode::GetYDis()
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftYDis)
      return FPrimaryVar(i);
  return 0;
}

float MkNode::GetZDis()
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftZDis)
      return FPrimaryVar(i);
  return 0;
}

float MkNode::GetXAng()
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftXAng)
      return FPrimaryVar(i);
  return 0;
}

float MkNode::GetYAng()
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftYAng)
      return FPrimaryVar(i);
  return 0;
}

float MkNode::GetZAng()
{
  for (int i=0;i<FDOFs.GetSize();i++)
    if(FDOFs[i].DOFType == doftZAng)
      return FPrimaryVar(i);
  return 0;
}

bool MkNode::isEq(MkNode &node)
{
  return FPoint == node.FPoint &&
    FNodeType == node.FNodeType &&
    FDOFs == node.FDOFs &&
    FPrimaryVar == node.FPrimaryVar;
}

bool MkNode::isNE(MkNode &node)
{
  return !isEq(node);
}

bool MkNode::operator==(MkNode &node)
{
  bool flag;
  flag = fabs(FPoint.X-node.FPoint.X) < EPS &&
         fabs(FPoint.Y-node.FPoint.Y) < EPS &&
         fabs(FPoint.Z-node.FPoint.Z) < EPS;
  return flag;
}

bool MkNode::operator!=(MkNode &node)
{
  return !(*this==node);
}

#ifdef __BCPLUSPLUS__
void MkNode::Draw(TObject *Sender)
{
  FPoint.Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkNode::Draw(MkPaint *pb)
{
  FPoint.Draw(pb);
}
#endif

void Swap(MkNode &n1, MkNode &n2)
{
  MkNode n;
  n = n1;
  n1 = n2;
  n2 = n;
}
//--------------------------------------------------------------------
MkNodes::MkNodes(int size,MkNode *nodes)
{

    if (size < 0) {
      MkDebug("::MkNodes - MkNodes(int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0) {
       FNode = NULL;
       return;
    }

    FNode = new MkNode[FSize];
    for (int i=0;i<FSize;i++) (*this)[i] = nodes[i];
}

MkNodes::MkNodes(int size)
{
    if (size < 0) {
      MkDebug("::MkNodes - MkNodes(int size)");;
      return;
    }

    FSize = FSizeOfArray = size;

    if (FSizeOfArray == 0) {
       FNode = NULL;
       return;
    }

    FNode = new MkNode[FSizeOfArray];
}

MkNodes::~MkNodes()
{
   FSizeOfArray = FSize = 0;
   if (FNode) {
      delete[] FNode;
      FNode = NULL;
   }
}

void MkNodes::Initialize(int size)
{
    if (size < 0) {
      MkDebug("::MkNodes - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;

    FSize = FSizeOfArray = size;

    if (FSizeOfArray == 0) {
       if (FNode!=NULL) delete[] (MkNode*)FNode;
       FNode = NULL;
       return;
    }

    if (FNode!=NULL) delete[] (MkNode*)FNode;
    FNode = new MkNode[FSizeOfArray];
}

void MkNodes::Initialize(int size,MkNode *nodes)
{

    if (size < 0 || nodes == NULL) {
      MkDebug("::MkNodes - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0) {
       if (FNode!=NULL) delete[] (MkNode*)FNode;
       FNode = NULL;
       return;
    }

    if (FNode!=NULL) delete[] (MkNode*)FNode;
    FNode = new MkNode[FSizeOfArray];
    for (int i=0;i<FSizeOfArray;i++) FNode[i] = nodes[i];
}

int MkNodes::Grow(int delta)
{
    int i;
    MkNode *node=NULL;

    if (!(node = new MkNode[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        node[i] = FNode[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        node[i] = NullNode;
    if (FNode) {
       delete[] (MkNode*)FNode;
       FNode = NULL;
    }
    FNode = node;
    FSizeOfArray = FSizeOfArray+delta;

    return FSizeOfArray;
}

int MkNodes::Shrink(int delta)
{
    int i;
    MkNode *node=NULL;

    if (!(node = new MkNode[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        node[i] = FNode[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        node[i] = NullNode;
    if (FNode) {
       delete[] (MkNode*)FNode;
       FNode = NULL;
    }
    FNode = node;
    FSizeOfArray = FSizeOfArray-delta;

    return FSizeOfArray;
}

bool MkNodes::Add(MkNode &node)
{
    int tmp=FSizeOfArray;
    bool flag=false;
    for (int i=0;i<FSize;i++) if (FNode[i]==node) flag=true;

    if(flag) return false;
    if (FSize>=FSizeOfArray) {
       Grow(FSize-FSizeOfArray+10);
       if (tmp==FSizeOfArray) return false;
    }
    FSize++;
    FNode[FSize-1] = node;

    return true;
}

bool MkNodes::Add(int index, MkNode &node)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) Grow(FSize-FSizeOfArray+1);
    if(tmp==FSizeOfArray) return false;

    for (int i=FSize-1;i>=index;i--)
      FNode[i+1] = FNode[i];
    FSize++;
    FNode[index] = node;
    return true;
}

bool MkNodes::Delete(MkNode &node)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FNode[i] == node) break;
    }
    if(i==FSize) return false;
    if(FNode[i] == node) {
      for (int j=i;j<FSize-1;j++)
        FNode[j] = FNode[j+1];
    }
    FSize--;
    FNode[FSize] = NullNode;
    return true;
}

bool MkNodes::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FNode[j] = FNode[j+1];

    FSize--;
    FNode[FSize] = NullNode;
    return true;
}

bool MkNodes::Clear()
{
   FSizeOfArray = FSize = 0;
   if (FNode) {
      delete[] FNode;
      FNode = NULL;
   }
   return true;
}

MkNode & MkNodes::operator[](int i)
{
    if (FSizeOfArray == 0) return NullNode;
    if (i >= FSize && i < FSizeOfArray) FSize = i+1;

    if (i >=0 && i < FSize) return FNode[i];
    else return NullNode;
}

MkNodes & MkNodes::operator=(MkNodes &nodes)
{
    int i;

    Clear();
    FSize = nodes.FSize;
    FSizeOfArray = nodes.FSizeOfArray;
    if (FSize == 0) {
       FNode = NULL;
       return *this;
    }
    this->FNode = new MkNode[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FNode[i] = nodes.FNode[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FNode[i] = NullNode;

    return *this;
}

bool MkNodes::operator==(MkNodes &nodes)
{
  int i;

  if (FSize != nodes.FSize) return false;
  for (i=0;i<FSize;i++)
    if (this->FNode[i] != nodes.FNode[i]) return false;

  return true;
}

#ifdef __BCPLUSPLUS__
void MkNodes::Draw(TObject *Sender)
{
  for(int i=0;i<FSize;i++)
    FNode[i].Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkNodes::Draw(MkPaint *pb)
{
  for(int i=0;i<FSize;i++)
    FNode[i].Draw(pb);
}
#endif
//--------------------------------------------------------------------
void MkElement::FindCenter()
{
  int i;
  Center.SetPoint(0,0,0);
  for (i=0;i<ElemNode.getSzX();i++) {
    Center = Center+GetElemNode(i).GetPoint();
  }
  Center = Center/(int)ElemNode.getSzX();
}

void MkElement::FindCentroid()
{

}

MkElement & MkElement::operator=(MkElement &elem)
{
  ElemNode = elem.ElemNode;
  Prop = elem.Prop;
  NodeRef = elem.NodeRef;
  Stiff = elem.Stiff;
  TranMat = elem.TranMat;
  Center=elem.Center;
  Centroid=elem.Centroid;
#ifdef __BCPLUSPLUS__
  Color=elem.Color;
#endif
  return *this;
}

bool MkElement::operator==(MkElement &elem)
{
  if (ElemNode.getSzX() !=elem.ElemNode.getSzX()) return false;
  for (int i=0;i<ElemNode.getSzX();i++)
    if ((*NodeRef)[ElemNode(i)] != (*elem.NodeRef)[elem.ElemNode(i)])
      return false;
  return true;
}

bool MkElement::operator!=(MkElement &elem)
{
  return !(*this==elem);
}

#ifdef __BCPLUSPLUS__
void MkElement::Out(TMemo *memo)
{
  int i;
  char str[256],s[256];
  MkNodes &node=*NodeRef;
  MkInt &elemnode=ElemNode;

  if(!memo) return;

  sprintf(str,"ElemNode:: \n");
  for(i=0;i<ElemNode.getSzX();i++) {
    sprintf(s,"elemnode(%d)=%d, ",i,ElemNode(i));
    strcat(str,s);
  }
  strcat(str,"\n");
  memo->Lines->Add(str);

  sprintf(str,"Node :: \n");
  for(i=0;i<ElemNode.getSzX();i++) {
    int n = elemnode(i);
    sprintf(s,"node(%d)=(%f,%f,%f), \n",n,node[n].GetPoint().X,
                                        node[n].GetPoint().Y,
                                        node[n].GetPoint().Z);
    strcat(str,s);
  }
  strcat(str,"\n");
  memo->Lines->Add(str);
}
#endif

void MkElement::Out(char *fname)
{
  int i;
  FILE *fp;

  if(!(fp=fopen(fname,"a"))) return;
  char str[256],s[256];
  MkNodes &node=*NodeRef;
  MkInt &elemnode=ElemNode;

  sprintf(str,"ElemNode:: \n");
  for(i=0;i<ElemNode.getSzX();i++) {
    sprintf(s,"elemnode(%d)=%d, ",i,ElemNode(i));
    strcat(str,s);
  }
  strcat(str,"\n");
  fputs(str,fp);

  sprintf(str,"Node :: \n");
  for(i=0;i<ElemNode.getSzX();i++) {
    int n = elemnode(i);
    sprintf(s,"node(%d)=(%f,%f,%f), \n",n,node[n].GetPoint().X,
                                        node[n].GetPoint().Y,
                                        node[n].GetPoint().Z);
    strcat(str,s);
  }
  strcat(str,"\n");
  fputs(str,fp);
  fclose(fp);
}

void MkElement::Out()
{
  int i;
  char str[256],s[256];
  MkNodes &node=*NodeRef;
  MkInt &elemnode=ElemNode;
  puts("MkElement's Out()\n");

#ifdef __BCPLUSPLUS__
  sprintf(str,"ElemNode:: (%s)\n",ClassName().c_str());
#else
  sprintf(str,"ElemNode:: (%s)\n",ClassName().c_str());
#endif
  for(i=0;i<ElemNode.getSzX();i++) {
    sprintf(s,"elemnode(%d)=%d, ",i,ElemNode(i));
    strcat(str,s);
  }
  strcat(str,"\n");
  puts(str);

  sprintf(str,"Node :: \n");
  for(i=0;i<ElemNode.getSzX();i++) {
    int n = elemnode(i);
    sprintf(s,"node(%d)=(%f,%f,%f), \n",n,node[n].GetPoint().X,
                                          node[n].GetPoint().Y,
                                          node[n].GetPoint().Z);
    strcat(str,s);
  }
  strcat(str,"\n");
  puts(str);
}

//--------------------------------------------------------------------
MkLineElement::MkLineElement(MkInt &elemnode, MkNodes &nodes)
{
  MkPoint sp,ep;
  ElemNode = elemnode;
  NodeRef = &nodes;
  assert(ElemNode.getSzX()==2);
  assert(NodeRef);
  sp = GetElemNode(0).GetPoint();
  ep = GetElemNode(1).GetPoint();
  Length = CalDist(sp,ep);
  Line.SetLine(sp,ep);
  className = "MkLineElement";
}

float MkLineElement::GetLength()
{
  assert(NodeRef);
  assert(ElemNode.getSzX()==2);
  MkPoint pnt1,pnt2;
  pnt1 = (*NodeRef)[ElemNode(0)].GetPoint();
  pnt2 = (*NodeRef)[ElemNode(1)].GetPoint();
  Length = CalDist(pnt1,pnt2);
  return Length;
}

bool MkLineElement::IsIn(MkNode &node)
{
  char str[256];
  assert(NodeRef);
  MkPoint p = node.GetPoint();
  MkLine l((*NodeRef)[ElemNode(0)].GetPoint(),(*NodeRef)[ElemNode(1)].GetPoint());
  l.SetFiniteness(true);
  bool flag = l.IsInLine(p);
  if(flag) {
    sprintf(str,"p(%f,%f,%f) is in sp(%f,%f,%f) - ep(%f,%f,%f)\n",p.X,p.Y,p.Z,l[0].X,l[0].Y,l[0].Z,l[1].X,l[1].Y,l[1].Z);
    //MkDebug(str);
  }
  return flag;
}

bool MkLineElement::IsInside(int node)  // is node located in this element
{
  return false;
}

bool MkLineElement::IsIncluded(int node)
{
  assert(ElemNode.getSzX()==2);
  for (int i=0;i<ElemNode.getSzX();i++)
    if(ElemNode(i)==node) return true;
  return false;
}

bool MkLineElement::IsCross(MkElement *elem)
{
  MkLineElement *le = dynamic_cast<MkLineElement*>(elem);
  if(!le) return false;

  MkLine line[2];
  MkPoint p[2];
  p[0] = (*NodeRef)[ElemNode(0)].GetPoint();
  p[1] = (*NodeRef)[ElemNode(1)].GetPoint();
  line[0].SetLine(p[0],p[1]);
  p[0] = (*le->NodeRef)[ElemNode(0)].GetPoint();
  p[1] = (*le->NodeRef)[ElemNode(1)].GetPoint();
  line[1].SetLine(p[0],p[1]);

  return line[0]&&line[1]; 
}

MkCrossType MkLineElement::GetCrossType(MkElement *elem)
{
  bool flag;
  MkLineElement *le = dynamic_cast<MkLineElement*>(elem);
  if(!le) return ctNot;

  MkLine line[2];
  MkPoint p[4];
  p[0] = (*NodeRef)[ElemNode(0)].GetPoint();
  p[1] = (*NodeRef)[ElemNode(1)].GetPoint();
  line[0].SetLine(p[0],p[1]);
  line[0].SetFiniteness(true);

  p[2] = (*le->NodeRef)[ElemNode(0)].GetPoint();
  p[3] = (*le->NodeRef)[ElemNode(1)].GetPoint();
  line[1].SetLine(p[2],p[3]);
  line[1].SetFiniteness(true);

  flag = line[0]&&line[1];

  if(!flag) {
    return ctNot;
  }
  else if(p[0]==p[2] || p[0]==p[3] || p[1]==p[2] || p[1]==p[3]) {
    return ctNodeToNode;
  }
  else if(line[0].IsInLine(p[2]) || line[0].IsInLine(p[3])) {
    return ctEdgeToNode;
  }
  else if(line[1].IsInLine(p[0]) || line[1].IsInLine(p[1])) {
    return ctNodeToEdge;
  }
  else {
    return ctIntersect;
  }
}

bool MkLineElement::SetupTrans()
{
  int i,j,h;
  int cnt,size;

  MkInt ti;
  assert(NodeRef);
  assert(ElemNode.getSzX()==2);
  int ndof=0;
  for (i=0;i<ElemNode.getSzX();i++)
      ndof += (*NodeRef)[ElemNode(i)].GetDOFs().GetSize();

  cnt=0;
  ti.Initialize(ndof);
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = (*NodeRef)[ElemNode(i)].GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      switch(dofs[j].DOFType) {
        case doftXDis:
        case doftXAng: ti[cnt] = 0;break;
        case doftYDis:
        case doftYAng: ti[cnt] = 1;break;
        case doftZDis:
        case doftZAng: ti[cnt] = 2;break;
        default: ti[cnt] = 3;break;
      }
      cnt++;
    }
  }

  assert(ndof==cnt);

  TranMat.Initialize(ndof,ndof);
  TranMat.LoadIdentity();

  MkMatrix4 Tran;
  Tran.Identity();
  MkPoint pnt1,pnt2,pnt3;
  pnt1 = (*NodeRef)[ElemNode(0)].GetPoint();
  pnt2 = (*NodeRef)[ElemNode(1)].GetPoint();

  float k,l,m,n;
  l = pnt2.X-pnt1.X;
  m = pnt2.Y-pnt1.Y;
  n = pnt2.Z-pnt1.Z;

  Length = sqrt(l*l+m*m+n*n);
  assert(fabs(Length)>0.001);
  l /= Length; m /= Length; n /= Length;
  k = sqrt(l*l+m*m);

  float theta, phi;

  if(fabs(k)>0.0001) {
    theta = acos(l/k);
    if(fabs(m-k*sin(theta)) > 0.0001) theta = -theta;
    if(theta<0) theta+= 2*(float)M_PI;
    if(theta>2*M_PI) theta-= 2*(float)M_PI;
    assert(fabs(m-k*sin(theta)) < (float)0.0001);
    phi = acos(n);
  }
  else {
    theta = 0;
    phi = 0;
  }

  Tran.RotateInY(phi*180/M_PI-90);
  Tran.RotateInZ(theta*180/M_PI);

  cnt=0;
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = (*NodeRef)[ElemNode(i)].GetDOFs();
    size = dofs.GetSize();
    for (j=0;j<size;j++) {
      for (h=j;h<size;h++) {
        if(dofs[h].isSameKindWith(dofs[j])) {
          TranMat(cnt+h,cnt+j) = Tran(ti(cnt+h),ti(cnt+j));
          TranMat(cnt+j,cnt+h) = Tran(ti(cnt+h),ti(cnt+j));
        }
      }
    }
    cnt+=size;
  }

  assert(cnt==ndof);
  return true;
}

MkLineElement & MkLineElement::operator=(MkLineElement &elem)
{
  MkElement::operator=((MkElement&)elem);
  Length = elem.Length;
  return *this;
}

#ifdef __BCPLUSPLUS__
void MkLineElement::Draw(TObject *Sender)
{
  TColor C=clBlack; 
  MkPoint p[2];
  p[0] = GetElemNode(0).GetPoint();
  p[1] = GetElemNode(1).GetPoint();
  MkLine l(p[0],p[1],C);

  l.Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkLineElement::Draw(MkPaint *pb)
{
  MkPoint p[2];
  p[0] = GetElemNode(0).GetPoint();
  p[1] = GetElemNode(1).GetPoint();
  MkLine l(p[0],p[1]);

  l.Draw(pb);
}
#endif
//--------------------------------------------------------------------
bool MkTrussElement::Post()
{
  return false;
}

bool MkTrussElement::SetupStiff()
{
  if(isStiffed) return true;
  int cnt;
  int ndof=0;
  int i,j;
  MkInt ti;
  MkMatrix km(12,12);
  MkMatrix t;

  assert(NodeRef);
  assert(ElemNode.getSzX()==2);

  if(TranMat.GetFI()==0) SetupTrans();

  assert(TranMat.GetFI()*TranMat.GetFJ());
  assert(TranMat.GetFI()==TranMat.GetFJ());

  for (i=0;i<ElemNode.getSzX();i++) {
    ndof += (*NodeRef)[ElemNode(i)].GetDOFs().GetSize();
  }

  assert(TranMat.GetFI()==ndof);
  Stiff.Initialize(ndof,ndof);

  cnt=0;
  ti.Initialize(ndof);
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = (*NodeRef)[ElemNode(i)].GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      switch(dofs[j].DOFType) {
        case doftXDis:ti[cnt] = 0+i*6; break;
        case doftYDis:ti[cnt] = 1+i*6; break;
        case doftZDis:ti[cnt] = 2+i*6; break;
        case doftXAng:ti[cnt] = 3+i*6; break;
        case doftYAng:ti[cnt] = 4+i*6; break;
        case doftZAng:ti[cnt] = 5+i*6; break;
        default:      ti[cnt] = -1;    break;
      }
      cnt++;
    }
  }

  float ea = Prop(0)*Prop(1);
  float a = ea/Length;

  km(0,0)= a;
  km(0,6)=-a;
  km(6,0)=-a;
  km(6,6)= a;

  t = TranMat;
  t.Transpose();
  t*= km;
  t*= TranMat;
  km = t;

  for (i=0;i<ndof;i++)
    for (j=i;j<ndof;j++) {
      Stiff(i,j) = km(ti(i),ti(j));
      Stiff(j,i) = km(ti(i),ti(j));
    }
  SetStiffed(true);
  return true;
}

float MkTrussElement::GetAxialForce(float xi)
{
  float ea=Prop(0)*Prop(1);
  return GetAxialStrain(xi)*ea+JackingForce+IniForceCorrection;
}

float MkTrussElement::GetAxialStress(float xi)
{
  float e = Prop(0),a = Prop(1);
  return GetAxialStrain(xi)*e+(fabs(a)>EPS?1/a:0)*JackingForce+(fabs(a)>EPS?1/a:0)*IniForceCorrection;
}

float MkTrussElement::GetAxialStrain(float xi)
{
  int i,j,ndof=0,cnt,cnt2;
  float u[2];
  MkVector ul,ug;
  MkMatrix t;

  ndof += GetElemNode(0).GetDOFs().GetSize();
  ndof += GetElemNode(1).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j) = GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j) = GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftXDis) {
	u[0] = ul(j);
          break;
      }
    }
    i = dofs.GetSize();
  } while(false);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftXDis) {
	u[1] = ul(i+j);
          break;
      }
    }
  } while(false);

  return dN1(xi)*u[0]+dN2(xi)*u[1];
}

float MkTrussElement::GetDisp(int s)
{
  if(s!=0 && s!=1) return 0;

  int i,j,ndof=0,cnt,cnt2;
  float u[2];
  MkVector ul,ug;
  MkMatrix t;

  ndof += GetElemNode(0).GetDOFs().GetSize();
  ndof += GetElemNode(1).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j) = GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j) = GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftXDis) {
	u[0] = ul(j);
          break;
      }
    }
    i = dofs.GetSize();
  } while(false);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftXDis) {
	u[1] = ul(i+j);
          break;
      }
    }
  } while(false);

  return u[s];
}

float MkTrussElement::N1(float xi)
{
  return 1-xi;
}

float MkTrussElement::N2(float xi)
{
  return xi;
}

float MkTrussElement::dN1(float xi)
{
  return Length<EPS?0:-1/Length;
}

float MkTrussElement::dN2(float xi)
{
  return Length<EPS?0:1/Length;
}

MkTrussElement & MkTrussElement::operator=(MkTrussElement &elem)
{
  MkLineElement::operator=(elem);
  JackingForce = elem.JackingForce;
  return *this;
}

#ifdef __BCPLUSPLUS__
void MkTrussElement::Out(TMemo *memo)
{
  if(!memo) return;
  memo->Lines->Add("MkTrussElement stiffness Output");
  Stiff.Out(memo);
}
#endif

void MkTrussElement::Out(char *fname)
{
  FILE *fp;
  fp = fopen(fname,"a");
  if(!strlen(fname)) return;
  fprintf(fp,"MkTrussElement stiffness Output");
  fclose(fp);

  Stiff.Out(fname);
}

void MkTrussElement::Out()
{
  int i;
  char str[256],s[256];
  MkNodes &node=*NodeRef;
  MkInt &elemnode=ElemNode;

  puts("MkTrussElement's Out()\n");

#ifdef __BCPLUSPLUS__
  sprintf(str,"ElemNode:: (%s)\n",ClassName().c_str());
#else
  sprintf(str,"ElemNode:: (%s)\n",ClassName().c_str());
#endif
  for(i=0;i<ElemNode.getSzX();i++) {
    sprintf(s,"elemnode(%d)=%d, ",i,ElemNode(i));
    strcat(str,s);
  }
  strcat(str,"\n");
  puts(str);

  sprintf(str,"Node :: \n");
  for(i=0;i<ElemNode.getSzX();i++) {
    int n = elemnode(i);
    sprintf(s,"node(%d)=(%f,%f,%f), \n",n,node[n].GetPoint().X,
                                          node[n].GetPoint().Y,
                                          node[n].GetPoint().Z);
    strcat(str,s);
  }
  strcat(str,"\n");
  puts(str);
}

#ifdef __BCPLUSPLUS__
void MkTrussElement::Draw(TObject *Sender)
{
  MkLineElement::Draw(Sender);
}

void MkTrussElement::DrawResult(TObject *Sender,MkResultType rt)
{

}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkTrussElement::Draw(MkPaint *pb)
{
  MkLineElement::Draw(pb);
}

void MkTrussElement::DrawResult(MkPaint *pb,MkResultType rt)
{

}
#endif
//--------------------------------------------------------------------
bool MkBeamElement::Post()
{
  int cnt;
  int ndof=0;
  int i,j;
  MkInt ti;
  MkVector ul,ug,fixend;
  MkMatrix km(12,12);
  MkMatrix t;

  assert(NodeRef);
  assert(ElemNode.getSzX()==2);

  for (i=0;i<ElemNode.getSzX();i++) {
    ndof += (*NodeRef)[ElemNode(i)].GetDOFs().GetSize();
  }

  assert(TranMat.GetFI()==ndof);

  cnt=0;
  ti.Initialize(ndof);
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = (*NodeRef)[ElemNode(i)].GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      switch(dofs[j].DOFType) {
        case doftXDis: ti[cnt] = 0+i*6; break;
        case doftYDis: ti[cnt] = 1+i*6; break;
        case doftZDis: ti[cnt] = 2+i*6; break;
        case doftXAng:  ti[cnt] = 3+i*6; break;
        case doftYAng:  ti[cnt] = 4+i*6; break;
        case doftZAng:  ti[cnt] = 5+i*6; break;
        default:        ti[cnt] = -1;    break;
      }
      cnt++;
    }
  }

  float ea = Prop(0)*Prop(1);
  float eiy = Prop(0)*Prop(3);
  float eiz = Prop(0)*Prop(4);
  float gj = Prop(2);

  float a1,a2,a3,a4,a5,a6,a7,a8;

  a1=ea/Length;
  a2=12.*eiz/(Length*Length*Length);
  a3=12.*eiy/(Length*Length*Length);
  a4=6.*eiz/(Length*Length);
  a5=6.*eiy/(Length*Length);
  a6=4.*eiz/Length;
  a7=4.*eiy/Length;
  a8=gj/Length;

  km(0,0)=km(6,6)=a1;
  km(0,6)=km(6,0)=-a1;
  km(1,1)=km(7,7)=a2;
  km(1,7)=km(7,1)=-a2;
  km(2,2)=km(8,8)=a3;
  km(2,8)=km(8,2)=-a3;
  km(3,3)=km(9,9)=a8;
  km(3,9)=km(9,3)=-a8;
  km(4,4)=km(10,10)=a7;
  km(4,10)=km(10,4)=.5*a7;
  km(5,5)=km(11,11)=a6;
  km(5,11)=km(11,5)=.5*a6;
  km(1,5)=km(5,1)=a4;
  km(1,11)=km(11,1)=a4;
  km(5,7)=km(7,5)=-a4;
  km(7,11)=km(11,7)=-a4;
  km(4,8)=km(8,4)=a5;
  km(8,10)=km(10,8)=a5;
  km(2,4)=km(4,2)=-a5;
  km(2,10)=km(10,2)=-a5;

  km(1,1) += (AxialLoad/30)*(36/Length);
  km(2,2) += (AxialLoad/30)*(36/Length);
  km(4,4) += (AxialLoad/30)*(4*Length);
  km(5,5) += (AxialLoad/30)*(4*Length);
  km(7,7) += (AxialLoad/30)*(36/Length);
  km(8,8) += (AxialLoad/30)*(36/Length);
  km(10,10) += (AxialLoad/30)*(4*Length);
  km(11,11) += (AxialLoad/30)*(4*Length);

  km(1,5) += (AxialLoad/30)*(3);
  km(5,1) += (AxialLoad/30)*(3);
  km(1,7) += (AxialLoad/30)*(-36/Length);
  km(7,1) += (AxialLoad/30)*(-36/Length);
  km(1,11) += (AxialLoad/30)*(3);
  km(11,1) += (AxialLoad/30)*(3);

  km(2,4) += (AxialLoad/30)*(-3);
  km(4,2) += (AxialLoad/30)*(-3);
  km(2,8) += (AxialLoad/30)*(-36/Length);
  km(8,2) += (AxialLoad/30)*(-36/Length);
  km(2,10) += (AxialLoad/30)*(-3);
  km(10,2) += (AxialLoad/30)*(-3);

  km(4,8) += (AxialLoad/30)*(3);
  km(8,4) += (AxialLoad/30)*(3);
  km(4,10) += (AxialLoad/30)*(-Length);
  km(10,4) += (AxialLoad/30)*(-Length);

  km(5,7) += (AxialLoad/30)*(-3);
  km(7,5) += (AxialLoad/30)*(-3);
  km(5,11) += (AxialLoad/30)*(-Length);
  km(11,5) += (AxialLoad/30)*(-Length);

  km(7,11) += (AxialLoad/30)*(-3);
  km(11,7) += (AxialLoad/30)*(-3);
  km(8,10) += (AxialLoad/30)*(3);
  km(10,8) += (AxialLoad/30)*(3);

  ul.Initialize(ndof);
  ug.Initialize(ndof);
  fixend.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j)= GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j)= GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;
  t = TranMat;
  t.Transpose();
  fixend=t*FixedEnd;

  Resultant = km*ul + fixend;

  return true;
}

bool MkBeamElement::SetupStiff()
{
  if(isStiffed) return true;
  int cnt;
  int ndof=0;
  int i,j;
  MkInt ti;
  MkMatrix km(12,12);
  MkMatrix t;

  assert(NodeRef);
  assert(ElemNode.getSzX()==2);

  if(TranMat.GetFI()==0) SetupTrans();

  assert(TranMat.GetFI()*TranMat.GetFJ());
  assert(TranMat.GetFI()==TranMat.GetFJ());

  for (i=0;i<ElemNode.getSzX();i++) {
    ndof += (*NodeRef)[ElemNode(i)].GetDOFs().GetSize();
  }

  assert(TranMat.GetFI()==ndof);
  Stiff.Initialize(ndof,ndof);

  cnt=0;
  ti.Initialize(ndof);
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = (*NodeRef)[ElemNode(i)].GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      switch(dofs[j].DOFType) {
        case doftXDis: ti[cnt] = 0+i*6; break;
        case doftYDis: ti[cnt] = 1+i*6; break;
        case doftZDis: ti[cnt] = 2+i*6; break;
        case doftXAng:  ti[cnt] = 3+i*6; break;
        case doftYAng:  ti[cnt] = 4+i*6; break;
        case doftZAng:  ti[cnt] = 5+i*6; break;
        default:        ti[cnt] = -1;    break;
      }
      cnt++;
    }
  }

  float ea = Prop(0)*Prop(1);
  float eiy = Prop(0)*Prop(3);
  float eiz = Prop(0)*Prop(4);
  float gj = Prop(2);

  float a1,a2,a3,a4,a5,a6,a7,a8;

  a1=ea/Length;
  a2=12.*eiz/(Length*Length*Length);
  a3=12.*eiy/(Length*Length*Length);
  a4=6.*eiz/(Length*Length);
  a5=6.*eiy/(Length*Length);
  a6=4.*eiz/Length;
  a7=4.*eiy/Length;
  a8=gj/Length;

  km(0,0)=km(6,6)=a1;
  km(0,6)=km(6,0)=-a1;
  km(1,1)=km(7,7)=a2;
  km(1,7)=km(7,1)=-a2;
  km(2,2)=km(8,8)=a3;
  km(2,8)=km(8,2)=-a3;
  km(3,3)=km(9,9)=a8;
  km(3,9)=km(9,3)=-a8;
  km(4,4)=km(10,10)=a7;
  km(4,10)=km(10,4)=.5*a7;
  km(5,5)=km(11,11)=a6;
  km(5,11)=km(11,5)=.5*a6;
  km(1,5)=km(5,1)=a4;
  km(1,11)=km(11,1)=a4;
  km(5,7)=km(7,5)=-a4;
  km(7,11)=km(11,7)=-a4;
  km(4,8)=km(8,4)=a5;
  km(8,10)=km(10,8)=a5;
  km(2,4)=km(4,2)=-a5;
  km(2,10)=km(10,2)=-a5;

  km(1,1) += (AxialLoad/30)*(36/Length);
  km(2,2) += (AxialLoad/30)*(36/Length);
  km(4,4) += (AxialLoad/30)*(4*Length);
  km(5,5) += (AxialLoad/30)*(4*Length);
  km(7,7) += (AxialLoad/30)*(36/Length);
  km(8,8) += (AxialLoad/30)*(36/Length);
  km(10,10) += (AxialLoad/30)*(4*Length);
  km(11,11) += (AxialLoad/30)*(4*Length);

  km(1,5) += (AxialLoad/30)*(3);
  km(5,1) += (AxialLoad/30)*(3);
  km(1,7) += (AxialLoad/30)*(-36/Length);
  km(7,1) += (AxialLoad/30)*(-36/Length);
  km(1,11) += (AxialLoad/30)*(3);
  km(11,1) += (AxialLoad/30)*(3);

  km(2,4) += (AxialLoad/30)*(-3);
  km(4,2) += (AxialLoad/30)*(-3);
  km(2,8) += (AxialLoad/30)*(-36/Length);
  km(8,2) += (AxialLoad/30)*(-36/Length);
  km(2,10) += (AxialLoad/30)*(-3);
  km(10,2) += (AxialLoad/30)*(-3);

  km(4,8) += (AxialLoad/30)*(3);
  km(8,4) += (AxialLoad/30)*(3);
  km(4,10) += (AxialLoad/30)*(-Length);
  km(10,4) += (AxialLoad/30)*(-Length);

  km(5,7) += (AxialLoad/30)*(-3);
  km(7,5) += (AxialLoad/30)*(-3);
  km(5,11) += (AxialLoad/30)*(-Length);
  km(11,5) += (AxialLoad/30)*(-Length);

  km(7,11) += (AxialLoad/30)*(-3);
  km(11,7) += (AxialLoad/30)*(-3);
  km(8,10) += (AxialLoad/30)*(3);
  km(10,8) += (AxialLoad/30)*(3);

  t = TranMat;
  t.Transpose();
  t*= km;
  t*= TranMat;
  km = t;

  for (i=0;i<ndof;i++)
    for (j=i;j<ndof;j++) {
      Stiff(i,j) = km(ti(i),ti(j));
      Stiff(j,i) = km(ti(i),ti(j));
    }
  SetStiffed(true);
  return true;
}

float MkBeamElement::N1(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return 1-3*xi*xi+2*xi*xi*xi;
}

float MkBeamElement::N2(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return Length*(-xi+2*xi*xi-xi*xi*xi);
}

float MkBeamElement::N3(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return 3*xi*xi-2*xi*xi*xi;
}

float MkBeamElement::N4(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return Length*(xi*xi-xi*xi*xi);
}

float MkBeamElement::dN1(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return (-6*xi+6*xi*xi)/Length;
}

float MkBeamElement::dN2(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return 4*xi-3*xi*xi;
}

float MkBeamElement::dN3(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return (6*xi-6*xi*xi)/Length;
}

float MkBeamElement::dN4(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return 2*xi-3*xi*xi;
}

float MkBeamElement::d2N1(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return (-6+12*xi)/Length/Length;
}

float MkBeamElement::d2N2(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return (4-6*xi)/Length;
}

float MkBeamElement::d2N3(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return (6-12*xi)/Length/Length;
}

float MkBeamElement::d2N4(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return (2-6*xi)/Length;
}

float MkBeamElement::d3N1(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return 12/Length/Length/Length;
}

float MkBeamElement::d3N2(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return -6/Length/Length;
}

float MkBeamElement::d3N3(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return -12/Length/Length/Length;
}

float MkBeamElement::d3N4(float xi)
{
  if (xi<0.0-EPS||xi>1.0+EPS) return 0;
  return -6/Length/Length;
}

float MkBeamElement::GetDisp(float xi)
{
  int i,j,ndof=0,cnt,cnt2;
  MkVector u(4);
  MkVector d(4);
  MkVector ul,ug;
  MkMatrix t;

  for (i=0;i<ElemNode.getSzX();i++)
    ndof += GetElemNode(i).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j)= GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j)= GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  cnt=0;
  cnt2=0;
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = GetElemNode(i).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftYDis) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      else if(dofs[j].DOFType == doftZAng) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      cnt++;
      if (cnt>ndof) break;
    }
  }

  d[0] = N1(xi);
  d[1] = N2(xi);
  d[2] = N3(xi);
  d[3] = N4(xi);

  d*=Prop(0)*Prop(4);

  return d.Dot(u);
}

float MkBeamElement::GetAngDis(float xi)
{
  int i,j,ndof=0,cnt,cnt2;
  MkVector u(4);
  MkVector d(4);
  MkVector ul,ug;
  MkMatrix t;

  for (i=0;i<ElemNode.getSzX();i++)
    ndof += GetElemNode(i).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j)= GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j)= GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  cnt=0;
  cnt2=0;
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = GetElemNode(i).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftYDis) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      else if(dofs[j].DOFType == doftZAng) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      cnt++;
      if (cnt>ndof) break;
    }
  }

  d[0] = dN1(xi);
  d[1] = dN2(xi);
  d[2] = dN3(xi);
  d[3] = dN4(xi);

  d*=Prop(0)*Prop(4);

  return d.Dot(u);
}

float MkBeamElement::GetMoment(float xi)
{
  int i,j,ndof=0,cnt,cnt2;
  MkVector u(4);
  MkVector d(4);
  MkVector ul,ug;
  MkMatrix t;

  for (i=0;i<ElemNode.getSzX();i++)
    ndof += GetElemNode(i).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j)= GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j)= GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  cnt=0;
  cnt2=0;
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = GetElemNode(i).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftYDis) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      else if(dofs[j].DOFType == doftZAng) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      cnt++;
      if (cnt>ndof) break;
    }
  }

  d[0] = d2N1(xi);
  d[1] = d2N2(xi);
  d[2] = d2N3(xi);
  d[3] = d2N4(xi);

  d*=Prop(0)*Prop(4);

  return d.Dot(u);
}

float MkBeamElement::GetShearForce(float xi)
{
  int i,j,ndof=0,cnt,cnt2;
  MkVector u(4);
  MkVector d(4);
  MkVector ul,ug;
  MkMatrix t;

  for (i=0;i<ElemNode.getSzX();i++)
    ndof += GetElemNode(i).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j)= GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j)= GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  cnt=0;
  cnt2=0;
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = GetElemNode(i).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftYDis) {
	if(cnt2>=4) MkDebug("MkTrussElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      else if(dofs[j].DOFType == doftZAng) {
	if(cnt2>=4) MkDebug("MkTrussElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      cnt++;
      if (cnt>ndof) break;
    }
  }

  d[0] = d3N1(xi);
  d[1] = d3N2(xi);
  d[2] = d3N3(xi);
  d[3] = d3N4(xi);

  d*=Prop(0)*Prop(4);

  return d.Dot(u);
}

float MkBeamElement::GetCompStress(float xi)
{
  float area,iz;
  area = Prop(1);
  iz = Prop(4);
  if(area<EPS || iz<EPS) return 0;
  return GetMoment(xi)*Height/(2*iz)+AxialLoad/area;
}

float MkBeamElement::GetShearStress(float xi)
{
  float iz,ht,h;
  iz = Prop(4);
  ht = Height;
  h = Height-Thickness;
  if(iz<EPS) return 0;
  return GetShearForce(xi)/(iz*Thickness)*(Width/2*(ht*ht/4-h*h/4)+Thickness/2*(h*h/4));
}

float MkBeamElement::GetTensStress(float xi)
{
  float area,iz;
  area = Prop(1);
  iz = Prop(4);
  if(area<EPS || iz<EPS) return 0;
  return -GetMoment(xi)*Height/(2*iz)+AxialLoad/area;
}

float MkBeamElement::GetDisp(int s)
{
  if(s!=0 && s!=1) {
    MkDebug("MkBeamElement::GetDisp(int s) takes s param only 1 and 0\n");
    return 0;
  }
  int i,j,ndof=0,cnt,cnt2;
  MkVector u(4);
  MkVector d(4);
  MkVector ul,ug;
  MkMatrix t;

  for (i=0;i<ElemNode.getSzX();i++)
    ndof += GetElemNode(i).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j)= GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j)= GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  cnt=0;
  cnt2=0;
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = GetElemNode(i).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftYDis) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      else if(dofs[j].DOFType == doftZAng) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      cnt++;
      if (cnt>ndof) break;
    }
  }

  return u[s*2];
}

float MkBeamElement::GetAngDisp(int s)
{
  if(s!=0 && s!=1) {
    MkDebug("MkBeamElement::GetAngDisp(int s) takes s param only 1 and 0\n");
    return 0;
  }
  int i,j,ndof=0,cnt,cnt2;
  MkVector u(4);
  MkVector d(4);
  MkVector ul,ug;
  MkMatrix t;

  for (i=0;i<ElemNode.getSzX();i++)
    ndof += GetElemNode(i).GetDOFs().GetSize();

  ul.Initialize(ndof);
  ug.Initialize(ndof);

  do {
    MkDOFs &dofs = GetElemNode(0).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(j)= GetElemNode(0)[j];
    }
    i=dofs.GetSize();
  }while(false);

  do {
    MkDOFs &dofs = GetElemNode(1).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      ug(i+j)= GetElemNode(1)[j];
    }
  }while(false);

  t = TranMat;
  t.Transpose();
  ul = t*ug;

  cnt=0;
  cnt2=0;
  for (i=0;i<ElemNode.getSzX();i++) {
    MkDOFs &dofs = GetElemNode(i).GetDOFs();
    for (j=0;j<dofs.GetSize();j++) {
      if(dofs[j].DOFType == doftYDis) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      else if(dofs[j].DOFType == doftZAng) {
	if(cnt2>=4) MkDebug("MkBeamElement::GetStrain  cnt2 exceeds.\n");
	u[cnt2] = ul(cnt);
	cnt2++;
      }
      cnt++;
      if (cnt>ndof) break;
    }
  }

  return u[s*2+1];
}


float MkBeamElement::GetMoment(int s) //s:0 first, s:1 last
{
  if(s!=0 && s!=1) {
    MkDebug("MkBeamElement::GetMoment(int s) takes s param only 1 and 0\n");
    return 0;
  }
  return Resultant[s*6+5];
}

float MkBeamElement::GetShearForce(int s) //s:0 first, s:1 last
{
  if(s!=0 && s!=1) {
    MkDebug("MkBeamElement::GetMoment(int s) takes s param only 1 and 0\n");
    return 0;
  }
  return Resultant[s*6+1];
}

MkBeamElement & MkBeamElement::operator=(MkBeamElement &elem)
{
  MkLineElement::operator=(elem);
  AxialLoad = elem.AxialLoad;
  Height = elem.Height;
  Width = elem.Width;
  Thickness = elem.Thickness;
  return *this;
}

#ifdef __BCPLUSPLUS__
void MkBeamElement::Out(TMemo *memo)
{
  if(!memo) return;
  memo->Lines->Add("MkBeamElement stiffness Output");
  Stiff.Out(memo);
}
#endif

void MkBeamElement::Out(char *fname)
{
  FILE *fp;
  fp = fopen(fname,"a");
  if(!strlen(fname)) return;
  fprintf(fp,"MkBeamElement stiffness Output");
  fclose(fp);

  Stiff.Out(fname);
}

void MkBeamElement::Out()
{
  int i;
  char str[256],s[256];
  MkNodes &node=*NodeRef;
  MkInt &elemnode=ElemNode;

  puts("MkBeamElement's Out()\n");

#ifdef __BCPLUSPLUS__
  sprintf(str,"ElemNode:: (%s)\n",ClassName().c_str());
#else
  sprintf(str,"ElemNode:: (%s)\n",ClassName().c_str());
#endif
  for(i=0;i<ElemNode.getSzX();i++) {
    sprintf(s,"elemnode(%d)=%d, ",i,ElemNode(i));
    strcat(str,s);
  }
  strcat(str,"\n");
  puts(str);

  sprintf(str,"Node :: \n");
  for(i=0;i<ElemNode.getSzX();i++) {
    int n = elemnode(i);
    sprintf(s,"node(%d)=(%f,%f,%f), \n",n,node[n].GetPoint().X,
                                          node[n].GetPoint().Y,
                                          node[n].GetPoint().Z);
    strcat(str,s);
  }
  strcat(str,"\n");
  puts(str);

  sprintf(str,"Axialload is %f\n",AxialLoad);
  puts(str);
}

#ifdef __BCPLUSPLUS__
void MkBeamElement::Draw(TObject *Sender)
{
  MkLineElement::Draw(Sender);
  if (haveResult) {
    DrawResult(Sender,rtCompFrc);
    DrawResult(Sender,rtCompStr);
  }
}

void MkBeamElement::DrawResult(TObject *Sender,MkResultType rt)
{
  int i,nseg = 9;
  float v[10]={0,0,0,0,0,0,0,0,0,0};
  float vmax=-1e20,vmin=1e20;
  MkPoint sp,ep,p;
  MkLine l,line;

  switch(rt) {
  case rtMoment:
    for(i=0;i<nseg;i++) v[i] = GetMoment(float((i-4)/4.0));
    break;
  case rtShearFrc:
    for(i=0;i<nseg;i++) v[i] = GetShearForce(float((i-4)/4.0));
    break;
  case rtCompStr:
    for(i=0;i<nseg;i++) v[i] = GetCompStress((i-4)/4.0);
    break;
  case rtShearStr:
    for(i=0;i<nseg;i++) v[i] = GetShearStress((i-4)/4.0);
    break;
  case rtTensStr:
    for(i=0;i<nseg;i++) v[i] = GetTensStress((i-4)/4.0);
    break;
  }
  for (i=0;i<nseg;i++) {
    vmax = vmax<v[i] ? v[i]:vmax;
    vmin = vmin>v[i] ? v[i]:vmin;
  }
  sp = GetElemNode(0).GetPoint();
  ep = GetElemNode(1).GetPoint();
  l.SetLine(sp,ep);

  for (i=0;i<nseg-1;i++) {
    sp = l.GetDivision(float(i)/nseg);
    sp.X += v[i];//*scale;
    ep = l.GetDivision(float(i+1)/nseg);
    ep.X += v[i+1];//*scale;

    line.SetLine(sp,ep);
    line.Draw(Sender);
  }

}
#endif


#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkBeamElement::Draw(MkPaint *pb)
{
  MkLineElement::Draw(pb);
  if (haveResult) {
    DrawResult(pb,rtCompFrc);
    DrawResult(pb,rtCompStr);
  }
}

void MkBeamElement::DrawResult(MkPaint *pb,MkResultType rt)
{
  int i,nseg = 9;
  float v[10]={0,0,0,0,0,0,0,0,0,0};
  float vmax=-(float)1e10,vmin=(float)1e10;
//  float scale;
  MkPoint sp,ep,p;
  MkLine l,line;

  switch(rt) {
  case rtMoment:
    for(i=0;i<nseg;i++) v[i] = GetMoment(float((i-4)/4.0));
    break;
  case rtShearFrc:
    for(i=0;i<nseg;i++) v[i] = GetShearForce(float((i-4)/4.0));
    break;
  case rtCompStr:
    for(i=0;i<nseg;i++) v[i] = GetCompStress(float((i-4)/4.0));
    break;
  case rtShearStr:
    for(i=0;i<nseg;i++) v[i] = GetShearStress(float((i-4)/4.0));
    break;
  case rtTensStr:
    for(i=0;i<nseg;i++) v[i] = GetTensStress(float((i-4)/4.0));
    break;
  }
  for (i=0;i<nseg;i++) {
    vmax = vmax<v[i] ? v[i]:vmax;
    vmin = vmin>v[i] ? v[i]:vmin;
  }
  sp = GetElemNode(0).GetPoint();
  ep = GetElemNode(1).GetPoint();
  l.SetLine(sp,ep);

  for (i=0;i<nseg-1;i++) {
    sp = l.GetDivision(float(i)/nseg);
    sp.X += v[i];//*scale;
    ep = l.GetDivision(float(i+1)/nseg);
    ep.X += v[i+1];//*scale;

    line.SetLine(sp,ep);
    line.Draw(pb);
  }
}
#endif
//--------------------------------------------------------------------
MkTriElement::MkTriElement()
{
  Radius = 0;
  className = "MkTriElement";
}

//--------------------------------------------------------------------
MkElements::MkElements(int size,MkElement **elems)
{
  int i;
  if (size < 0) {
    MkDebug("::MkElements - MkElements(int size)");
    return;
  }

  FSizeOfArray = FSize = size;
  if (FSize == 0) {
     FElement = NULL;
     return;
  }

  FElement = new MkElement*[FSize];
  assert(FElement);

  for (i=0;i<FSize;i++)
    dyn_cp(FElement[i],elems[i]);
  for (i=FSize;i<FSizeOfArray;i++)
    FElement[i]=NULL;
}

MkElements::~MkElements()
{
  Clear();
}

void MkElements::Initialize(int size,MkElement **elems)
{
  int i;
  if (size < 0) {
    MkDebug("::MkElements - MkElements(int size)");
    return;
  }

  Clear();

  FSizeOfArray = FSize = size;
  if (FSize == 0) {
     FElement = NULL;
     return;
  }

  FElement = new MkElement*[FSize];
  assert(FElement);

  for (i=0;i<FSize;i++)
    dyn_cp(FElement[i],elems[i]);
  for (i=FSize;i<FSizeOfArray;i++)
    FElement[i]=NULL;
}

bool MkElements::Add(MkElement *e)
{
  int i,j,k,l,size=FSize;
  bool flag=false;

  MkElement *elem;
  MkElement *elem_t=NULL;
  char str[256];

  if(FSizeOfArray==0) {
    FElement = new MkElement*[1];
    if(!FElement) return false;
    dyn_cp(FElement[0],e);
    FSizeOfArray = FSize = 1;
    return true;
  }

  dyn_cp(elem,e);
                   
  MkNodes &nodes  = FElement[0]->GetNodes();
  MkNodes &node   = elem->GetNodes();
  MkInt &elemnode = elem->GetElemNode();

  for (i=0;i<FSize;i++)
    if (*FElement[i]==*elem) flag = true;

  if(flag) return false;

  for (k=0;k<elemnode.getSzX();k++) {
    flag = true;
    for (j=0;j<FSize;j++) {
      if(FElement[j]->IsIn(node[elemnode(k)]) &&
         FElement[j]->GetElemNode(0)!=node[elemnode(k)] &&
         FElement[j]->GetElemNode(1)!=node[elemnode(k)]) {
         node[elemnode(k)].SetMovable(false);
         flag = false;
         FElement[j]->Out();
         node[elemnode(k)].Out();
      }
      if (!flag) {
        dyn_cp(elem_t,FElement[j]);

        MkInt &en1 = FElement[j]->GetElemNode();
        MkInt &en2 = elem_t->GetElemNode();

        nodes.Add(node[elemnode(k)]);

        for(l=0;l<nodes.GetSize();l++) {
          if(nodes[l]==node[elemnode(k)]) break;
        }

        if (l < nodes.GetSize()) {
          en1(1) = l;
          sprintf(str,"en1(1) = %d, l = %d\n",en1(1),l);
          MkDebug(str);
          en2(0) = l;
          sprintf(str,"en2(0) = %d, l = %d\n",en2(0),l);
          MkDebug(str);
        }
        else MkDebug("MkElements::Add(MkElement *) something wrong!!!\n");

        size = FSize;
	if(FSize==FSizeOfArray) Grow(10);
	if(FSize==FSizeOfArray) {
	  DeleteAElement(elem);
	  DeleteAElement(elem_t);
	  return false;
	}
        FSize = size+1;
        dyn_cp(FElement[size],elem_t);

        DeleteAElement(elem_t);
        break;
      }
    }
  }

  for (k=0;k<elemnode.getSzX();k++) {
    flag = true;
    for (j=0;j<nodes.GetSize();j++) {
      if(nodes[j]==node[elemnode(k)]) {
        nodes[j].SetMovable(false);
        elemnode(k) = j;
        flag = false;
        break;
      }
    }
    if(flag && k<elemnode.getSzX()) {
      nodes.Add(node[elemnode(k)]);
      for (j=nodes.GetSize()-1;j>=0;j++) {
        if(nodes[j]==node[elemnode(k)]) {
          elemnode(k) = j;
          break;
        }
      }
    }
  }

  elem->SetNodes(nodes);

  size = FSize;
  if(FSize==FSizeOfArray) Grow(10);
  if(FSize==FSizeOfArray) {
    DeleteAElement(elem);
    return false;
  }
  FSize = size+1;

  dyn_cp(FElement[size],elem);

  DeleteAElement(elem);
  return true;
}

// do not check the nodes, cross type. assume that the nodes are already in Nodes repositoty...
bool MkElements::AddNew(MkElement *e)  
{
  int i,size;
  bool flag=false;
  MkElement *elem;

  if(FSize==0) {
    FElement = new MkElement*[1];
    if(!FElement) return false;
    dyn_cp(FElement[0],e);
    FSizeOfArray = FSize = 1;
    return true;
  }

  dyn_cp(elem,e);

  for (i=0;i<FSize;i++)
    if (*FElement[i]==*elem) flag = true;

  if(flag) return false;

  size = FSize;
  if(FSize==FSizeOfArray) Grow(10);
  if(FSize==FSizeOfArray) {
    DeleteAElement(elem);
    return false;
  }

  FSize = size+1;
  dyn_cp(FElement[size],elem);

  DeleteAElement(elem);
  return true;
}

bool MkElements::Delete(MkElement *elem)
{
  int i,j,k,size=FSize;
  bool flag=false;
  MkElement **elems;

  for (i=0;i<FSize;i++)
    if (*FElement[i]==*elem) {flag = true; break;}

  if(flag) {
    elems = new MkElement*[FSize-1];
    assert(elems);

    for (j=0;j<i;j++) dyn_cp(elems[j],FElement[j]);
    for (j=i+1;j<FSize;j++) dyn_cp(elems[j-1],FElement[j]);

    Clear();
    Initialize(size-1,elems);
  }
  return true;
}

bool MkElements::DeleteAElement(MkElement *elem)
{
  if(elem) {
    if(elem->isGeneralEmenent()) {
      MkGenericElement *l=(MkGenericElement*)elem;
	  delete l;		
	}
    else if(elem->isCubicElement()) {
      MkCubicElement *l=(MkCubicElement*)elem;
	  delete l;		
	}
    else if(elem->isPrismElement()) {
      MkPrismElement *l=(MkPrismElement*)elem;
	  delete l;		
	}
    else if(elem->isTetraElement()) {
      MkTetraElement *l=(MkTetraElement*)elem;
	  delete l;		
	}
    else if(elem->isQuadElement()) {
      MkQuadElement *l=(MkQuadElement*)elem;
	  delete l;		
	}
    else if(elem->isTriElement()) {
      MkTriElement *l=(MkTriElement*)elem;
	  delete l;		
	}
    else if(elem->isBeamElement()) {
      MkBeamElement *l=(MkBeamElement*)elem;
	  delete l;		
	}
    else if(elem->isTrussElement()) {
      MkTrussElement *l=(MkTrussElement*)elem;
	  delete l;
	}
    else if(elem->isLineElement()) {
      MkLineElement *l=(MkLineElement*)elem;
      delete l;
	}
	else delete elem;
    elem = NULL;
	return true;
  }
  return false;
}


int MkElements::Grow(int delta)
{
    int i,size,sizearray;
    MkElement **element=NULL;

    if (!(element = new MkElement*[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
      dyn_cp(element[i],FElement[i]);
    for (i=FSize; i<FSizeOfArray+delta;i++)
      element[i]=NULL;

    size = FSize;
    sizearray = FSizeOfArray;

    Clear();

    FElement = element;
	FSize = size;
    FSizeOfArray = sizearray+delta;

    return FSizeOfArray;
}

int MkElements::Shrink(int delta)
{
    int i,size,sizearray;
    MkElement **element=NULL;

    if (!(element = new MkElement*[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
      dyn_cp(element[i],FElement[i]);
    for (i=FSize; i<FSizeOfArray-delta;i++)
      element[i]=NULL;

	size = FSize;
	sizearray = FSizeOfArray;

	Clear();

    FElement = element;
    FSizeOfArray = sizearray-delta;
	FSize = size;

    return FSizeOfArray;
}

bool MkElements::Clear()
{
 if (FElement) {
   for(int i=0;i<FSize;i++) {
     if(FElement[i]) {
		DeleteAElement(FElement[i]);
     }
   }
   delete[] FElement;
   FElement = NULL;
   FSizeOfArray = 0;
   FSize = 0;
   return true;
 }
 else return false;
}

MkElements & MkElements::operator=(MkElements &elems)
{
    int i;

    Clear();
    FSizeOfArray = elems.FSizeOfArray;
    FSize = elems.FSize;
    if (FSize == 0) {
       FElement = NULL;
       return *this;
    }
    FElement = new MkElement*[FSizeOfArray];
    assert(FElement);

    for (i=0;i<FSize;i++)
      dyn_cp(FElement[i],&elems[i]);
    for (i=FSize;i<FSizeOfArray;i++)
      FElement[i] = NULL;

    return *this;
}

MkElement & MkElements::operator[](int i)
{
    if (FSizeOfArray == 0) return NullElement;
    if (i >= FSize && i < FSizeOfArray) {FSize = i+1;return *FElement[i];}
    else if (i >=0 && i < FSize) return *FElement[i];
    else return NullElement;
}

bool MkElements::operator==(MkElements &e)
{
  int i;
  bool flag=true;
  if (FSize!=e.FSize) return false;
  for (i=0;i<FSize;i++) {
    flag = flag && *FElement[i]==e[i];
    if(!flag) return flag;
  }
  return flag;
}

bool MkElements::operator!=(MkElements &e)
{
  return !(*this==e);
}

#ifdef __BCPLUSPLUS__
void MkElements::Out(TObject *Sender)
{
  int i;
  if (FSize==0) return;
  if (String(Sender->ClassName()) == String("TMemo")) {
    TMemo *memo = dynamic_cast<TMemo*>(Sender);
    if(!memo) return;
    for (i=0;i<FSize;i++) FElement[i]->Out(memo);
  }
}
#else
void MkElements::Out(char *fname)
{
  int i;
  if (FSize==0) return;
  FILE *fp=fopen(fname,"w");
  fclose(fp);
  for (i=0;i<FSize;i++) FElement[i]->Out(fname);
}
#endif
void MkElements::Out()
{
  int i;
  if (FSize==0) return;
  for (i=0;i<FSize;i++) FElement[i]->Out();
}

#ifdef __BCPLUSPLUS__
void MkElements::Draw(TObject *Sender)
{
  for(int i=0;i<FSize;i++) 
    FElement[i]->Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkElements::Draw(MkPaint *pb)
{
  for(int i=0;i<FSize;i++)
    FElement[i]->Draw(pb);
}
#endif

//--------------------------------------------------------------------
MkBound::MkBound()
{
  FMesh = NULL;
  FBoundType = btOneDimension;
  FBoundGroup= 0;
  NextBound = NULL;
  PrevBound = NULL;
  className = "MkBound";
}

MkBound::MkBound(int n)
{
  FMesh = NULL;
  FBoundType = btOneDimension;
  FBoundGroup= 0;
  NextBound = NULL;
  PrevBound = NULL;
  className = "MkBound";
}

MkBound::MkBound(MkBoundType bt)
{
  FMesh = NULL;
  FBoundType = bt;
  FBoundGroup= 0;
  NextBound = NULL;
  PrevBound = NULL;
  className = "MkBound";
}

MkBound::~MkBound()
{
  if(FMesh) delete FMesh;
}

bool MkBound::IsNeighborWith(MkBound &bnd)
{
  return false;
}

//bool MkBound::BoundMesh()
//{
//}

MkBounds::MkBounds(int size,MkBound *bnds)
{

    if (size < 0) {
      MkDebug("::MkBounds - MkBounds(int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0) {
       FBound = NULL;
       return;
    }

    FBound = new MkBound[FSize];
    for (int i=0;i<FSize;i++) (*this)[i] = bnds[i];
}

MkBounds::MkBounds(int size)
{
    if (size < 0) {
      MkDebug("::MkBounds - MkBounds(int size)");;
      return;
    }

    FSizeOfArray = size;
    FSize = 0;
    if (FSizeOfArray == 0) {
       FBound = NULL;
       return;
    }

    FBound = new MkBound[FSizeOfArray];
}

MkBounds::~MkBounds()
{
   FSizeOfArray = FSize = 0;
   if (FBound) {
      delete[] FBound;
      FBound = NULL;
   }
}

void MkBounds::Initialize(int size)
{
    if (size < 0) {
      MkDebug("::MkBounds - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;

    FSizeOfArray = size;
    FSize = 0;
    
    if (FSizeOfArray == 0) {
       if (FBound!=NULL) delete[] (MkBound*)FBound;
       FBound = NULL;
       return;
    }

    if (FBound!=NULL) delete[] (MkBound*)FBound;
    FBound = new MkBound[FSizeOfArray];
}

void MkBounds::Initialize(int size,MkBound *bnds)
{

    if (size < 0 || bnds == NULL) {
      MkDebug("::MkBounds - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0) {
       if (FBound!=NULL) delete[] (MkBound*)FBound;
       FBound = NULL;
       return;
    }

    if (FBound!=NULL) delete[] (MkBound*)FBound;
    FBound = new MkBound[FSizeOfArray];
    for (int i=0;i<FSizeOfArray;i++) FBound[i] = bnds[i];
}

int MkBounds::Grow(int delta)
{
    int i;
    MkBound *bnd=NULL;

    if (!(bnd = new MkBound[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        bnd[i] = FBound[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        bnd[i] = NullBound;
    if (FBound) {
       delete[] (MkBound*)FBound;
       FBound = NULL;
    }
    FBound = bnd;
    FSizeOfArray = FSizeOfArray+delta;
    return FSizeOfArray;
}

int MkBounds::Shrink(int delta)
{
    int i;
    MkBound *bnd=NULL;

    if (!(bnd = new MkBound[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        bnd[i] = FBound[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        bnd[i] = NullBound;
    if (FBound) {
       delete[] (MkBound*)FBound;
       FBound = NULL;
    }
    FBound = bnd;
    FSizeOfArray = FSizeOfArray-delta;
    return FSizeOfArray;
}

bool MkBounds::Add(MkBound &bnd)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) Grow(FSize-FSizeOfArray+1);
    if(tmp==FSizeOfArray) return false;

    FSize++;
    FBound[FSize-1] = bnd;
    return true;
}

bool MkBounds::Add(int index, MkBound &bnd)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) Grow(FSize-FSizeOfArray+1);
    if(tmp==FSizeOfArray) return false;

    for (int i=FSize-1;i>=index;i--)
      FBound[i+1] = FBound[i];
    FSize++;
    FBound[index] = bnd;
    return true;
}

bool MkBounds::Delete(MkBound &bnd)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FBound[i] == bnd) break;
    }
    if(i==FSize) return false;
    if(FBound[i] == bnd) {
      for (int j=i;j<FSize-1;j++)
        FBound[j] = FBound[j+1];
    }
    FSize--;
    FBound[FSize] = NullBound;
    return true;
}

bool MkBounds::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FBound[j] = FBound[j+1];

    FSize--;
    FBound[FSize] = NullBound;
    return true;
}

bool MkBounds::Clear()
{
   FSizeOfArray = FSize = 0;
   if (FBound) {
      delete[] FBound;
      FBound = NULL;
   }
   return true;
}

MkBound & MkBounds::operator[](int i)
{
    if (FSize == 0) return NullBound;
    else if (i >=0 && i < FSize) return FBound[i];
    else return NullBound;
}

MkBounds & MkBounds::operator=(MkBounds &bnds)
{
    int i;

    Clear();
    FSize = bnds.FSize;
    FSizeOfArray = bnds.FSizeOfArray;
    if (FSize == 0) {
       FBound = NULL;
       return *this;
    }
    this->FBound = new MkBound[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FBound[i] = bnds.FBound[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FBound[i] = NullBound;

    return *this;
}

bool MkBounds::operator==(MkBounds &bnds)
{
    int i;

    if (FSize != bnds.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FBound[i] != bnds.FBound[i]) return false;

    return true;
}

#ifdef __BCPLUSPLUS__
void MkBounds::Draw(TObject *Sender)
{
}
#endif



//--------------------------------------------------------------------
MkMesh::MkMesh()
{
  className = "MkMesh";
}

MkMesh::MkMesh(int n)
{
  className = "MkMesh";
}

MkMesh::~MkMesh()
{

}

bool MkMesh::SetBoundary(MkBounds &bounds)
{
  return false;
}

bool MkMesh::Generate()
{
  int mt;
#ifdef __BCPLUSPLUS__
  AnsiString className;
  className=ClassName();
  if (className=="MkMesh") mt = 0;
  else if (className=="MkOneDimMesh") mt = 1;
  else if (className=="MkTwoDimMesh") mt = 2;
  else if (className=="MkThreeDimMesh") mt = 3;
#else
  if (!className.compare("MkMesh")) mt = 0;
  else if (!className.compare("MkOneDimMesh")) mt = 1;
  else if (!className.compare("MkTwoDimMesh")) mt = 2;
  else if (!className.compare("MkThreeDimMesh")) mt = 3;
#endif
  if(mt==0) {

  }
  else if(mt==1) {

  }
  else if(mt==2) {

  }
  else if(mt==3) {

  }
  return true;
}

bool CheckOverlap(MkLine &line, MkLines &lines)
{
  int i,j;
  for (i=0;i<lines.GetSize();i++) {
    if ((line[0]!=lines[i][0]) && (line[0]!=lines[i][1]) &&
        (line[1]!=lines[i][0]) && (line[1]!=lines[i][1]) &&
        (line && lines[i])) return true;
  }
  return false;
}

bool CheckOverlap(MkLine &line, MkTriangles &tri)
{
  int i,j;
  for (i=0;i<tri.GetSize();i++) {
    for (j=0;j<3;j++) {
      if ((line[0]!=tri[i](j)[0]) && (line[0]!=tri[i](j)[1]) &&
          (line[1]!=tri[i](j)[0]) && (line[1]!=tri[i](j)[1]) &&
          (line && tri[i](j))) return true;
    }
  }
  return false;
}

bool CheckOverlap(MkPoint &pnt, MkTriangles &tri)
{
  int i;
  for (i=0;i<tri.GetSize();i++) {
   if(tri[i].isInside(pnt)) return true;
  }
  return false;
}

bool ResetLines(MkLines &lines, MkTriangle &tri)
{
  int i,j,k,index[3];
  MkLine l[2];
  index[0] = index[1] = index[2] = -1;  
  for(i=0;i<lines.GetSize();i++) {
    l[0] = lines[i];
    l[1].SetLine(l[0][1],l[0][0]);
    for(j=0;j<3;j++) {
      if((tri(j)==l[0]) || (tri(j)==l[1])) {
        lines[i].Select();
        index[j] = i;
      }
    }
  }
  lines.DeleteSelected();
  for(j=0;j<3;j++) {
    if(index[j]==-1) {
      MkLine line;
      line.SetLine(tri(j)[1],tri(j)[0]);
      lines.Add(line);
    }
  }
}

bool ResetPoints(MkPoints &pnts, MkLines &lines)
{
  int i,j;
  bool found;
  for (i=0;i<pnts.GetSize();i++) {
    found = false;
    for (j=0;j<lines.GetSize();j++) {
      if(pnts[i]==lines[j][0] || pnts[i]==lines[j][1]) {
        found=true;
        break;
      }
    }
    if(!found) {
      MkPoint pnt = pnts[i];
      pnts.Delete(pnt);
      i--;
    }
  }
  return true;
}

bool MkMesh::GenerateWithPolygon(MkPolygon &poly)
{
  int i, j, k, step;
  bool flag=true, isfirst;
  MkPoints pnts, nodepnt;
  MkLine *line, l, li[2];
  MkLines lines;
  MkTriangles triangles;
  MkTriangle tri, t;
  MkVector norm;

  line = new MkLine[poly.GetSize()];
  if (!line) return false;
  pnts.Initialize(poly.GetSize(),poly.GetPoints());

  for (i=0;i<poly.GetSize();i++) {
    line[i].SetLine(poly[i],poly[i+1==poly.GetSize()?0:i+1]);
    if(fabs(norm[0]+norm[1]+norm[2])<EPS) line[0].GetVector().Cross(line[i].GetVector(),norm);
  }
  norm.Normalize();

  lines.Initialize(poly.GetSize(),line);

  delete[] line;
  line = NULL;

  step = 0;
  while (flag) {
    if(lines.GetSize()==0) break;
    l = lines[0];
    isfirst = true;
    for(i=0;i<pnts.GetSize();i++) {
      if(l[0]==pnts[i] || l[1]==pnts[i]) continue;
      li[0].SetLine(l[0],pnts[i]);
      li[1].SetLine(l[1],pnts[i]);

      tri.Reset(l[0],l[1],pnts[i]);
      if(isfirst) {t = tri;isfirst = false;continue;}
      if(!CheckOverlap(li[0],lines) && !CheckOverlap(li[1],lines) &&
         !CheckOverlap(li[0],triangles) && !CheckOverlap(li[1],triangles) &&
         tri.GetRadius() < t.GetRadius() && (tri.GetNormal()*norm) > 0) t = tri;
    }

    triangles.Add(t);
    ResetLines(lines,t);
    ResetPoints(pnts,lines);

    flag = flag && lines.GetSize()>3;
    if (step>1000) flag = false;
    step++;
  }
  if(lines.GetSize()==3) {
    t.Reset(lines[0][0],lines[1][0],lines[2][0]);
    triangles.Add(t);
  }
  else return false;
  if(step>1000) return false;

  pnts.Clear();
  for(i=0;i<triangles.GetSize();i++) {
    for(j=0;j<3;j++) {
      bool flag=true;
      for(k=0;k<pnts.GetSize();k++) {
        if(triangles[i][j]==pnts[k]) flag = false;
      }
      if(flag) pnts.Add(triangles[i][j]);
    }
  }

  MkDOFs dof(6);
  dof[0].SetType(doftXDis,bndtFree); // temp
  dof[1].SetType(doftYDis,bndtFree);
  dof[2].SetType(doftZDis,bndtFree);
  dof[3].SetType(doftXAng,bndtFree);
  dof[4].SetType(doftYAng,bndtFree);
  dof[5].SetType(doftZAng,bndtFree);

  Nodes.Initialize(pnts.GetSize());
  for(i=0;i<pnts.GetSize();i++) {
    Nodes[i].SetPoint(pnts[i]);
    Nodes[i].SetDOF(dof);
  }

  MkElement **telem;
  telem = new MkElement*[triangles.GetSize()];
  if(!telem) {
    Nodes.Clear();
    return false;
  }
  for(i=0;i<triangles.GetSize();i++)
    telem[i] = new MkTriElement();

  MkInt elemnode(3);
  for(i=0;i<triangles.GetSize();i++) {
    for(j=0;j<3;j++) {
      elemnode(j) = -1;
      for(k=0;k<pnts.GetSize();k++) {
        if(triangles[i][j] == pnts[k]) elemnode(j) = k;
      }
      if(elemnode(j) == -1) {
        for(i=0;i<triangles.GetSize();i++)
          if(telem[i]) {delete telem[i];telem[i] = NULL;}
        if(telem) {delete[] telem;telem = NULL;}
        Nodes.Clear();
        return false;//strange;
      }
    }
    telem[i]->SetElemNode(elemnode);
    telem[i]->SetNodes(Nodes);    
  }

  Elements.Initialize(triangles.GetSize(),telem);

  for (i=0;i<triangles.GetSize();i++) if(telem[i]) {delete telem[i];telem[i]=NULL;}
  if(telem) {delete[] telem; telem = NULL;}
  
  return true;
}

MkMesh & MkMesh::operator=(MkMesh &mesh)
{
  Elements = mesh.Elements;
  Nodes=mesh.Nodes;
  Layers=mesh.Layers;
  Bounds=mesh.Bounds;
  BoundNodes=mesh.BoundNodes;
  InternalNodes=mesh.InternalNodes;
  Front=mesh.Front;
  isGenerated=mesh.isGenerated;
  strcpy(FName,mesh.FName);
  Division=mesh.Division;
  return *this;
}

MkMeshes::MkMeshes(int size,MkMesh *meshes)
{
    if (size < 0) {
      MkDebug("::MkMeshes - MkMeshes(int size)");;
      return;
    }

    FSize = size;
    if (FSize == 0) {
       FMesh = NULL;
       return;
    }

    FMesh = new MkMesh[FSize];
    for (int i=0;i<FSize;i++) (*this)[i] = meshes[i];
}

MkMeshes::MkMeshes(int size)
{
    if (size < 0) {
      MkDebug("::MkMeshes - MkMeshes(int size)");;
      return;
    }

    FSize = size;
    if (FSize == 0) {
       FMesh = NULL;
       return;
    }

    FMesh = new MkMesh[FSize];
}

MkMeshes::~MkMeshes()
{
   FSize = 0;
   if (FMesh) {
      delete[] FMesh;
      FMesh = NULL;
   }
}

void MkMeshes::Initialize(int size)
{
    if (size < 0) {
      MkDebug("::MkMeshes - Initialize(int size)");;
      return;
    }
    if (FSize == size) return;

    FSize = size;

    if (FSize == 0) {
       if (FMesh!=NULL) delete[] (MkMesh*)FMesh;
       FMesh = NULL;
       return;
    }

    if (FMesh!=NULL) delete[] (MkMesh*)FMesh;
    FMesh = new MkMesh[FSize];
}

void MkMeshes::Initialize(int size,MkMesh *meshes)
{

    if (size < 0 || meshes == NULL) {
      MkDebug("::MkMeshes - Initialize(int size)");;
      return;
    }
    if (FSize == size) return;
    FSize = size;
    if (FSize == 0) {
       if (FMesh!=NULL) delete[] (MkMesh*)FMesh;
       FMesh = NULL;
       return;
    }

    if (FMesh!=NULL) delete[] (MkMesh*)FMesh;
    FMesh = new MkMesh[FSize];
    for (int i=0;i<FSize;i++) FMesh[i] = meshes[i];
}

bool MkMeshes::Clear()
{
   FSize = FSize = 0;
   if (FMesh) {
      delete[] FMesh;
      FMesh = NULL;
   }
   return true;
}

MkMesh & MkMeshes::operator[](int i)
{
    if (FSize == 0) return NullMesh;
    else if (i >=0 && i < FSize) return FMesh[i];
    else return NullMesh;
}

MkMeshes & MkMeshes::operator=(MkMeshes &meshes)
{
    int i;

    Clear();
    FSize = meshes.FSize;
    if (FSize == 0) {
       FMesh = NULL;
       return *this;
    }
    this->FMesh = new MkMesh[FSize];

    for (i=0;i<FSize;i++)
      FMesh[i] = meshes.FMesh[i];

    return *this;
}

bool MkMeshes::operator==(MkMeshes &meshes)
{
    int i;

    if (FSize != meshes.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FMesh[i] != meshes.FMesh[i]) return false;

    return true;
}

#ifdef __BCPLUSPLUS__
void MkMeshes::Draw(TObject *Sender)
{

}
#endif
//--------------------------------------------------------------------
