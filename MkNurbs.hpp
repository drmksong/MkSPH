//---------------------------------------------------------------------------
#ifndef MkNurbsH
#define MkNurbsH

#include <nurbh.h>
#include <nurbtype.h>
#include <nurbfns.h>
//#include <tran.h>
#include <trantype.h>
#include <tranfns.h>
#include "MkFloat.h"
#include "MkInt.h"
#include "MkPoint.h"

#include <Inventor/SbVec3f.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoSphere.h>

enum NurbsShape {nsNone, nsArc, nsCircle, nsPolyline, nsTorus, nsTriangle,
                 nsSquare, nsSphere, nsVase};

enum NurbsDirection {ndU, ndV };

extern bool isNurbInitiaized;

class MkNurbs {
private:
  PR_nurb FNurbs;
  
private:
  void AllocateKnots(NurbsDirection nd);
  void AllocatePoints();
  void Boehmc(PR_dir dir, PR_nurb &new_nurb);
  void Bvalue(float u, int deriv, PR_nurb &new_nurb);
  void Clear();
  void Copy(PR_nurb &nurb);
  void CopyKnots(NurbsDirection nd, PR_nurb &nurb);
  void DeallocateKnots(NurbsDirection nd);
  void DeallocatePoints(void);
  void DifferenceKnots(NurbsDirection nd, PR_dir &k1, PR_dir &k2);
  void Du(float x, MkPoints &rps);
  void Dump(char *fname);
  void DumpKnots(char *fname, NurbsDirection nd, char *str);
  void DumpLists(char *fname);
  void DumpSummary(char *fname);
  void Elevate();
  MkPoint Evaluate(float u);
  MkPoint Evaluate(float u, float v);
  void Extrude(MkPoint rp);
  void Init(int kx,int ntx,int ky,int nty);
  int InqErrorCount();
  void Interchange();
  void Interpolate(PR_nurb &nurb, float tol);
  void MakeKnots(NurbsDirection nd,int order, int cpnum);
  void NumSubDivs(NurbsDirection nd, float tol);
  MkPoint PartialDU(float u, float v);
  MkPoint PartialDV(float u, float v);
  MkPoint PartialDUV(float u, float v);
  void Revolve();
  void Ruled(PR_nurb &nurb1, PR_nurb &nurb2, int &err);
  void ScaleKnots(float min, float max, NurbsDirection nd);
  void Split(float f, PR_nurb &nurb1, PR_nurb &nurb2);
  void Stats(char *fname);
  void Tessalate(int);
  void Transpose();
  void UnionKnots(PR_dir &k1, NurbsDirection nd);
  void Write(char *fname);
  void Xform(Pmatrix3 &mat);

protected:
  int FNumControlPoint;
  MkFloat FKnots;
  MkPoints FControlPoints;

public:
  MkNurbs();
  ~MkNurbs();

  MkPoints & GetControlPoint();
  void SetControlPoint(MkPoints &rps);
  MkFloat & GetKnot();
  void SetKnot(MkFloat & knots);

  void SetOrderNumber(NurbsDirection nd, int k){nd==ndU?FNurbs.pf_u.pf_k=k:FNurbs.pf_v.pf_k=k;}
  void SetControlPointNumber(NurbsDirection nd, int n){nd==ndU?FNurbs.pf_u.pf_n=n:FNurbs.pf_v.pf_n=n;}
  void SetKnotNumber(NurbsDirection nd, int nt){nd == ndU ? FNurbs.pf_u.pf_nt = nt:FNurbs.pf_v.pf_nt=nt;}

  // Shape
  void MakeItArc(float startAng, float endAng);
  void MakeItCircle();
  void MakeItPolyline(MkPoints &rps);
  void MakeItTorus();
  void MakeItTriangle();
  void MakeItSquare();
  void MakeItSphere();
  void MakeItVase();
  SoSeparator * MakeFaceSet();

  MkPoint operator()(float u, float v){return Evaluate(u,v);}
  MkPoint operator()(float u){return Evaluate(u);}
};
//---------------------------------------------------------------------------
#endif
