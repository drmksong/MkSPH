//---------------------------------------------------------------------------
#ifndef MkLayerH
#define MkLayerH

#ifdef __BCPLUSPLUS__
#include "GlobalVarUnit.h"
#endif

#include "MkObject.h"
#include "MkEntity.h"
#include "MkCut.h"
#include "MkFill.h"
#include "MkProfile.h"
//---------------------------------------------------------------------------
class MkLayer : public MkEntity {
protected:
  MkSide Side;
  MkRect Rect;
  MkProfile Profile;
  float R_depth;
  float WetUnitWeight[2];
  float SubUnitWeight[2];
  float Cohesion[2];
  float Friction[2];
  float HorSubReact[2]; // 수평 지반반력계수(modulus of subgrade reaction
  float VerSubReact[2]; // 수직 지반반력계수(modulus of subgrade reaction
  float Burden;
  float K0;
  float Ka;
  float Kp;
  float Bearing;

//  AnsiString layername;    MkEntity::Name
//  double depth,R_depth;    MkEntity::Depth
//  double Rt,Rsub,C,Ks,Pi;  WetUnitWeight, SubUnitWeight, Cohesion, HorSubReact, Friction
//  double bearing;          Bearing
//  MkFloat Profile;         //Polygon

public:
  MkLayer();
  MkLayer(int n);
  ~MkLayer(){};
public:
  void SetSide(MkSide side){Side = side;}
  void SetRect(MkRect &rect){Rect=rect;}
  void SetK0(float k0){K0 = k0;}
  void SetKa(float ka){Ka = ka;}
  void SetKp(float kp){Kp = kp;}
  void SetBurden(float burden){Burden=burden;}
  void SetR_depth(float d){R_depth=d;}
  void SetWetUnitWeight(int n,float w){if(n>=0&&n<2) WetUnitWeight[n] = w;}
  void SetSubUnitWeight(int n,float w){if(n>=0&&n<2) SubUnitWeight[n] = w;}
  void SetCohesion(int n,float w){if(n>=0&&n<2) Cohesion[n] = w;}
  void SetFriction(int n,float w){if(n>=0&&n<2) Friction[n] = w;}
  void SetHorSubReact(int n,float w){if(n>=0&&n<2) HorSubReact[n] = w;}
  void SetVerSubReact(int n,float w){if(n>=0&&n<2) VerSubReact[n] = w;}
  void SetBearing(float b){Bearing = b;}
//  void SetProfile(MkFloat &prof){Profile.CopyFrom(prof);}
  void SetProfile(MkProfile &prof){Profile = prof;}
public:
  MkSide GetSide(){return Side;}
  float GetActivPress(MkPoint &pnt);
  float GetPassivPress(MkPoint &pnt);
  float GetStopPress(MkPoint &pnt);
  float GetWetUnitWeight(MkPoint &pnt);
  float GetCohesion(MkPoint &pnt);
  float GetFriction(MkPoint &pnt);
  float GetHorSubReact(MkPoint &pnt);
  float GetVerSubReact(MkPoint &pnt);
  float GetK0(MkPoint &pnt);
  float GetKa(MkPoint &pnt);
  float GetKp(MkPoint &pnt);
  float GetR_depth(){return R_depth;}
  float GetWetUnitWeight(MkLine &line);
  float GetCohesion(MkLine &line);
  float GetFriction(MkLine &line);
  float GetHorSubReact(MkLine &line);
  float GetVerSubReact(MkLine &line);
  float GetK0(MkLine &line);
  float GetKa(MkLine &line);
  float GetKp(MkLine &line);
  float GetBearing(){return Bearing;}
//  MkFloat &GetProfile(){return Profile;}
  MkProfile &GetProfile(){return Profile;}

  float GetWetUnitWeight(int i){return ((i<0||i>2) ? 0: WetUnitWeight[i]);}
  float GetSubUnitWeight(int i){return ((i<0||i>2) ? 0: SubUnitWeight[i]);}
  float GetCohesion(int i){return ((i<0||i>2) ? 0: Cohesion[i]);}
  float GetFriction(int i){return ((i<0||i>2) ? 0: Friction[i]);}
  float GetHorSubReact(int i){return ((i<0||i>2) ? 0: HorSubReact[i]);}
  float GetVerSubReact(int i){return ((i<0||i>2) ? 0: VerSubReact[i]);}

public: //key function
  MkRect &GetRect(){return Rect;}
//  float GetWeight(){return (WetUnitWeight[0]+WetUnitWeight[1])/2.0*Rect.GetHeight();}
  bool IsIn(MkPoint &pnt){return Rect.IsInSurface(pnt,(float)EPS) || Rect.IsIn(pnt);}
  bool IsIn(MkLine &line){return Rect.IsIn(line[0]) || Rect.IsIn(line[1]);}
  float InLen(MkLine &line);
  MkLine InLine(MkLine &line);
  void CalcCoeff();

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkLayer");}
  bool UpdateFrom(TStringGrid *grid){Grid = grid;return UpdateFrom();};
  bool UpdateTo(TStringGrid *grid){Grid = grid;return UpdateTo();};
  bool UpdateFrom();
  bool UpdateTo();
  void Out(TObject *);
  void Import(MkGlobalVar &globalvar, int sec, MkSide ls, int lay);
  void Import(MkGlobalVar &globalvar, int sec, int lay);
  void Export(MkGlobalVar &globalvar, int sec, int lay);
#else
  std::string ClassName(){return className;}
#endif
  void Out(char *fname);

#ifdef __BCPLUSPLUS__
  void  Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif

public: // operator function
  bool operator==(MkLayer&);
  bool operator!=(MkLayer&);
  MkLayer & operator=(MkLayer&);
};

class MkLayers {
protected:
  MkLayer *FLayer;
  int FSize;//Actual size of entities
  int FSizeOfArray;
#ifdef __BCPLUSPLUS__
  TColor Color;
#endif
public:
  MkLayers(int size,MkLayer *ent);
  MkLayers(int size);
  MkLayers(){FSizeOfArray = FSize = 0;FLayer = NULL;}
   ~MkLayers();
  virtual void Initialize(int size);
  virtual void Initialize(int size,MkLayer *);

  bool Add(MkLayer &layer);  // change of size of ground
  bool Add(int index,MkLayer &layer);
  bool Delete(MkLayer &layer);  // change of size of ground
  bool Delete(int index);
  int Grow(int Delta);            // change of size of array
  int Shrink(int Delta);          // change of size of array

  void SetupPress();
  int GetSize(){return FSize;};
  int GetNumber(){return FSize;};

//  bool Apply(MkCut &);
//  bool Apply(MkFill &);

  bool Clear();
  virtual MkLayer & operator[](int);


  MkLayers & operator=(MkLayers &layers);
  bool operator==(MkLayers &layers);
#ifdef __BCPLUSPLUS__
  TColor GetColor(){return Color;};
  void SetColor(TColor c){Color = c;}
  void Out(TObject *);
  void Import(MkGlobalVar &globalvar, int sec);
  void Export(MkGlobalVar &globalvar, int sec);
#else
#endif

  void Out(char *fname);

  float GetStopPress(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetStopPress(pnt);
      return 0;
    }
  float GetActivPress(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetActivPress(pnt);
      return 0;
    }
  float GetPassivPress(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetPassivPress(pnt);
      return 0;
    }
  float GetWetUnitWeight(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetWetUnitWeight(pnt);
      return 0;
    }
  float GetCohesion(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetCohesion(pnt);
      return 0;
    }
  float GetFriction(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetFriction(pnt);
      return 0;
    }
  float GetHorSubReact(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetHorSubReact(pnt);
      return 0;
    }
  float GetVerSubReact(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetVerSubReact(pnt);
      return 0;
    }
  float GetK0(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetK0(pnt);
      return 0;
    }
  float GetKa(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetKa(pnt);
      return 0;
    }
  float GetKp(MkPoint &pnt)
    {
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(pnt))
          return FLayer[i].GetKp(pnt);
      return 0;
    }
  float GetWetUnitWeight(MkLine &line)
    {
      float w=0,l;
      if(line.GetLength()<0.001) return GetWetUnitWeight(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          w+=FLayer[i].GetWetUnitWeight(line)*l;
        }
      return w/line.GetLength();
    }
  float GetCohesion(MkLine &line)
    {
      float c=0,l;
      if(line.GetLength()<0.001) return GetCohesion(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          c+= FLayer[i].GetCohesion(line)*l;
        }
       return c/line.GetLength();
    }
  float GetFriction(MkLine &line)
    {
      float f=0,l;
      if(line.GetLength()<0.001) return GetFriction(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          f+=FLayer[i].GetCohesion(line)*l;
        }
      return f/line.GetLength();
    }
  float GetHorSubReact(MkLine &line)
    {
      float h=0,l;
      if(line.GetLength()<0.001) return GetHorSubReact(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          h += FLayer[i].GetHorSubReact(line)*l;
        }
      return h/line.GetLength();
    }
  float GetVerSubReact(MkLine &line)
    {
      float v=0,l;
      if(line.GetLength()<0.001) return GetVerSubReact(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          v += FLayer[i].GetHorSubReact(line)*l;
        }
      return v/line.GetLength();
    }
  float GetK0(MkLine &line)
    {
      float v=0,l;
      if(line.GetLength()<0.001) return GetVerSubReact(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          v += FLayer[i].GetK0(line)*l;
        }
      return v/line.GetLength();
    }
  float GetKa(MkLine &line)
    {
      float v=0,l;
      if(line.GetLength()<0.001) return GetVerSubReact(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          v += FLayer[i].GetKa(line)*l;
        }
      return v/line.GetLength();
    }
  float GetKp(MkLine &line)
    {
      float v=0,l;
      if(line.GetLength()<0.001) return GetVerSubReact(line[0]);
      for (int i=0;i<FSize;i++)
        if(FLayer[i].IsIn(line)) {
          l = FLayer[i].InLen(line);
          v += FLayer[i].GetKp(line)*l;
        }
      return v/line.GetLength();
    }

#ifdef __BCPLUSPLUS__
  void  Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};
#endif
//---------------------------------------------------------------------------
extern MkLayer NullLayer;
extern MkLayers NullLayers;
//---------------------------------------------------------------------------
/*
class MkLayer : public MkObject {
public:
  AnsiString layername;
  double depth,R_depth;
  double Rt,Rsub,C,Ks,Pi;
  double bearing;
  MkFloat Profile;
public:
  MkLayer();
  MkLayer(int);
  ~MkLayer(){}
  void Clear();
  void Import(MkGlobalVar &globalvar, int sec, MkSide ls, int lay);
  void Export(MkGlobalVar &globalvar, int sec, MkSide ls, int lay);
};

class MkLayers {
protected:
  MkLayer *FLayer;
  int FSize;//Actual size of nodes
  int FSizeOfArray;
public:
  MkLayers(int size,MkLayer *layer);
  MkLayers(int size);
  MkLayers(){FSizeOfArray = FSize = 0;FLayer = NULL;}
  ~MkLayers();

  virtual void Initialize(int size);
  virtual void Initialize(int size,MkLayer *);
  int GetSize(){return FSize;};
  int GetNumber(){return FSize;};
  bool Add(MkLayer &layer);  // change of size of layer
  bool Add(int index,MkLayer &layer);
  bool Delete(MkLayer &layer);  // change of size of layer
  bool Delete(int index);
  int Grow(int Delta);            // change of size of array
  int Shrink(int Delta);          // change of size of array
  bool Clear();
  virtual MkLayer & operator[](int);
  MkLayers & operator=(MkLayers &layers);
  bool operator==(MkLayers &layers);

  void Import(MkGlobalVar &globalvar, int sec,MkSide ls);
  void Export(MkGlobalVar &globalvar, int sec,MkSide ls);
};
//---------------------------------------------------------------------------
extern MkLayer NullMkLayer;
*/
