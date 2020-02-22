//---------------------------------------------------------------------------
#ifndef MkFillH
#define MkFillH
#include "MkEntity.h"
#include "MkWall.h"
#include "MkProfile.h"
//---------------------------------------------------------------------------
class MkFill : public MkEntity {
protected:
  float WetUnitWeight;
  float SubUnitWeight;
  float Cohesion;
  float Friction;
  float HorSubReact; // 수평 지반반력계수(modulus of subgrade reaction
  float VerSubReact; // 수직 지반반력계수(modulus of subgrade reaction
  float Burden;
  float K0;
  float Ka;
  float Kp;

  MkWall *Wall[2];
  MkProfile Profile;
//  MkPolygon FillLine;
public:
  MkFill();
    MkFill(int );
  ~MkFill(){}
  void SetWall(MkWall *wall){Wall[0] = wall;}
  void SetWall(MkWall *wall, int n){if (n<0||n>2) return;Wall[n]=wall;}
  void SetK0(float k0){K0 = k0;}
  void SetKa(float ka){Ka = ka;}
  void SetKp(float kp){Kp = kp;}
  void SetWetUnitWeight(float w){WetUnitWeight = w;}
  void SetSubUnitWeight(float w){SubUnitWeight = w;}
  void SetCohesion(float w){Cohesion = w;}
  void SetFriction(float w){Friction = w;}
  void SetHorSubReact(float w){HorSubReact = w;}
  void SetVerSubReact(float w){VerSubReact = w;}
  void SetProfile(MkProfile &prof){Profile = prof;}
//  void SetFillLine(MkPolygon &fillline){FillLine = fillline;}

public:
  float GetK0(){return K0;}
  float GetKa(){return Ka;}
  float GetKp(){return Kp;}
  float GetWetUnitWeight(){return WetUnitWeight;}
  float GetSubUnitWeight(){return SubUnitWeight;}
  float GetCohesion(){return Cohesion;}
  float GetFriction(){return Friction;}
  float GetHorSubReact(){return HorSubReact;}
  float GetVerSubReact(){return VerSubReact;}
  float GetDepth(void){return Depth;}
  float GetDepth(float x);
  MkProfile &GetProfile(){return Profile;}
//  MkPolygon &GetFillLine(){return FillLine;}

#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkFill");};
  bool UpdateFrom(TStringGrid* grid){Grid=grid;return UpdateFrom();}
  bool UpdateTo(TStringGrid* grid){Grid=grid;return UpdateTo();}
  bool UpdateFrom();
  bool UpdateTo();
  void Out(TObject *);
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

  bool operator==(MkFill &fill);
  bool operator!=(MkFill &fill);
  MkFill& operator=(MkFill &fill);
};
//---------------------------------------------------------------------------
class MkFills {
protected:
    MkFill *FFill;
    int FSize;//Actual size of entities
    int FSizeOfArray;
#ifdef __BCPLUSPLUS__
    TColor Color;
#endif
public:
    MkFills(int size,MkFill *ent);
    MkFills(int size);
    MkFills(){FSizeOfArray = FSize = 0;FFill = NULL;}
     ~MkFills();
    virtual void Initialize(int size);
    virtual void Initialize(int size,MkFill *);
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    bool Add(MkFill &fill);  // change of size of fill
    bool Add(int index,MkFill &fill);
    bool Delete(MkFill &fill);  // change of size of fill
    bool Delete(int index);
    int Grow(int Delta);            // change of size of array
    int Shrink(int Delta);          // change of size of array
    bool Clear();
#ifdef __BCPLUSPLUS__
    TColor GetColor(){return Color;};
    void SetColor(TColor c){Color = c;}
    void Out(TObject *);
#else
#endif

    void Out(char *fname);

    virtual MkFill & operator[](int);
    MkFills & operator=(MkFills &fills);
    bool operator==(MkFills &fills);

#ifdef __BCPLUSPLUS__
  void  Draw(TObject *);
#endif
    
#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};
//---------------------------------------------------------------------------
extern MkFill NullFill;
//---------------------------------------------------------------------------
#endif
