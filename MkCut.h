//---------------------------------------------------------------------------
#ifndef MkCutH
#define MkCutH
#include "MkEntity.h"
#include "MkWall.h"
#include "MkProfile.h"
//---------------------------------------------------------------------------
class MkCut : public MkEntity {
protected:
  MkWall *Wall[2];
//  MkPolygon CutLine;
  MkProfile Profile;
public:
  MkCut();
  MkCut(int );
  ~MkCut(){}
  void SetWall(MkWall *wall){Wall[0] = wall;}
  void SetWall(MkWall *wall, int n){if (n<0||n>2) return;Wall[n]=wall;}
  void SetProfile(MkProfile &prof){Profile = prof;}
//  void SetCutLine(MkPolygon &cutline){CutLine = cutline;}
public:
  MkWall *GetWall(int n){if (0<=n && n < 2) return Wall[n];else return NULL;}
  MkProfile &GetProfile(){return Profile;}
  float GetDepth(void){return Depth;}
  float GetDepth(float x);
//  MkPolygon &GetCutLine(){return CutLine;}
#ifdef __BCPLUSPLUS__
  AnsiString ClassName(){return AnsiString("MkCut");};
  bool UpdateFrom(TStringGrid*){return false;}
  bool UpdateTo(TStringGrid*){return false;}
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

  bool operator==(MkCut &);
  bool operator!=(MkCut &);
  MkCut& operator=(MkCut& cut);
};

class MkCuts {
protected:
    MkCut *FCut;
    int FSize;//Actual size of entities
    int FSizeOfArray;
#ifdef __BCPLUSPLUS__
    TColor Color;
#endif
public:
    MkCuts(int size,MkCut *ent);
    MkCuts(int size);
    MkCuts(){FSizeOfArray = FSize = 0;FCut = NULL;}
     ~MkCuts();
    virtual void Initialize(int size);
    virtual void Initialize(int size,MkCut *);
    int GetSize(){return FSize;};
    int GetNumber(){return FSize;};
    bool Add(MkCut &cut);  // change of size of cut
    bool Add(int index,MkCut &cut);
    bool Delete(MkCut &cut);  // change of size of cut
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

#ifdef __BCPLUSPLUS__
  void  Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif

    virtual MkCut & operator[](int);
    MkCuts & operator=(MkCuts &cuts);
    bool operator==(MkCuts &cuts);
};
//---------------------------------------------------------------------------
extern MkCut  NullCut;
//---------------------------------------------------------------------------
#endif
