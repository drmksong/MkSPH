//---------------------------------------------------------------------------
#ifndef MkExcavStepH
#define MkExcavStepH
#include "MkEntity.h"
//---------------------------------------------------------------------------
enum MkExcavType {etNone=0, etInstall, etUninstall, etCut, etFill};
enum MkEntityType {ettNone=0, ettPile, ettStrut, ettAnchor, ettBolt, ettSlab, ettPanel, ettWale, ettCut, ettFill};

struct AStep {
  int          StepNum;
  MkExcavType  ExcavType;
  MkEntityType EntityType;
  MkSide       Side;
  int          EntityNum;
  float        Depth;  // valid for cut and fill only, any other support will refer the KoSupports.
  char         Script[256];
public: // constructor

#ifdef __BCPLUSPLUS__
  AStep(){StepNum=-1;ExcavType=MkExcavType(0);EntityType=ettNone;Side=mkNone;EntityNum=-1;Depth=0;memset(Script,'\0',255);}
#else
  AStep(){StepNum=-1;ExcavType=MkExcavType(0);EntityType=ettNone;Side=mkNone;EntityNum=-1;Depth=0;memset(Script,'\0',255);}
#endif

  AStep(int n, MkExcavType et, MkEntityType ett, int en)
    {
      StepNum = n;
      ExcavType = et;
      EntityType = ett;
      Side = mkNone;
      EntityNum = en;
      Depth = 0;
      UnParse();
    }

  AStep(int n, MkExcavType et, MkEntityType ett, int en, float depth)
    {
      StepNum = n;
      ExcavType = et;
      EntityType = ett;
      Side = mkNone;
      EntityNum = en;
      Depth = depth;
      UnParse();
    }
public: // setting function
  void SetScript(char *str){strcpy(Script,str);}
  void SetStepNum(int n){StepNum = n;}
  void SetExcavType(MkExcavType et){ExcavType = et;}
  void SetSide(MkSide side){Side = side;}
  void SetEntityType(MkEntityType ett){EntityType = ett;}
  void SetEntityNum(int en){EntityNum = en;}
  void SetDepth(float depth){Depth = depth;}
public:
  bool Parse();
  bool UnParse();
  char * Out(){UnParse();return Script;}
  bool In(char *str){strcpy(Script,str);return Parse();}
  AStep &operator=(AStep &step)
    {
      StepNum = step.StepNum;
      ExcavType = step.ExcavType;
      EntityType=step.EntityType;
      EntityNum=step.EntityNum;
      Depth=step.Depth;
      strcpy(Script,step.Script);
      return *this;
    }
  bool operator==(AStep &step)
    {
      return StepNum==step.StepNum&&
             ExcavType==step.ExcavType&&     
             EntityType==step.EntityType&&
             EntityNum==step.EntityNum&&
             fabs(Depth-step.Depth)<EPS&&
             !strcmp(Script,step.Script);
    }
};

class MkSteps {
public:
  AStep *Step;
  int Size;
public: // constructor
  MkSteps(){Step = NULL;Size = 0;}
  ~MkSteps();
  MkSteps(AStep *steps, int n);
  void Initialize(int n);
  void Initialize(AStep *steps, int n);
  void Clear(){if(Step) {delete Step;Step=NULL;} Size=0;}
  int GetSize(){return Size;}
  bool Sort();
  bool Add(AStep &step);
  bool Del(AStep &step);
  AStep &operator[](int n){return Step[n];}
  bool operator==(MkSteps &);
  bool operator!=(MkSteps &);
  MkSteps &operator=(MkSteps &);
};

class MkExcavStep {
protected:
  float DepthUnderDan;
  MkSteps Steps;
  MkSteps OneStep;
public: // constructors
  MkExcavStep(){}
  ~MkExcavStep(){}
public: // key functions
  bool Install(int step, MkEntityType ett, int en) {return Add(step,etInstall,ett,en);}
  bool Uninstall(int step, MkEntityType ett, int en) {return Add(step,etUninstall,ett,en);}
  bool Install(int step, MkEntityType ett, int en, float depth) {return Add(step,etInstall,ett,en,depth);}
  bool Uninstall(int step, MkEntityType ett, int en, float depth) {return Add(step,etUninstall,ett,en,depth);}
  bool Cut(int step, MkEntityType ett, int en, float depth) {return Add(step,etCut,ett,en,depth);}
  bool Fill(int step, MkEntityType ett, int en, float depth) {return Add(step,etFill,ett,en,depth);}
public: // setting
  void SetDepthUnderDan(float dud){DepthUnderDan=dud;}
  void SetSteps(MkSteps &steps){Steps = steps;}
  float GetDepthUnderDan(){return DepthUnderDan;}
  MkSteps &GetSteps(){return Steps;}
  int GetMaxExcavStep();
public: // manipulating member
  bool Sort(){return Steps.Sort();}
  bool Add(int step, MkExcavType exca_t, MkEntityType ett, int en)
    {
      AStep astep(step,exca_t,ett,en); 
      return Steps.Add(astep);
    }
  bool Del(int step, MkExcavType exca_t, MkEntityType ett, int en)
    {
      AStep astep(step,exca_t,ett,en);
      return Steps.Del(astep);
    }
  bool Add(int step, MkExcavType exca_t, MkEntityType ett, int en, float depth)
    {
      AStep astep(step,exca_t,ett,en,depth);
      return Steps.Add(astep);
    }
  bool Del(int step, MkExcavType exca_t, MkEntityType ett, int en, float depth)
    {
      AStep astep(step,exca_t,ett,en,depth);
      return Steps.Del(astep);
    }
  void Clear(){Steps.Clear();}
  bool Out(char *fname);
  bool In(char *fname);

  bool Out(FILE *fp);
  bool In(FILE *fp);

#ifdef __BCPLUSPLUS__
  bool Out(TObject *);
  bool In(TObject *);
#endif

  MkSteps &operator()(int n);
  bool operator==(MkExcavStep &);
  bool operator!=(MkExcavStep &);
  MkExcavStep &operator=(MkExcavStep &);
};
//---------------------------------------------------------------------------
#endif
