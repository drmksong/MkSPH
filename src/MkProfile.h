//---------------------------------------------------------------------------
#ifndef MkProfileH
#define MkProfileH

#include <stdio.h>
#include "MkPolygon.hpp"
//---------------------------------------------------------------------------
enum MkProfType {pfNone, pfLayer, pfCut, pfFill};

class MkProfile {
protected:
  MkProfType ProfType;
  int ProfNum;
  MkPolygon Polygon;
  char *Script;

public:
  MkProfile(){ProfType = pfNone; ProfNum=0; Script=NULL;}
  MkProfile(int ){ProfType = pfNone; ProfNum=0; Script=NULL;}
  MkProfile(MkProfType pf, MkPolygon &poly);
  ~MkProfile();

public:
  void Initialize(int n){Polygon.Initialize(n);}

  void SetProfile(MkPolygon &poly){Polygon = poly;MakeScript();}
  void SetProfType(MkProfType &pf){ProfType = pf;MakeScript();}
  void SetProfNum(int n){ProfNum = n;MakeScript();}
  void SetBound(float lb, float rb);
  void SetDivision(MkPoints &pnts);

  MkPolygon &GetProfile(){return Polygon;}
  MkProfType &GetProfType(){return ProfType;}
  int GetProfNum(){return ProfNum;}
  int GetSize(){return Polygon.GetSize();}
public:
  void Clear(){ClearScript(); ProfType=pfNone;ProfNum=0;Polygon.Clear();}
  void ClearScript();
  bool MakeScript();
  bool ParseScript(char *str);
  char *GetScript(){return Script;}

public:
  MkPoint &operator[](int n) {return Polygon[n];}

  bool operator==(MkProfile &pf);
  bool operator!=(MkProfile &pf);
  MkProfile &operator=(MkProfile &pf);
};

class MkProfiles {
protected:
  char FileName[512];
  MkProfile *FProfile;
  int FSize;//Actual size of profs
  int FSizeOfArray;
public:
  MkProfiles(int size,MkProfile *prof);
  MkProfiles(int size);
  MkProfiles(){FSizeOfArray = FSize = 0;FProfile = NULL;}
  ~MkProfiles();
  virtual void Initialize(int size);
  virtual void Initialize(int size,MkProfile *);
  int GetSize(){return FSize;};
  int GetNumber(){return FSize;};
  bool Add(MkProfile &prof);  // change of size of prof
  bool Add(MkProfiles &prof) {
       bool flag=true;
       MkDebug("MkProfiles::Add(Profiles) is called.\n");
       for (int i=0;i<prof.GetSize();i++)
         flag = Add(prof[i]) && flag;
       return flag;
  }  // change of size of prof
  bool Add(int index,MkProfile &prof);
  bool Delete(MkProfile &prof);  // change of size of prof
  bool Delete(int index);
  int Grow(int Delta);            // change of size of array
  int Shrink(int Delta);          // change of size of array
  bool Clear();
  void Out();
  virtual MkProfile & operator[](int);
  MkProfiles & operator=(MkProfiles &profs);
  bool operator==(MkProfiles &profs);
public:
  void SetFileName(char *fname){strcpy(FileName,fname);}
  char *GetFileName(){return FileName;}
public:
  bool Open();
  bool Save();
  bool Open(char *fname){SetFileName(fname);return Open();}
  bool Save(char *fname){SetFileName(fname);return Save();}
  bool SyncDivision();
};
extern MkProfile NullProfile;
#endif

