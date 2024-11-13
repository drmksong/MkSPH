//---------------------------------------------------------------------------
#pragma hdrstop
#include "MkExcavStep.h"

#ifdef __BCPLUSPLUS__
#include <vcl.h>
#include "AdvMemo.hpp"
#endif

//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif
//---------------------------------------------------------------------------
/*
   step 0 cut 0
   step 0 install strut 0
   step 1 cut 1
   step 1 install strut 1
   step 1 install anchor 0
   step 2 fill 0
   step 3 uninstall anchor 0
*/

bool AStep::Parse()
{
  char str[256],s[256];
  int n,npile=0;
  float depth;

  MkDebug("Parse::");
  MkDebug(Script);

  strcpy(str,Script);
  ExtractStr(s,0,str);

  MkDebug(" 0 ::");MkDebug(s);MkDebug("\n");

  ToLower(s);

  if(strcmp(s,"step")) return false;

  ExtractStr(s,1,str);
  MkDebug(" 1 ::");MkDebug(s);MkDebug("\n");
  sscanf(s,"%d",&StepNum);

  ExtractStr(s,2,str);
  MkDebug(" 2 ::");MkDebug(s);MkDebug("\n");
  ToLower(s);

  if(!strcmp(s,"install")) ExcavType = etInstall;
  else if(!strcmp(s,"uninstall")) ExcavType = etUninstall;
  else if(!strcmp(s,"cut")) ExcavType = etCut;
  else if(!strcmp(s,"fill")) ExcavType = etFill;

  if(ExcavType == etInstall || ExcavType == etUninstall) {
    ExtractStr(s,3,str);
    MkDebug(" 3 ::");MkDebug(s);MkDebug("\n");
    ToLower(s);

    if(!strcmp(s,"pile")) {EntityType=ettPile;npile++;}
    else if(!strcmp(s,"strut")) EntityType=ettStrut;
    else if(!strcmp(s,"anchor")) EntityType=ettAnchor;
    else if(!strcmp(s,"bolt")) EntityType=ettBolt;
    else if(!strcmp(s,"slab")) EntityType=ettSlab;
    else if(!strcmp(s,"panel")) EntityType=ettPanel;
    else if(!strcmp(s,"wale")) EntityType=ettWale;
    ExtractStr(s,4,str);
    MkDebug(" 4 ::");MkDebug(s);MkDebug("\n");
    sscanf(s,"%d",&n);
    EntityNum=n;
    ExtractStr(s,5,str);
    if(strlen(s)) {
      sscanf(s,"%f",&depth);
      Depth = depth;
    }
    else Depth=0;
  }
  else if(ExcavType == etCut) {
    ExtractStr(s,3,str);
    MkDebug(" 3 ::");MkDebug(s);MkDebug("\n");
    sscanf(s,"%d",&n);
    EntityType = ettCut;
    EntityNum = n;

    ExtractStr(s,4,str);
    if(strlen(s)) {
      sscanf(s,"%f",&depth);
      Depth = depth;
    }
    else Depth=0;
  }

  else if(ExcavType == etFill) {
    ExtractStr(s,3,str);
    MkDebug(" 3 ::");MkDebug(s);MkDebug("\n");
    sscanf(s,"%d",&n);
    EntityType = ettFill;
    EntityNum = n;

    ExtractStr(s,4,str);
    if(strlen(s)) {
      sscanf(s,"%f",&depth);
      Depth = depth;
    }
    else Depth=0;
  }
  MkDebug("Out::");MkDebug(Out());MkDebug("\n");
  return true;
}

bool AStep::UnParse()
{
  char str[256],s[256];
  memset(str,'\0',255);
  memset(s,'\0',255);

  sprintf(str,"step %d ",StepNum);
  switch(ExcavType) {
    case etInstall:
      strcat(str,"Install ");break;
    case etUninstall:
      strcat(str,"Uninstall ");break;
    case etCut:
      strcat(str,"Cut ");break;
    case etFill:
      strcat(str,"Fill ");break;
    default:
      return false;
  }

  switch(EntityType) {
    case ettPile:
      strcat(str,"Pile ");break;
    case ettStrut:
      strcat(str,"Strut ");break;
    case ettAnchor:
      strcat(str,"Anchor ");break;
    case ettBolt:
      strcat(str,"Bolt ");break;
    case ettSlab:
      strcat(str,"Slab ");break;
    case ettPanel:
      strcat(str,"Panel ");break;
    case ettWale:
      strcat(str,"Wale ");break;
    case ettCut:
      /*strcat(str,"Cut ");*/break;
    case ettFill:
      /*strcat(str,"Fill ");*/break;
    default:
      return false;
  }

  sprintf(s," %d",EntityNum);
  strcat(str,s);
  if(Depth>EPS) sprintf(s," %.5g\n",Depth);
  else strcpy(s,"\n");
  strcat(str,s);

  strcpy(Script,str);
  return true;
}
//---------------------------------------------------------------------------
MkSteps::MkSteps(AStep *steps, int size)
{
  int i;
    if (size < 0) {
      MkDebug("::MkSteps - Initialize(int size)");
      return;
    }
    if (Size == size) return;

    Size = size;
    
    if (Size == 0) {
       if (Step!=NULL) delete[] Step;
       Step = NULL;
       return;
    }

    if (Step!=NULL) delete[] Step;
    Step = new AStep[Size];

    for(i=0;i<Size;i++) Step[i] = steps[i];
}

MkSteps::~MkSteps()
{
  if(Step) delete[] Step;
  Step = NULL;
  Size = 0;
}

void MkSteps::Initialize(int size)
{
  int i;
  AStep astep;
  if (size < 0) {
    MkDebug("::MkSteps - Initialize(int size)");
    return;
  }
  if (Size == size) return;

  Size = size;

  if (Size == 0) {
     if (Step!=NULL) delete[] Step;
     Step = NULL;
     return;
  }

  if (Step!=NULL) delete[] Step;
  Step = new AStep[Size];

  for(i=0;i<Size;i++) Step[i] = astep;
}

void MkSteps::Initialize(AStep *steps, int size)
{
  int i;
    if (size < 0) {
      MkDebug("::MkSteps - Initialize(int size)");
      return;
    }
    if (Size == size) {
      for(i=0;i<Size;i++) Step[i] = steps[i];
      return;
    }

    Size = size;
    
    if (Size == 0) {
       if (Step!=NULL) delete[] Step;
       Step = NULL;
       return;
    }

    if (Step!=NULL) delete[] Step;
    Step = new AStep[Size];

    for(i=0;i<Size;i++) Step[i] = steps[i];
}

bool MkSteps::Sort()
{
  bool flag;
  AStep step;
  do {
    flag = false;
    for(int i=0;i<Size-1;i++) {
      if(Step[i].StepNum>Step[i+1].StepNum) {
        step = Step[i+1];
        Step[i+1] = Step[i];
        Step[i] = step;
        flag = true;
      }
    }
  } while(flag);
  return true;
}

bool MkSteps::Add(AStep &astep)
{
  AStep *step;
  int size;
  int i,index=0;

  if (Size==0) {
    if(!(Step = new AStep[1])) return false;
    Step[0] = astep;
    Size = 1;
    return true;
  }

  size = Size+1;
  step = new AStep[size];
  if(!step) return false;

  for (i=1;i<Size;i++)
    if(Step[i-1].StepNum <= astep.StepNum && astep.StepNum < Step[i].StepNum) {
      index = i;
      break;
    }
  if(astep.StepNum < Step[0].StepNum) index = 0;
  else if (Step[Size-1].StepNum <= astep.StepNum) index = Size;

  for (i=0;i<size;i++)
    step[i] = i<index? Step[i] : (i==index)? astep : Step[i-1];

  if(Step) delete[] Step;
  Step = step;
  Size = size;

  return true;
}

bool MkSteps::Del(AStep &astep)
{
  AStep *step;
  int size,i,j;
  int index=-1;

  for (i=0;i<Size;i++) if(Step[i] == astep) index = i;
  if(index==-1) return false;

  size = Size-1;
  step = new AStep[size];
  if(!step) return false;

  for(i=0;i<Size;i++) {
    j=i<=index?i:i-1;
    step[j] = Step[i];
  }

  if(Step) delete[] Step;
  Step = step;
  Size = size;
  return true;
}

bool MkSteps::operator==(MkSteps &steps)
{
  int i;
  bool flag=true;
  for (i=0;i<Size;i++) {
    flag = flag && Step[i] == steps[i];
    if (!flag) return flag;
  }
  return flag;
}

bool MkSteps::operator!=(MkSteps &steps)
{
  return !(*this==steps);
}

// this routine seems to have some bugs...
// try to find possible errors.
MkSteps &MkSteps::operator=(MkSteps &steps)
{
  Initialize(steps.Step,steps.Size);
  return *this;
}
//---------------------------------------------------------------------------
bool MkExcavStep::Out(char *fname)
{
  int i;
  FILE *fp;
  fp = fopen(fname,"a");
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
  fprintf(fp,"\n");
  fprintf(fp,"                           <Information of Excavation Step>\n");
  fprintf(fp,"\n");
  fprintf(fp,"%d : number of steps\n",Steps.Size);
  for(i=0;i<Steps.Size;i++) {
    fprintf(fp,Steps[i].Out());
  }
  fprintf(fp,"DepthUnderDan %.5g\n",DepthUnderDan);
  fclose(fp);
  return true;
}

bool MkExcavStep::In(char *fname)
{
  int i,size;
  char str[256];
  AStep step;
  FILE *fp;

  fp = fopen(fname,"r");
  if(!fp) {
    MkDebug("Strange :: the file is not found.");
    return false;
  }

  Steps.Clear();
  fgets(str,255,fp);
  fgets(str,255,fp);
  fgets(str,255,fp);
  fgets(str,255,fp);
  sscanf(str,"%d",&size);

  for(i=0;i<size;i++) {
    fgets(str,255,fp);
    step.In(str);
    MkDebug(step.Out());
    Steps.Add(step);
  }
  fclose(fp);
  return true;
}

//---------------------------------------------------------------------------
int MkExcavStep::GetMaxExcavStep()
{
  int m=0;
  for(int i=0;i<Steps.GetSize();i++) m=(Steps[i].StepNum > m)?Steps[i].StepNum:m;
  return m;
}

bool MkExcavStep::Out(FILE *fp)
{
  int i;

  if(!fp) {
    MkDebug("Strange :: the file is not found.");
    return false;
  }
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
  fprintf(fp,"\n");
  fprintf(fp,"                           <Information of Excavation Step>\n");
  fprintf(fp,"\n");
  fprintf(fp,"%d : number of steps\n",Steps.Size);
  for(i=0;i<Steps.Size;i++) {
    fprintf(fp,Steps[i].Out());
  }
  fprintf(fp,"DepthUnderDan %.5g",DepthUnderDan);
  return true;
}

bool MkExcavStep::In(FILE *fp)
{
  int i,size;
  char str[256],dummy[256];
  AStep step;

  if(!fp) {
    MkDebug("Strange :: the file is not found.");
    return false;
  }

  Steps.Clear();
  fgets(str,255,fp);
  fgets(str,255,fp);
  fgets(str,255,fp);
  fgets(str,255,fp);
  sscanf(str,"%d",&size);

  for(i=0;i<size;i++) {
    fgets(str,255,fp);
    step.In(str);
    MkDebug(step.Out());
    Steps.Add(step);
  }
  fgets(str,255,fp);
  if(strlen(str)<0) {
    DepthUnderDan = 0.5;
    return false;
  }
  else sscanf(str,"%s %f",dummy,&DepthUnderDan);

  return true;
}

#ifdef __BCPLUSPLUS__
bool MkExcavStep::Out(TObject *Sender)
{
  int i;
  AnsiString str;

  TAdvMemo *memo=dynamic_cast<TAdvMemo*>(Sender);
  if(!memo) {
    MkDebug("Strange :: this object is not a memo component.");
    return false;
  }
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
  memo->Lines->Add(" ");
  memo->Lines->Add(" <Information of Excavation Step>");
  memo->Lines->Add(" ");
  memo->Lines->Add(AnsiString(Steps.Size)+" : number of steps");
  for(i=0;i<Steps.Size;i++) {
    str = Steps[i].Out();
    str = str.SubString(1,str.Length()-1);
    memo->Lines->Add(str);
  }
  return true;
}

bool MkExcavStep::In(TObject *Sender)
{
  int i,n,size=0;
  char str[256];
  AStep step;

  TAdvMemo *memo=dynamic_cast<TAdvMemo*>(Sender);
  if(!memo) {
    MkDebug("Strange :: this object is not a memo component.");
    return false;
  }

  Steps.Clear();

  n=0;
  while(size==0 && n<100) {
    if(!memo->Lines->Strings[n].Length()){n++;continue;}
    sscanf(memo->Lines->Strings[n].c_str(),"%d",&size);
    n++;
  }
  if(n>=100) return false;

  for(i=0;i<size&&n<100;n++) {
    if(!memo->Lines->Strings[n].Length()) {continue;}
    sscanf(memo->Lines->Strings[n].c_str(),"%s",str);
    if(AnsiString("Step").UpperCase()!=AnsiString(str).UpperCase()) {continue;}
    step.In(memo->Lines->Strings[n].c_str());
    MkDebug(step.Out());
    Steps.Add(step);
    i++;
  }
  return true;
}
#endif

MkSteps &MkExcavStep::operator()(int n)
{
  int i;
  OneStep.Clear();
  for (i=0;i<Steps.Size;i++)
    if (Steps[i].StepNum<=n)
      OneStep.Add(Steps[i]);
  return OneStep;
}

bool MkExcavStep::operator==(MkExcavStep &es)
{
  return Steps == es.Steps;
}

bool MkExcavStep::operator!=(MkExcavStep &es)
{
  return !(Steps == es.Steps);
}

MkExcavStep &MkExcavStep::operator=(MkExcavStep &es)
{
  Steps = es.Steps;
  DepthUnderDan = es.DepthUnderDan;
  return *this;
}

