//---------------------------------------------------------------------------
#include "MkLayer.h"
#ifdef __BCPLUSPLUS__
#include <vcl.h>
#endif
#pragma hdrstop
//---------------------------------------------------------------------------
#ifdef __BCPLUSPLUS__
#pragma package(smart_init)
#endif
//---------------------------------------------------------------------------
MkLayer NullLayer(0);
MkLayers NullLayers(0);
MkLayer::MkLayer()
{
  Number=0;
#ifdef __BCPLUSPLUS__
  Name=" ";
#else
  memset(Name,'\0',255);
#endif

  Side = mkLeft;
  WetUnitWeight[0]=0;
  SubUnitWeight[0]=0;
  Cohesion[0]=0;
  Friction[0]=0;
  HorSubReact[0]=0;
  VerSubReact[0]=0;
  WetUnitWeight[1]=0;
  SubUnitWeight[1]=0;
  Cohesion[1]=0;
  Friction[1]=0;
  HorSubReact[1]=0;
  VerSubReact[1]=0;
  Burden=0;
  K0 = 1;
  Ka = 0;
  Kp = 0;
  Bearing=0;
  className = "MkLayer";
}

MkLayer::MkLayer(int n)
{
  Number=0;
#ifdef __BCPLUSPLUS__
  Name=" ";
#else
  memset(Name,'\0',255);
#endif

  Side            =mkLeft;
  WetUnitWeight[0]=0;
  SubUnitWeight[0]=0;
  Cohesion[0]     =0;
  Friction[0]     =0;
  HorSubReact[0]  =0;
  VerSubReact[0]  =0;
  WetUnitWeight[1]=0;
  SubUnitWeight[1]=0;
  Cohesion[1]     =0;
  Friction[1]     =0;
  HorSubReact[1]  =0;
  VerSubReact[1]  =0;
  Burden          =0;
  K0              =1;
  Ka              =0;
  Kp              =0;
  Bearing         =0;
  className = "MkLayer";
}
#ifdef __BCPLUSPLUS__
bool MkLayer::UpdateFrom()
{
  if(!Grid) return false;

  Number          =Grid->Cells[1][0].ToInt();
  Name            =Grid->Cells[1][1];
  WetUnitWeight[0]=Grid->Cells[1][2].ToDouble();
  SubUnitWeight[0]=Grid->Cells[1][3].ToDouble();
  Cohesion[0]     =Grid->Cells[1][4].ToDouble();
  Friction[0]     =Grid->Cells[1][5].ToDouble();
  HorSubReact[0]  =Grid->Cells[1][6].ToDouble();
  VerSubReact[0]  =Grid->Cells[1][7].ToDouble();
  WetUnitWeight[1]=Grid->Cells[1][8].ToDouble();
  SubUnitWeight[1]=Grid->Cells[1][9].ToDouble();
  Cohesion[1]     =Grid->Cells[1][10].ToDouble();
  Friction[1]     =Grid->Cells[1][11].ToDouble();
  HorSubReact[1]  =Grid->Cells[1][12].ToDouble();
  VerSubReact[1]  =Grid->Cells[1][13].ToDouble();

  return true;
}

bool MkLayer::UpdateTo()
{
  if(!Grid) return false;

  for(int i=0;i<20;i++) {Grid->Cells[0][i] = "";Grid->Cells[1][i]="";}

  Grid->Cells[0][0] =  "Number";
  Grid->Cells[0][1] =  "Name";
  Grid->Cells[0][2] =  "WetUnitWeight[0]";
  Grid->Cells[0][3] =  "SubUnitWeight[0]";
  Grid->Cells[0][4] =  "Cohesion[0]";
  Grid->Cells[0][5] =  "Friction[0]";
  Grid->Cells[0][6] =  "HorSubReact[0]";
  Grid->Cells[0][7] =  "VerSubReact[0]";
  Grid->Cells[0][8] =  "WetUnitWeight[1]";
  Grid->Cells[0][9] =  "SubUnitWeight[1]";
  Grid->Cells[0][10] = "Cohesion[1]";
  Grid->Cells[0][11] = "Friction[1]";
  Grid->Cells[0][12] = "HorSubReact[1]";
  Grid->Cells[0][13] = "VerSubReact[1]";

  Grid->Cells[1][0] =  Number;
  Grid->Cells[1][1] =  Name;
  Grid->Cells[1][2] =  WetUnitWeight[0];
  Grid->Cells[1][3] =  SubUnitWeight[0];
  Grid->Cells[1][4] =  Cohesion[0];
  Grid->Cells[1][5] =  Friction[0];
  Grid->Cells[1][6] =  HorSubReact[0];
  Grid->Cells[1][7] =  VerSubReact[0];
  Grid->Cells[1][8] =  WetUnitWeight[1];
  Grid->Cells[1][9] =  SubUnitWeight[1];
  Grid->Cells[1][10] = Cohesion[1];
  Grid->Cells[1][11] = Friction[1];
  Grid->Cells[1][12] = HorSubReact[1];
  Grid->Cells[1][13] = VerSubReact[1];

  return true;
}
#endif

float MkLayer::GetActivPress(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float y,ymax,p,t;
  float f,c;

  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;
  f = Friction[0]+(Friction[1]-Friction[0])*y/Rect.GetHeight();
  c = Cohesion[0]+(Cohesion[1]-Cohesion[0])*y/Rect.GetHeight();

  if(Ka>0&&Kp>0) {
    p = (WetUnitWeight[0]*y+(WetUnitWeight[1]-WetUnitWeight[0])*y*y/2/Rect.GetHeight()+Burden);
    assert(Ka>0);
    return p*Ka-2*c*sqrt(Ka);
  }
  else {
    t = tan((45-f/2)*M_PI/180.0);
    p = (WetUnitWeight[0]*y+(WetUnitWeight[1]-WetUnitWeight[0])*y*y/2/Rect.GetHeight()+Burden)*K0;
    return p*t*t-2*c*t;
  }
}

float MkLayer::GetPassivPress(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float y,ymax,p,t;
  float f,c;

  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;

  f = Friction[0]+(Friction[1]-Friction[0])*y/Rect.GetHeight();
  c = Cohesion[0]+(Cohesion[1]-Cohesion[0])*y/Rect.GetHeight();

  if(Ka>0&&Kp>0) {
    p = (WetUnitWeight[0]*y+(WetUnitWeight[1]-WetUnitWeight[0])*y*y/2/Rect.GetHeight()+Burden);
    assert(Kp>0);
    return p*Kp+2*c*sqrt(Kp);
  }
  else {
    t = tan((45+f/2)*M_PI/180.0);
    p = (WetUnitWeight[0]*y+(WetUnitWeight[1]-WetUnitWeight[0])*y*y/2/Rect.GetHeight()+Burden)*K0;
    return p*t*t+2*c*t;
  }
}

float MkLayer::GetStopPress(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float ymax,y;
  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;
  return (WetUnitWeight[0]*y+(WetUnitWeight[1]-WetUnitWeight[0])*y*y/2/Rect.GetHeight()+Burden)*K0;
}

float MkLayer::GetWetUnitWeight(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float ymax,y,value;
  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;
  value = (WetUnitWeight[1]-WetUnitWeight[0])*y/Rect.GetHeight()+WetUnitWeight[0];
  return value;
}

float MkLayer::GetCohesion(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float ymax,y;
  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;
  return (Cohesion[1]-Cohesion[0])*y/Rect.GetHeight()+Cohesion[0];
}

float MkLayer::GetFriction(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float ymax,y;
  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;
  return (Friction[1]-Friction[0])*y/Rect.GetHeight()+Friction[0];
}

float MkLayer::GetHorSubReact(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float ymax,y;
  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;
  return (HorSubReact[1]-HorSubReact[0])*y/Rect.GetHeight()+HorSubReact[0];
}

float MkLayer::GetVerSubReact(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  float ymax,y;
  ymax = max(Rect.GetOrigin().Y,Rect.GetOrigin().Y+Rect.GetHeight());
  y = ymax - pnt.Y;
  return (VerSubReact[1]-VerSubReact[0])*y/Rect.GetHeight()+VerSubReact[0];
}

float MkLayer::GetK0(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  return K0;
}

float MkLayer::GetKa(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  return Ka;
}

float MkLayer::GetKp(MkPoint &pnt)
{
  if(!IsIn(pnt)) return -1;
  return Kp;
}

float MkLayer::GetWetUnitWeight(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  float uw,lw;
  uw = GetWetUnitWeight(l[1]);
  lw = GetWetUnitWeight(l[0]);

  return (uw+lw)/2;
}

float MkLayer::GetCohesion(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  return (GetCohesion(l[1])+GetCohesion(l[0]))/2;
}

float MkLayer::GetFriction(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  return (GetFriction(l[1])+GetFriction(l[0]))/2;
}

float MkLayer::GetHorSubReact(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  return (GetHorSubReact(l[1])+GetHorSubReact(l[0]))*l.GetLength()/2;
}

float MkLayer::GetVerSubReact(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  return (GetVerSubReact(l[1])+GetVerSubReact(l[0]))*l.GetLength()/2;
}

float MkLayer::GetK0(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  return (GetK0(l[1])+GetK0(l[0]))*l.GetLength()/2;
}

float MkLayer::GetKa(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  return (GetKa(l[1])+GetKa(l[0]))*l.GetLength()/2;
}

float MkLayer::GetKp(MkLine &line)
{
  if(InLen(line)<0.0001) return -1;
  MkLine l=InLine(line);

  return (GetKp(l[1])+GetKp(l[0]))*l.GetLength()/2;
}

float MkLayer::InLen(MkLine &line)
{
  MkPoints pnts;
  Rect.GetCross(line,pnts);

  if(IsIn(line[0])&&IsIn(line[1])) return line.GetLength();
  if(!pnts.GetSize()) return 0;
  if(pnts.GetSize()==2) return MkLine(pnts[0],pnts[1]).GetLength();
  if(pnts.GetSize()==1 && IsIn(line[0])) return MkLine(line[0],pnts[0]).GetLength();
  if(pnts.GetSize()==1 && IsIn(line[1])) return MkLine(line[1],pnts[0]).GetLength();
  return 0;
}

MkLine MkLayer::InLine(MkLine &line)
{
  MkPoints pnts;
  Rect.GetCross(line,pnts);

  if(IsIn(line[0])&&IsIn(line[1])) return line;
  if(!pnts.GetSize()) return NullLine;
  if(pnts.GetSize()==2) return MkLine(pnts[0],pnts[1]);
  if(pnts.GetSize()==1 && IsIn(line[0])) return MkLine(line[0],pnts[0]);
  if(pnts.GetSize()==1 && IsIn(line[1])) return MkLine(line[1],pnts[0]);
  return NullLine;
}

void MkLayer::CalcCoeff()  // need to be update to consider variation in a layer
{
  K0 = 1-sin(Friction[0]*3.141592/180.0);
  Ka = (1+sin(Friction[0]*3.141592/180.0))/(1-sin(Friction[0]*3.141592/180.0));
  Kp = (1-sin(Friction[0]*3.141592/180.0))/(1+sin(Friction[0]*3.141592/180.0));
}

//  AnsiString layername;    MkEntity::Name
//  double depth,R_depth;    MkEntity::Depth
//  double Rt,Rsub,C,Ks,Pi;  WetUnitWeight, SubUnitWeight, Cohesion, HorSubReact, Friction
//  double bearing;          Bearing
//  MkFloat Profile;         Polygon

#ifdef __BCPLUSPLUS__
void MkLayer::Import(MkGlobalVar &globalvar, int sec,MkSide side,int lay)
{
  SetSide(side);
  Import(globalvar,sec,lay);
}
#endif

void MkLayer::Import(MkGlobalVar &globalvar, int sec,int lay)
{
#ifdef __BCPLUSPLUS__
  Name=(LPCTSTR)globalvar.layername[lay+1];
#else
  strncpy(Name,(LPCTSTR)globalvar.layername[lay+1]);
#endif
  Depth=((Side==mkLeft)? globalvar.layer_depth_L[sec+1][lay+1]:globalvar.layer_depth_R[sec+1][lay+1]);
  R_depth=((Side==mkLeft)? R_layer_depth_L[sec+1][lay+1]:R_layer_depth_R[sec+1][lay+1]);
  WetUnitWeight[0]=WetUnitWeight[1]=globalvar.Rt[sec+1][lay+1];
  SubUnitWeight[0]=SubUnitWeight[1]=globalvar.Rsub[sec+1][lay+1];
  Cohesion[0]=Cohesion[1]=globalvar.C[sec+1][lay+1];
  HorSubReact[0]=HorSubReact[1]=globalvar.Ks[sec+1][lay+1];
  Friction[0]=Friction[1]=globalvar.Pi[sec+1][lay+1];
  Bearing=::bearing[lay+1];
}

void MkLayer::Export(MkGlobalVar &globalvar, int sec,int lay)
{
#ifdef __BCPLUSPLUS__
  globalvar.layername[lay+1]=Name.c_str();
#else
  globalvar.layername[lay+1]=Name;
#endif
  ((Side==mkLeft)? globalvar.layer_depth_L[sec+1][lay+1]:globalvar.layer_depth_R[sec+1][lay+1])=Depth;
  ((Side==mkLeft)? R_layer_depth_L[sec+1][lay+1]:R_layer_depth_R[sec+1][lay+1])=R_depth;
  globalvar.Rt[sec+1][lay+1]=WetUnitWeight[0];
  globalvar.Rsub[sec+1][lay+1]=SubUnitWeight[0];
  globalvar.C[sec+1][lay+1]=Cohesion[0];
  globalvar.Ks[sec+1][lay+1]=HorSubReact[0];
  globalvar.Pi[sec+1][lay+1]=Friction[0];
  ::bearing[lay+1]=Bearing;
}

#ifdef __BCPLUSPLUS__
void MkLayer::Out(TObject *)
{

}
#endif

void MkLayer::Out(char *fname)
{
  FILE *fp;
  fp = fopen(fname,"a");
  if(!fp) {
    MkDebug(fname); MkDebug(" is not found, so fp is null and return false\n");
    return ;
  }
            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
//fprintf(fp,"  Layer  G.L.    rt     rsub Cohesion  Friction    Ks\n");
//fprintf(fp,"   No.   (m)   (t/m3)  (t/m3) (t/m2)   (degree)  (t/m3)\n");
  fprintf(fp,"  %3d Top:%5.2f %5.2f  %5.2f   %5.2f     %5.2f    %10.3f\n",Number,Rect.GetTop(),WetUnitWeight[0]*MPa2Tonf,SubUnitWeight[0]*MPa2Tonf,Cohesion[0]*MPa2Tonf,Friction[0],HorSubReact[0]*MPa2Tonf);
  fprintf(fp,"      Bot:%5.2f %5.2f  %5.2f   %5.2f     %5.2f    %10.3f\n",       Rect.GetBot(),WetUnitWeight[1]*MPa2Tonf,SubUnitWeight[1]*MPa2Tonf,Cohesion[1]*MPa2Tonf,Friction[1],HorSubReact[1]*MPa2Tonf);
  fprintf(fp,"\n");
  fclose(fp);
}

bool MkLayer::operator==(MkLayer& lay)
{
  bool flag;
  flag = MkEntity::operator==((MkEntity &)lay);
  flag = flag && Rect == lay.Rect;
  flag = flag && Profile == lay.Profile;
  flag = flag && WetUnitWeight[0] == lay.WetUnitWeight[0];
  flag = flag && SubUnitWeight[0] == lay.SubUnitWeight[0];
  flag = flag && Cohesion[0] == lay.Cohesion[0];
  flag = flag && Friction[0] == lay.Friction[0];
  flag = flag && HorSubReact[0] == lay.HorSubReact[0];
  flag = flag && VerSubReact[0] == lay.VerSubReact[0];
  flag = flag && WetUnitWeight[1] == lay.WetUnitWeight[1];
  flag = flag && SubUnitWeight[1] == lay.SubUnitWeight[1];
  flag = flag && Cohesion[1] == lay.Cohesion[1];
  flag = flag && Friction[1] == lay.Friction[1];
  flag = flag && HorSubReact[1] == lay.HorSubReact[1];
  flag = flag && VerSubReact[1] == lay.VerSubReact[1];
  flag = flag && Burden == lay.Burden;
  flag = flag && K0 == lay.K0;
  flag = flag && Ka == lay.Ka;
  flag = flag && Kp == lay.Kp;
  return flag;
}

bool MkLayer::operator!=(MkLayer& lay)
{
  return !operator==(lay);
}

MkLayer & MkLayer::operator=(MkLayer &lay)
{
  MkEntity::operator=((MkEntity &)lay);
  Rect = lay.Rect;
  Profile = lay.Profile;
  WetUnitWeight[0] = lay.WetUnitWeight[0];
  SubUnitWeight[0] = lay.SubUnitWeight[0];
  Cohesion[0] = lay.Cohesion[0];
  Friction[0] = lay.Friction[0];
  HorSubReact[0] = lay.HorSubReact[0];
  VerSubReact[0] = lay.VerSubReact[0];
  WetUnitWeight[1] = lay.WetUnitWeight[1];
  SubUnitWeight[1] = lay.SubUnitWeight[1];
  Cohesion[1] = lay.Cohesion[1];
  Friction[1] = lay.Friction[1];
  HorSubReact[1] = lay.HorSubReact[1];
  VerSubReact[1] = lay.VerSubReact[1];
  Burden = lay.Burden;
  K0 = lay.K0;
  Ka = lay.Ka;
  Kp = lay.Kp;
  return *this;
}

#ifdef __BCPLUSPLUS__
void  MkLayer::Draw(TObject *Sender)
{

}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void  MkLayer::Draw(MkPaint *pb)
{

}
#endif
//---------------------------------------------------------------------------
MkLayers::MkLayers(int size,MkLayer *layers)
{
    if (size < 0) {
      MkDebug("::MkLayers - MkLayers (int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSizeOfArray == 0) {
       FLayer= NULL;
       return;
    }

    FLayer= new MkLayer[FSizeOfArray];
    for (int i=0;i<FSizeOfArray;i++) (*this)[i] = layers[i];
}

MkLayers::MkLayers(int size)
{
    if (size < 0) {
      MkDebug("::MkLayers - MkLayers (int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSizeOfArray == 0) {
       FLayer= NULL;
       return;
    }

    FLayer= new MkLayer[FSizeOfArray];
}

MkLayers::~MkLayers()
{
   FSizeOfArray = FSize = 0;
   if (FLayer) {
      delete[] FLayer;
      FLayer = NULL;
   }
}

void MkLayers::SetupPress()
{
  float burden,xm,ym,xmin,ymin,xmax,ymax;
  for(int i=0;i<FSize;i++) {
    MkRect &rect = FLayer[i].GetRect();
    burden = 0;
    xmin = min(rect.GetOrigin().X,rect.GetOrigin().X+rect.GetWidth());
    xmax = max(rect.GetOrigin().X,rect.GetOrigin().X+rect.GetWidth());
    ymin = min(rect.GetOrigin().Y,rect.GetOrigin().Y+rect.GetHeight());
    ymax = max(rect.GetOrigin().Y,rect.GetOrigin().Y+rect.GetHeight());

    for(int j=0;j<FSize;j++) {
      MkRect &rec = FLayer[j].GetRect();
      xm = rec.GetOrigin().X+rec.GetWidth()/2.0;
      ym = rec.GetOrigin().Y+rec.GetHeight()/2.0;
      if((xmin<xm&&xm<xmax)&&(ymax<ym))
        burden+=FLayer[j].GetWetUnitWeight(rec.GetOrigin());// should be reimplemented.2004.08.09
    }
    FLayer[i].SetBurden(burden);
  }
}

void MkLayers::Initialize(int size)
{
    int i;
    if (size < 0) {
      MkDebug("::MkLayers - Initialize(int size)");
      return;
    }
    if (FSizeOfArray== size) return;

    MkLayer *layer = new MkLayer[size];
    if(layer) {
      for(i=0;i<min(size,FSizeOfArray);i++)
        layer[i] = FLayer[i];
    }

    FSizeOfArray= FSize=size;

    if (FSizeOfArray== 0) {
       if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
       FLayer = NULL;
       return;
    }

    if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
    FLayer = (layer!=NULL) ? layer : new MkLayer[FSizeOfArray];
    for (i=0;i<FSize;i++) FLayer[i].Number = i;
}

void MkLayers::Initialize(int size,MkLayer *layers)
{
    int i;
    if (size < 0 || layers == NULL) {
      MkDebug("::MkLayers - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSizeOfArray = size;
    FSize = size;

    if (FSizeOfArray== 0) {
       if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
       FLayer = NULL;
       return;
    }

    if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
    FLayer = new MkLayer[FSizeOfArray];
    for (i=0;i<FSizeOfArray;i++) FLayer[i] = layers[i];
    for (i=0;i<FSize;i++) FLayer[i].Number = i;
}

int MkLayers::Grow(int delta)
{
    int i;
    MkLayer *layer=NULL;

    if (!(layer = new MkLayer[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        layer[i] = FLayer[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        layer[i] = NullLayer;
    if (FLayer) {
       delete[] (MkLayer*)FLayer;
       FLayer = NULL;
    }
    FLayer = layer;
    FSizeOfArray = FSizeOfArray+delta;
    return FSizeOfArray;
}

int MkLayers::Shrink(int delta)
{
    int i;
    MkLayer *layer=NULL;

    if (!(layer = new MkLayer[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        layer[i] = FLayer[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        layer[i] = NullLayer;
    if (FLayer) {
       delete[] (MkLayer*)FLayer;
       FLayer = NULL;
    }
    FLayer = layer;
    FSizeOfArray = FSizeOfArray-delta;
    return FSizeOfArray;
}

bool MkLayers::Add(MkLayer &layer)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }
    
    FSize++;
    FLayer[FSize-1] = layer;
    return true;
}

bool MkLayers::Add(int index, MkLayer &layer)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) {
      Grow(FSize-FSizeOfArray+1);
      if(tmp==FSizeOfArray) return false;
    }

    for (int i=FSize-1;i>=index;i--)
      FLayer[i+1] = FLayer[i];
    FSize++;
    FLayer[index] = layer;
    return true;
}

bool MkLayers::Delete(MkLayer &layer)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FLayer[i] == layer) break;
    }
    if(i==FSize) return false;
    if(FLayer[i] == layer) {
      for (int j=i;j<FSize-1;j++)
        FLayer[j] = FLayer[j+1];
    }
    FSize--;
    FLayer[FSize] = NullLayer;
    return true;
}

bool MkLayers::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FLayer[j] = FLayer[j+1];

    FSize--;
    FLayer[FSize] = NullLayer;
    return true;
}

bool MkLayers::Clear()
{
   FSize = 0;
   if (FLayer) {
      delete[] FLayer;
      FLayer = NULL;
   }
   return true;
}
/*
bool MkLayers::Apply(MkCut &cut)
{
  return false;
}

bool MkLayers::Apply(MkFill &fill)
{
  return false;
}
*/

#ifdef __BCPLUSPLUS__
void MkLayers::Import(MkGlobalVar &globalvar, int sec)
{
  int i;
  Initialize(globalvar.layer_ea);
  for (i=0;i<globalvar.layer_ea && i<10;i++) {
    FLayer[i].SetSide(mkLeft);
    FLayer[i].Import(globalvar,sec,i);
  }
}

void MkLayers::Export(MkGlobalVar &globalvar, int sec)
{
  int i;
  for (i=0;i<FSize && i<10;i++)
    FLayer[i].Export(globalvar,sec,i);
}
#endif

MkLayer & MkLayers::operator[](int i)
{
    if (0<=i && i<FSize) return FLayer[i];
    else if(FSize<=i && i<FSizeOfArray) {
      FSize=i+1;
      return FLayer[i];
    }
    else if (FSizeOfArray <= i) {
      Grow(i-FSizeOfArray+5);
      return NullLayer;
    }
    else return NullLayer;
}

MkLayers & MkLayers::operator=(MkLayers &layers)
{
    int i;

    Clear();
    FSize = layers.FSize;
    FSizeOfArray = layers.FSizeOfArray;

    if (FSizeOfArray == 0) {
       FLayer = NULL;
       return *this;
    }
    this->FLayer = new MkLayer[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FLayer[i] = layers.FLayer[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FLayer[i] = NullLayer;


    return *this;
}

bool MkLayers::operator==(MkLayers &layers)
{
    int i;

    if (FSize != layers.FSize) return false;
    for (i=0;i<FSize;i++)
      if (this->FLayer[i] != layers.FLayer[i]) return false;

    return true;
}

#ifdef __BCPLUSPLUS__
void MkLayers::Out(TObject *Sender)
{

}
#endif

void MkLayers::Out(char *fname)
{
  FILE *fp;
  fp = fopen(fname,"a");
  if(!fp) {
    return;
  }
  if(FSize==0) return;

            //12345678901234567890123456789012345678901234567890123456789012345678901234567890
  fprintf(fp,"\n");
  fprintf(fp,"                           <Information of Layers>\n");
  fprintf(fp,"\n");
  fprintf(fp,"  Layer  G.L.  rt    rsub    Cohesion    Friction    Ks\n");
  fprintf(fp,"   No.   (m) (t/m3) (t/m3)    (t/m2)     (degree)  (t/m3)\n");
  fclose(fp);

  for(int i=0;i<FSize;i++)
    FLayer[i].Out(fname);
}

#ifdef __BCPLUSPLUS__
void MkLayers::Draw(TObject *Sender)
{
    TColor C;
    float Offset;
    if (FSize == 0) return;
    if (String(Sender->ClassName()) == String("MkPaintBox")) {
       MkPaintBox *pb=(MkPaintBox*)Sender;
       C = pb->Canvas->Pen->Color;
       pb->Canvas->Pen->Color = Color;
       for (int i = 0; i < FSize;i++) {
           Offset = pb->Offset(3);
           FLayer[i].Draw(pb);
       }
       pb->Canvas->Pen->Color = C;
    }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkLayers::Draw(MkPaint *pb)
{

}
#endif
//---------------------------------------------------------------------------
/*
MkLayer NullMkLayer(0);
MkLayer::MkLayer()
{
  Clear();
}

MkLayer::MkLayer(int)
{
  Clear();
}

void MkLayer::Clear()
{
  layername=0;
  depth=0;
  R_depth=0;
  Rt=0;
  Rsub=0;
  C=0;
  Ks=0;
  Pi=0;
  bearing=0;
}

void MkLayer::Import(MkGlobalVar &globalvar, int sec,MkSide side,int lay)
{
  layername=(LPCTSTR)globalvar.layername[lay+1];
  depth=((side==mkLeft)? globalvar.layer_depth_L[sec+1][lay+1]:globalvar.layer_depth_R[sec+1][lay+1]);
  R_depth=((side==mkLeft)? R_layer_depth_L[sec+1][lay+1]:R_layer_depth_R[sec+1][lay+1]);
  Rt=globalvar.Rt[sec+1][lay+1];
  Rsub=globalvar.Rsub[sec+1][lay+1];
  C=globalvar.C[sec+1][lay+1];
  Ks=globalvar.Ks[sec+1][lay+1];
  Pi=globalvar.Pi[sec+1][lay+1];
  bearing=::bearing[lay+1];
}

void MkLayer::Export(MkGlobalVar &globalvar, int sec,MkSide side,int lay)
{
    globalvar.layername[lay+1]=layername.c_str();
    ((side==mkLeft)? globalvar.layer_depth_L[sec+1][lay+1]:globalvar.layer_depth_R[sec+1][lay+1])=depth;
    ((side==mkLeft)? R_layer_depth_L[sec+1][lay+1]:R_layer_depth_R[sec+1][lay+1])=R_depth;
    globalvar.Rt[sec+1][lay+1]=Rt;
    globalvar.Rsub[sec+1][lay+1]=Rsub;
    globalvar.C[sec+1][lay+1]=C;
    globalvar.Ks[sec+1][lay+1]=Ks;
    globalvar.Pi[sec+1][lay+1]=Pi;
    ::bearing[lay+1]=bearing;
}
//---------------------------------------------------------------------------
MkLayers::MkLayers(int size,MkLayer *layers)
{
    if (size < 0) {
      MkDebug("::MkLayers - MkLayers(int size)");;
      return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0) {
       FLayer = NULL;
       return;
    }

    FLayer = new MkLayer[FSize];
    for (int i=0;i<FSize;i++) (*this)[i] = layers[i];
}

MkLayers::MkLayers(int size)
{
    if (size < 0) {
      MkDebug("::MkLayers - MkLayers(int size)");;
      return;
    }

    FSize = FSizeOfArray = size;

    if (FSizeOfArray == 0) {
       FLayer = NULL;
       return;
    }

    FLayer = new MkLayer[FSizeOfArray];
}

MkLayers::~MkLayers()
{
   FSizeOfArray = FSize = 0;
   if (FLayer) {
      delete[] FLayer;
      FLayer = NULL;
   }
}

void MkLayers::Initialize(int size)
{
    if (size < 0) {
      MkDebug("::MkLayers - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;

    FSize = FSizeOfArray = size;

    if (FSizeOfArray == 0) {
       if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
       FLayer = NULL;
       return;
    }

    if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
    FLayer = new MkLayer[FSizeOfArray];
}

void MkLayers::Initialize(int size,MkLayer *layers)
{

    if (size < 0 || layers == NULL) {
      MkDebug("::MkLayers - Initialize(int size)");;
      return;
    }
    if (FSizeOfArray == size) return;
    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0) {
       if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
       FLayer = NULL;
       return;
    }

    if (FLayer!=NULL) delete[] (MkLayer*)FLayer;
    FLayer = new MkLayer[FSizeOfArray];
    for (int i=0;i<FSizeOfArray;i++) FLayer[i] = layers[i];
}

int MkLayers::Grow(int delta)
{
    int i;
    MkLayer *layer=NULL;

    if (!(layer = new MkLayer[FSizeOfArray+delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        layer[i] = FLayer[i];
    for (i=FSize; i<FSizeOfArray+delta;i++)
        layer[i] = NullMkLayer;
    if (FLayer) {
       delete[] (MkLayer*)FLayer;
       FLayer = NULL;
    }
    FLayer = layer;
    FSizeOfArray = FSizeOfArray+delta;

    return FSizeOfArray;
}

int MkLayers::Shrink(int delta)
{
    int i;
    MkLayer *layer=NULL;

    if (!(layer = new MkLayer[FSizeOfArray-delta])) return FSizeOfArray;

    for (i = 0; i < FSize;i++)
        layer[i] = FLayer[i];
    for (i=FSize; i<FSizeOfArray-delta;i++)
        layer[i] = NullMkLayer;
    if (FLayer) {
       delete[] (MkLayer*)FLayer;
       FLayer = NULL;
    }
    FLayer = layer;
    FSizeOfArray = FSizeOfArray-delta;

    return FSizeOfArray;
}

bool MkLayers::Add(MkLayer &layer)
{
    int tmp=FSizeOfArray;
    bool flag=false;
    for (int i=0;i<FSize;i++) if (FLayer[i]==layer) flag=true;

    if(flag) return false;
    if (FSize>=FSizeOfArray) {
       Grow(FSize-FSizeOfArray+10);
       if (tmp==FSizeOfArray) return false;
    }
    FSize++;
    FLayer[FSize-1] = layer;

    return true;
}

bool MkLayers::Add(int index, MkLayer &layer)
{
    int tmp=FSizeOfArray;

    if(FSize>=FSizeOfArray) Grow(FSize-FSizeOfArray+1);
    if(tmp==FSizeOfArray) return false;

    for (int i=FSize-1;i>=index;i--)
      FLayer[i+1] = FLayer[i];
    FSize++;
    FLayer[index] = layer;
    return true;
}

bool MkLayers::Delete(MkLayer &layer)
{
    int i;
    for (i=0;i<FSize;i++) {
      if(FLayer[i] == layer) break;
    }
    if(i==FSize) return false;
    if(FLayer[i] == layer) {
      for (int j=i;j<FSize-1;j++)
        FLayer[j] = FLayer[j+1];
    }
    FSize--;
    FLayer[FSize] = NullMkLayer;
    return true;
}

bool MkLayers::Delete(int index)
{
    for (int j=index;j<FSize-1;j++)
        FLayer[j] = FLayer[j+1];

    FSize--;
    FLayer[FSize] = NullMkLayer;
    return true;
}

bool MkLayers::Clear()
{
   FSizeOfArray = FSize = 0;
   if (FLayer) {
      delete[] FLayer;
      FLayer = NULL;
   }
   return true;
}

MkLayer & MkLayers::operator[](int i)
{
    if (FSizeOfArray == 0) return NullMkLayer;
    if (i >= FSize && i < FSizeOfArray) FSize = i+1;

    if (i >=0 && i < FSize) return FLayer[i];
    else return NullMkLayer;
}

MkLayers & MkLayers::operator=(MkLayers &layers)
{
    int i;

    Clear();
    FSize = layers.FSize;
    FSizeOfArray = layers.FSizeOfArray;
    if (FSize == 0) {
       FLayer = NULL;
       return *this;
    }
    this->FLayer = new MkLayer[FSizeOfArray];

    for (i=0;i<FSize;i++)
      FLayer[i] = layers.FLayer[i];
    for (i=FSize;i<FSizeOfArray;i++)
      FLayer[i] = NullMkLayer;

    return *this;
}

bool MkLayers::operator==(MkLayers &layers)
{
  int i;

  if (FSize != layers.FSize) return false;
  for (i=0;i<FSize;i++)
    if (this->FLayer[i] != layers.FLayer[i]) return false;

  return true;
}

void MkLayers::Import(MkGlobalVar &globalvar, int sec,MkSide side)
{
  int i;
  Initialize(globalvar.layer_ea);
  for (i=0;i<FSize && i<10;i++)
    FLayer[i].Import(globalvar,sec,side,i);
}

void MkLayers::Export(MkGlobalVar &globalvar, int sec,MkSide side)
{
  int i;
  for (i=0;i<FSize && i<10;i++)
    FLayer[i].Export(globalvar,sec,side,i);
}
*/
//---------------------------------------------------------------------------
