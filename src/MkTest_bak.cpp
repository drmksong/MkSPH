//---------------------------------------------------------------------------
#pragma hdrstop
#include "MkTest.h"
//---------------------------------------------------------------------------

void nurbtest()
{
  printf("%f",sin(1.0));

  MkPoint rp; 
  MkNurbs mn;
  mn.MakeItSphere();

  for (float u=0;u<1.0;u+=0.05) {
    for (float v=0;v<1.0;v+=0.05) {
      rp = mn(u,v);
      printf("Point at (%f,%f) is (%f, %f, %f)\n",u,v,rp.X,rp.Y,rp.Z);
    }
  }
}

void exceptiontest()
{
  try {
    MkFloat a(0);
  }
  catch (MkFloat::Size s)
  {
    printf("Exception catched : Size of MkFloat is %d,%d,%d\n",s.X,s.Y,s.Z);
  }
  try {MkInt b(0,1,5);}
  catch (MkInt::Size s)
  {
    printf("Exception catched : Size is MkInt %d,%d,%d\n",s.X,s.Y,s.Z);
  }

  try {
    MkMatrix m(1,0);
  }
  catch (MkFloat::Size s)
  {
    printf("Exception catched : Size of Matrix is %d,%d,%d\n",s.X,s.Y,s.Z);
  }

  printf("End of test\n");
}

void polytest()
{
  int i,j,n=17;
  MkPolygon pa,pb,pc,sa,sb;
  MkPolygons ps;
  MkPoints pnts;
  pa.Initialize(n);
  pb.Initialize(n);
  pc.Initialize(n);
  pa.SetCloseness(true);
  pb.SetCloseness(true);
  pc.SetCloseness(true);

  for(i=0;i<n;i++) {
    pa[i].X = cos(2*i*3.14/n);
    pa[i].Y = sin(2*i*3.14/n);
    pb[i].X = cos(2*i*3.14/n)+1.9;
    pb[i].Y = sin(2*i*3.14/n);
    pc[i].X = cos(2*i*3.14/n)+0.95;
    pc[i].Y = sin(2*i*3.14/n)+1.5;
  }
  pa.FindInter(pb,pnts);

  puts(" pa");
  for(i=0;i<pa.GetSize();i++) {
    printf("%f %f \n",pa[i].X, pa[i].Y);
  }

  puts("\n pb");
  for(i=0;i<pb.GetSize();i++) {
    printf("%f %f \n",pb[i].X, pb[i].Y);
  }

  puts("\n pc");
  for(i=0;i<pc.GetSize();i++) {
    printf("%f %f \n",pc[i].X, pc[i].Y);
  }

  puts("\n pnts");

  for(i=0;i<pnts.GetSize();i++) {
    printf("%f %f \n",pnts[i].X, pnts[i].Y);
  }
  puts("\n ps");
  /*
  if(pnts.GetSize()==2) {
    pa.Extract(pnts,sa);
    pb.Extract(pnts,sb);
    sb.InverseDirection();
    sa.Merge(sb);
  }
  */
  pa.BoolUni(pb,ps);
  //  if(ps.GetSize()!=1) exit(-1);

  for(i=0;i<ps.GetSize();i++) {
    for(j=0;j<ps[i].GetSize();j++) {
      printf("%f %f \n",ps[i][j].X, ps[i][j].Y);
    }
    puts("\n");
  }

  puts("\n ps2");

  pa = ps[0];
  pa.BoolSub(pc,ps);

  for(i=0;i<ps.GetSize();i++) {
    for(j=0;j<ps[i].GetSize();j++) {
      printf("%f %f \n",ps[i][j].X, ps[i][j].Y);
    }
    puts("\n");
  }
  /*
  puts("\n npoly");
  pa = ps[0];
  int npoly = pa.FindPoly(pc,btInt);
  printf("number of polygon %d  \n",npoly);
  */
}

void setlayer(MkLayers &layer)
{
  MkRect rect;
  MkPoint ori;

  layer.Initialize(2);

  ori.SetPoint(-2.0,-11.0);
  rect.SetOrigin(ori);
  rect.SetWidth(5.0);
  rect.SetHeight(11.0);
  layer[0].SetRect(rect);
  layer[0].SetWetUnitWeight(0,3);
  layer[0].SetSubUnitWeight(0,2);
  layer[0].SetCohesion(1,0);
  layer[0].SetFriction(0,45.0);
  layer[0].SetHorSubReact(0,10000);
  layer[0].SetWetUnitWeight(1,3);
  layer[0].SetSubUnitWeight(1,2);
  layer[0].SetCohesion(1,0);
  layer[0].SetFriction(1,45.0);
  layer[0].SetHorSubReact(1,10000);
  layer[0].CalcCoeff();

  ori.SetPoint(-2.0,-13.0);
  rect.SetOrigin(ori);
  rect.SetWidth(5.0);
  rect.SetHeight(2.0);
  layer[1].SetRect(rect);
  layer[1].SetWetUnitWeight(0,2.3);
  layer[1].SetSubUnitWeight(0,1.3);
  layer[1].SetCohesion(0,5);
  layer[1].SetFriction(0,35.0);
  layer[1].SetHorSubReact(0,60000);
  layer[1].SetWetUnitWeight(1,2.3);
  layer[1].SetSubUnitWeight(1,1.3);
  layer[1].SetCohesion(1,5);
  layer[1].SetFriction(1,35.0);
  layer[1].SetHorSubReact(1,60000);
  layer[1].CalcCoeff();

  layer.SetupPress();
}

void setlayer4pap(MkLayers &layer)
{
  int i;
  MkRect rect;
  MkPoint ori;

  layer.Initialize(5);

  for(i=0;i<5;i++) {
    ori.SetPoint(-2.0,-i-1-((i==4)?6:0));
    rect.SetOrigin(ori);
    rect.SetWidth(5.0);
    rect.SetHeight(1.0+((i==4)?6:0));
    layer[i].SetRect(rect);
    layer[i].SetWetUnitWeight(0,2);
    layer[i].SetSubUnitWeight(0,1);
    layer[i].SetCohesion(1,0);
    layer[i].SetFriction(0,30.0);
    layer[i].SetHorSubReact(0,200*(i+1));
    layer[i].SetWetUnitWeight(1,2);
    layer[i].SetSubUnitWeight(1,1);
    layer[i].SetCohesion(1,0);
    layer[i].SetFriction(1,30.0);
    layer[i].SetHorSubReact(1,200*(i+1));
    layer[i].CalcCoeff();
  }
  layer.SetupPress();
}

void setpile(MkPiles &pile)
{
  MkLine l;

  pile.Initialize(1);

  l.SetLine(0,-10,0,0);
  pile[0].SetLine(l);
  pile[0].SetArea(0.01348);
  //  pile[0].SetYoungMod(2.1e7);
  //  pile[0].SetSecMomentY(0.000204);
  //  pile[0].SetSecMomentZ(0.000204);
  pile[0].SetYoungMod(10000);
  pile[0].SetSecMomentY(1);
  pile[0].SetSecMomentZ(1);
  pile[0].SetSpacing(2.0);
  pile[0].SetT1(0.1);
  pile[0].SetT2(0.1);
  pile[0].SetHeight(0.5);
  pile[0].SetWidth(0.5);
  pile[0].SetWeight(3);
  pile[0].SetShearTor(1);
  pile[0].SetYieldMom(1);
}

void setpile4pap(MkPiles &pile)
{
  MkLine l;

  pile.Initialize(1);

  l.SetLine(0,-10,0,0);
  pile[0].SetLine(l);
  pile[0].SetArea(0.01348);
  //  pile[0].SetYoungMod(2.1e7);
  //  pile[0].SetSecMomentY(0.000204);
  //  pile[0].SetSecMomentZ(0.000204);
  pile[0].SetYoungMod(10000);
  pile[0].SetSecMomentY(1);
  pile[0].SetSecMomentZ(1);
  pile[0].SetSpacing(2.0);
  pile[0].SetT1(0.1);
  pile[0].SetT2(0.1);
  pile[0].SetHeight(0.5);
  pile[0].SetWidth(0.5);
  pile[0].SetWeight(3);
  pile[0].SetShearTor(1);
  pile[0].SetYieldMom(1);
}

void setcut(MkCuts &cut,MkPiles &pile)
{
  cut.Initialize(4);

  cut[0].SetDepth(-1.0);
  cut[0].SetWall(&pile[0],0);

  cut[1].SetDepth(-2.0);
  cut[1].SetWall(&pile[0],0);
//  cut[1].SetWall(&pile[1],1);

  cut[2].SetDepth(-3.0);
  cut[2].SetWall(&pile[0],0);
//  cut[2].SetWall(&pile[1],1);

  cut[3].SetDepth(-4.0);
  cut[3].SetWall(&pile[0],0);
//  cut[3].SetWall(&pile[1],1);
}

void setbc(MkBndConds &bc)
{
  bc.Initialize(1);
  MkRangeTree rt;
  MkDOFs dof(6);

  setbcrange(rt);
  dof[0].SetType(doftXDis,bndtFix); // temp
  dof[1].SetType(doftYDis,bndtFix);
  dof[2].SetType(doftZDis,bndtFix);
  dof[3].SetType(doftXAng,bndtFix);
  dof[4].SetType(doftYAng,bndtFix);
  dof[5].SetType(doftZAng,bndtFix);

  bc[0].SetDOFs(dof);
  bc[0].SetRangeTree(rt);
}

void setbcrange(MkRangeTree &rt)
{
  MkPoint p;
  MkRect r;

  MkRangeShape rs;
  p.SetPoint(-0.01,-10.01);
  r.SetOrigin(p);
  r.SetWidth(1.02);
  r.SetHeight(0.02);
  rs.SetShape(&r);

  rt.SetRoot(&rs);
}
/*

void setload(MkLoads &load,MkWall &w1,MkWall &w2)
{
}

void setlloadrange(MkRangeTree &rt)
{
}

void setrloadrange(MkRangeTree &rt)
{
}

void setrankine(MkLoads &load,MkLayers &lay, MkCuts &cut, MkFills &fill, MkWall &w0)
{
}

void setsubreact(MkSubreacts &sub, MkLayers &lay, MkCuts &cut, MkFills &fill, MkPile &pile0)
{
}
*/

/*
simptest:: is for comprison of results, SUNEX, SAP and ESCOT.
ks : 100~10000
EI : 2500~160000
*/
void simpletest()
{
  int i,j;
  float space=0.25;
  FILE *fp;

  MkSection sec;
  MkAnalyticSection &ana = sec.GetAnalyticSection();
  ana.AddAnalysis(atBeamSpring);
  //  ana.AddAnalysis(atMori); //original

  MkPiles &pile=sec.GetPiles();
  MkExcavStep &es=sec.GetExcavStep();
  MkLayers &lay=sec.GetLayers();
  MkCuts &cut=sec.GetCuts();
  MkFills &fill=sec.GetFills();
  MkLoads &load=sec.GetAnalyticSection().GetLoads();
  MkSubreacts &sub=sec.GetAnalyticSection().GetSubreacts();
  MkBndConds &bc=sec.GetBndConds();

  setpile(pile);
  setcut(cut,pile);
  setlayer(lay);
  setbc(bc);
  //  setrankine(load,lay,cut,fill,pile[0]);
  //  setsubreact(sub,lay,cut,fill,pile[0]);

  es.Install(0,ettPile,0);
  es.Cut(0,ettCut,0,cut[0].GetDepth());
  es.Cut(1,ettCut,1,cut[1].GetDepth());
  es.Cut(2,ettCut,2,cut[2].GetDepth());

  fp = fopen("test0912.tot","w");
  fclose(fp);

  sec.SetSpacing(space);
  sec.SetupAllNodes();
  sec.GetAnalyticSection().GetAnalysis()[0]->SetFileName("test0912.tot");

  sec.Solve(2,atBeamSpring);
  //  sec.Solve(2,atMori); // original
  sec.Out("test0912_bs.out");
}

/*
forpaper:: is for the KSRM annual symposium.
case01: case nodal force, nodal spring
case02: case nodal force, dist spring
case03: case dist force, nodal spring
case04: case dist force, dist spring
contents : 4 different cases of loading and ground spring calculation method.
*/

void case01()
{
  int i,j;
  float space=0.5;
  FILE *fp;

  MkSection sec;
  MkAnalyticSection &ana = sec.GetAnalyticSection();
  ana.AddAnalysis(atMori); //original
  ana.GetAnalysis()[0]->SetLoadingType(ltNodal);
  ana.GetAnalysis()[0]->SetSpringType(stNodal);

  MkPiles &pile=sec.GetPiles();
  MkExcavStep &es=sec.GetExcavStep();
  MkLayers &lay=sec.GetLayers();
  MkCuts &cut=sec.GetCuts();
  MkFills &fill=sec.GetFills();
  MkLoads &load=sec.GetAnalyticSection().GetLoads();
  MkSubreacts &sub=sec.GetAnalyticSection().GetSubreacts();
  MkBndConds &bc=sec.GetBndConds();

  //  setpile4pap(pile);
  setpile(pile);
  setcut(cut,pile);
  setlayer4pap(lay);
  //setlayer(lay);
  setbc(bc);
  //  setrankine(load,lay,cut,fill,pile[0]);
  //  setsubreact(sub,lay,cut,fill,pile[0]);

  es.Install(0,ettPile,0);
  es.Cut(0,ettCut,0,cut[0].GetDepth());
  //es.Cut(1,ettCut,1,cut[1].GetDepth());
  //es.Cut(2,ettCut,2,cut[2].GetDepth());
  //es.Cut(3,ettCut,3,cut[3].GetDepth());

  fp = fopen("test1009_a.tot","w");
  fclose(fp);

  sec.SetSpacing(space);
  sec.SetupAllNodes();
  sec.GetAnalyticSection().GetAnalysis()[0]->SetFileName("test1009_a.tot");

  sec.Solve(2,atMori);
  sec.Out("test1009_a.out");
}

void case02()
{
  int i,j;
  float space=0.25;
  FILE *fp;

  MkSection sec;
  MkAnalyticSection &ana = sec.GetAnalyticSection();
  ana.AddAnalysis(atMori); //original
  ana.GetAnalysis()[0]->SetLoadingType(ltNodal);
  ana.GetAnalysis()[0]->SetSpringType(stDistributal);

  MkPiles &pile=sec.GetPiles();
  MkExcavStep &es=sec.GetExcavStep();
  MkLayers &lay=sec.GetLayers();
  MkCuts &cut=sec.GetCuts();
  MkFills &fill=sec.GetFills();
  MkLoads &load=sec.GetAnalyticSection().GetLoads();
  MkSubreacts &sub=sec.GetAnalyticSection().GetSubreacts();
  MkBndConds &bc=sec.GetBndConds();

  //  setpile4pap(pile);
  setpile(pile);
  setcut(cut,pile);
  setlayer4pap(lay);
  //setlayer(lay);
  setbc(bc);
  //  setrankine(load,lay,cut,fill,pile[0]);
  //  setsubreact(sub,lay,cut,fill,pile[0]);

  es.Install(0,ettPile,0);
  es.Cut(0,ettCut,0,cut[0].GetDepth());
  es.Cut(1,ettCut,1,cut[1].GetDepth());
  //es.Cut(2,ettCut,2,cut[2].GetDepth());
  //es.Cut(3,ettCut,3,cut[3].GetDepth());

  fp = fopen("test1009_b.tot","w");
  fclose(fp);

  sec.SetSpacing(space);
  sec.SetupAllNodes();
  sec.GetAnalyticSection().GetAnalysis()[0]->SetFileName("test1009_b.tot");

  sec.Solve(2,atMori);
  sec.Out("test1009_b.out");
}

void case03()
{
   int i,j;
  float space=0.5;
  FILE *fp;

  MkSection sec;
  MkAnalyticSection &ana = sec.GetAnalyticSection();
  ana.AddAnalysis(atMori); //original
  ana.GetAnalysis()[0]->SetLoadingType(ltDistributal);
  ana.GetAnalysis()[0]->SetSpringType(stNodal);

  MkPiles &pile=sec.GetPiles();
  MkExcavStep &es=sec.GetExcavStep();
  MkLayers &lay=sec.GetLayers();
  MkCuts &cut=sec.GetCuts();
  MkFills &fill=sec.GetFills();
  MkLoads &load=sec.GetAnalyticSection().GetLoads();
  MkSubreacts &sub=sec.GetAnalyticSection().GetSubreacts();
  MkBndConds &bc=sec.GetBndConds();

  //  setpile4pap(pile);
  setpile(pile);
  setcut(cut,pile);
  setlayer4pap(lay);
  //setlayer(lay);
  setbc(bc);
  //  setrankine(load,lay,cut,fill,pile[0]);
  //  setsubreact(sub,lay,cut,fill,pile[0]);

  es.Install(0,ettPile,0);
  es.Cut(0,ettCut,0,cut[0].GetDepth());
  //es.Cut(1,ettCut,1,cut[1].GetDepth());
  //es.Cut(2,ettCut,2,cut[2].GetDepth());
  //es.Cut(3,ettCut,3,cut[3].GetDepth());

  fp = fopen("test1009_c.tot","w");
  fclose(fp);

  sec.SetSpacing(space);
  sec.SetupAllNodes();
  sec.GetAnalyticSection().GetAnalysis()[0]->SetFileName("test1009_c.tot");

  sec.Solve(2,atMori);
  sec.Out("test1009_c.out");
}

void case04()
{
  int i,j;
  float space=0.5;
  FILE *fp;

  MkSection sec;
  MkAnalyticSection &ana = sec.GetAnalyticSection();
  ana.AddAnalysis(atMori); //original
  ana.GetAnalysis()[0]->SetLoadingType(ltDistributal);
  ana.GetAnalysis()[0]->SetSpringType(stDistributal);

  MkPiles &pile=sec.GetPiles();
  MkExcavStep &es=sec.GetExcavStep();
  MkLayers &lay=sec.GetLayers();
  MkCuts &cut=sec.GetCuts();
  MkFills &fill=sec.GetFills();
  MkLoads &load=sec.GetAnalyticSection().GetLoads();
  MkSubreacts &sub=sec.GetAnalyticSection().GetSubreacts();
  MkBndConds &bc=sec.GetBndConds();

  //  setpile4pap(pile);
  setpile(pile);
  setcut(cut,pile);
  setlayer4pap(lay);
  //setlayer(lay);
  setbc(bc);
  //  setrankine(load,lay,cut,fill,pile[0]);
  //  setsubreact(sub,lay,cut,fill,pile[0]);

  es.Install(0,ettPile,0);
  es.Cut(0,ettCut,0,cut[0].GetDepth());
  //es.Cut(1,ettCut,1,cut[1].GetDepth());
  //es.Cut(2,ettCut,2,cut[2].GetDepth());
  //es.Cut(3,ettCut,3,cut[3].GetDepth());

  fp = fopen("test1009_d.tot","w");
  fclose(fp);

  sec.SetSpacing(space);
  sec.SetupAllNodes();
  sec.GetAnalyticSection().GetAnalysis()[0]->SetFileName("test1009_d.tot");

  sec.Solve(2,atMori);
  sec.Out("test1009_d.out");
}

void forpaper()
{
  //case01(); // nodal force, nodal spring
  case02(); // nodal force, dist spring
  //case03(); //dist force, nodal spring
  //case04(); // dist force, dist spring
}

int main(void)
{
  //simpletest(); // comparison test SUNEX,ESCOT,SAP
  forpaper();
  printf("hello\n");
  return 0;
}
#pragma package(smart_init)
