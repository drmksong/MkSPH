//---------------------------------------------------------------------------
#pragma hdrstop
#include "MkTest.h"

#define GI 50
#define GJ 50
#define TI 10
#define TJ 10
//---------------------------------------------------------------------------
GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};  /* Gray diffuse light. */
GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};  /* Infinite light location. */
GLfloat n[6][3] = {  /* Normals for the 6 faces of a cube. */
  {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
  {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0} };
GLint faces[6][4] = {  /* Vertex indices for the 6 faces of a cube. */
  {0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
  {4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3} };
GLfloat v[8][3];  /* Will be filled in with X,Y,Z vertexes. */
MkParticles Particles(3500);
MkSmoothFunc SmoothFunc, SmoothFuncTemp;
MkSPHGrid SPHGrids[GI*GJ];
MkSPHGrid TGrids[TI*TJ];


float angleX=0.5,angleY = 0.5;
float deltaAngle = 0.5;
int specialKey;
float red, green, blue;

int particles()
{
  int i,j,k;
  float R, w,dw,d2w;
  int I,J;
  float X,Y,Z;

  SmoothFunc.SetupSmoothFunc(smtGaussian, 2, 0.04, 1);
  SmoothFuncTemp.SetupSmoothFunc(smtGaussian, 2, 0.1, 1);

  srand(time(NULL));
  /*  for (i=0;i<100;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = -1;
    Particles[i].YY = Particles[i].Y = (i%10)*0.2222222-1.0;//rand()%10000/10000.0-0.5;
    Particles[i].ZZ = Particles[i].Z = ((i/10)%10)*0.2222222-1.0;//rand()%10000/10000.0-0.5;
    Particles[i].SetTemp(0);
    Particles[i].SetTemp1(0);
    Particles[i].SetupTempColor();
  }

  for (;i<200;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = 1;
    Particles[i].YY = Particles[i].Y = (i%10)*0.2222222-1.0;//rand()%10000/10000.0-0.5;
    Particles[i].ZZ = Particles[i].Z = ((i/10)%10)*0.2222222-1.0;//rand()%10000/10000.0-0.5;
    Particles[i].SetTemp(100);
    Particles[i].SetTemp1(100);
    Particles[i].SetupTempColor();
  }
  */
  /*
  for (i=0;i<40;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = -1;
    Particles[i].YY = Particles[i].Y = i*0.05-1;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(0);
    Particles[i].SetTemp1(0);
    Particles[i].SetupTempColor();
  }

  for (i=40;i<80;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = 1;
    Particles[i].YY = Particles[i].Y = (i-40)*0.05-1;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(100);
    Particles[i].SetTemp1(100);
    Particles[i].SetupTempColor();
  }

  for (i=80;i<120;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = (i-80)*0.05-1;
    Particles[i].YY = Particles[i].Y = -1;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp((i-80)*2.5);
    Particles[i].SetTemp1((i-80)*2.5);
    Particles[i].SetupTempColor();
  }

  for (i=120;i<160;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = (i-120)*0.05-1;
    Particles[i].YY = Particles[i].Y = 1;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp((i-120)*2.5);
    Particles[i].SetTemp1((i-120)*2.5);
    Particles[i].SetupTempColor();
  }
  */
  for (i=0;i<Particles.GetSize();i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X =0.9* (rand()%10000/5000.0-1.0);// (i%10)*0.2222222-1.0;
    Particles[i].YY = Particles[i].Y =0.9* (rand()%10000/5000.0-1.0);// ((i/10)%10)*0.2222222-1.0;
    Particles[i].ZZ = Particles[i].Z = 0;// rand()%10000/5000.0-1.0;//((i/100)%10)*0.2222222-1.0;
    Particles[i].SetTemp(0);
    Particles[i].SetTemp1(0);
    Particles[i].SetRadius(fabs(Normal(0.018,0.003)));//0.01+0.04*(rand()%10000/10000.0));
    /*    if (Particles[i].X < -0.99){
      Particles[i].XX = Particles[i].X = -1.0;
      Particles[i].SetTemp(0);
      Particles[i].SetTemp1(0);
    }
    else if (Particles[i].X > 0.99){
      Particles[i].XX = Particles[i].X = 1.0;
      Particles[i].SetTemp(100);
      Particles[i].SetTemp1(100);
      } */
    Particles[i].SetupTempColor();
  }

  for (i=0;i<GI*GJ;i++) {
    I = i%GI;
    J = i/GJ;
    SPHGrids[i].Setup(-1+2.0/GI*I,-1+2.0/GJ*J,0,-1+2.0/GI*(I+1),-1+2.0/GJ*(J+1),0);
    SPHGrids[i].Setup(15);
  }
  for (i=0;i<Particles.GetSize();i++){
    X = Particles[i].X;
    Y = Particles[i].Y;
    I = int((X+1)*GI/2); if (I==GI) I = GI-1;
    J = int((Y+1)*GJ/2); if (J==GJ) J = GJ-1;
    k = I+J*GI;
    SPHGrids[k].Register(i);
    Particles[i].SetGridNum(k);
  }
  for (i=0;i<GI*GJ;i++) {
    printf("%d : ",i);
    SPHGrids[i].Out();
  }
  /*
  for (i=0;i<Particles.GetSize();i++) {
    printf("%3i-particle, (%10.3f, %10.3f, %10.3f) \n",i,Particles[i].X,Particles[i].Y,Particles[i].Z); 
  }

  for (j=1;j<Particles.GetSize();j++) {
    printf("dist = %f, Wij,ab(dist) = %f\n",CalDist(Particles[0],Particles[j]), Wijab(Particles[0],Particles[j]));
  }*/
  //getch();
}

float Wijab(MkParticle &mp1, MkParticle &mp2)
{
  float wij,dist,R,dx,dy,dz,drdx,drdy,drdz;
  dist = CalDist(mp1, mp2);
  R = dist/SmoothFunc.GetSmoothLen();
  dx = mp2.X - mp1.X;
  dy = mp2.Y - mp1.Y;
  dz = mp2.Z - mp1.Z;
  drdx = SmoothFunc.dRdX(dx,R);
  drdy = SmoothFunc.dRdX(dy,R);
  drdz = SmoothFunc.dRdX(dz,R);
  wij = drdx*(SmoothFunc.d2WdR2(R)-1/R*SmoothFunc.dWdR(R))*drdx
    +drdy*(SmoothFunc.d2WdR2(R)-1/R*SmoothFunc.dWdR(R))*drdy
    +drdz*(SmoothFunc.d2WdR2(R)-1/R*SmoothFunc.dWdR(R))*drdz
    +1/(R*SmoothFunc.GetSmoothLen()*SmoothFunc.GetSmoothLen())*SmoothFunc.dWdR(R);
  //  printf("R is %f, drdx is %f, drdy is %f, drdz is %f, wij is %f\n",R, drdx, drdy, drdz, wij);
  return wij;
}

void pos_bak(void) 
{
  static int i=0,iter = 0;
  static float mu=10,eps=0.01,constant=0.01,hi=SmoothFunc.GetSmoothLen(),hij=(hi+hi)/2.0;
  int j;
  float unval;
  float rho_j=1, m_j=1;
  static bool flag=true;
  MkVector Xi(3), rij(3), rhs(3);

for (i=0;i<Particles.GetSize();i++) {
  for (j=0;j<Xi.GetSize();j++) Xi[j] = 0;
  for (j=0;j<Particles.GetSize();j++){
    if (i==j) continue;
    float dist = CalDist(Particles[i], Particles[j]);
    if (dist > 2*hi) continue;
    rij.SetVector(Particles[j].X - Particles[i].X,Particles[j].Y - Particles[i].Y,Particles[j].Z - Particles[i].Z );

    rhs = rij*(mu*hi*(pow(hi/(dist+eps),2)+constant));
    Xi += rhs;
  }
  //  printf("X is %f, Y is %f, Z is %f, dx is %f, dy is %f, dz is %f\n",Particles[i].X, Particles[i].Y, Particles[i].Z, Xi[0],Xi[1], Xi[2]);
  Particles[i].X -= Xi[0];
  Particles[i].Y -= Xi[1];
  Particles[i].Z -= Xi[2];

  
  if (Particles[i].X<-0.99 || Particles[i].X>0.99 || Particles[i].Y<-0.99 || Particles[i].Y>0.99) {
    //Particles[i].X = rand()%10000/5000.0-1.00;
    //Particles[i].Y = rand()%10000/5000.0-1.00;
    //i--;
  }
}
  
  i++;
  //  if (i>=Particles.GetSize()) {
  //    i=0;
    if (mu>2.5) mu = mu/1.04;
    else if (mu<2.5) mu = mu/1.01;
    //    if (mu<0.5 && flag) {mu = 1.0;flag =false;} 
    printf("mu =  %f\n",mu);
    //  }
}

bool pos(void)
{
  static int i=0,iter = 0;
  static float mu=40,eps=0.01,constant=0.01,hi,hij;
  int j,k,l;
  float X,Y;
  float dx=0,dy=0;
  float unval,dist;
  float rho_j=1, m_j=1;
  static bool flag=true;
  MkVector Xi(3), rij(3), rhs(3);
  int gridnum;
  MkInt grids(9);
  int I,J,II,JJ;
  /*
  if ((18 < mu && mu < 18.2) || (15 < mu && mu < 15.2) || (12 < mu && mu < 12.2) || (9 < mu && mu < 9.1)|| (6 < mu && mu < 6.1) || (4 < mu && mu < 4.1) || (2 < mu && mu < 2.1)) {
    for (i=0;i<Particles.GetSize();i++) {
      for (j=0;j<Particles.GetSize();j++) {
	if (i==j) continue;
	dist = CalDist(Particles[i],Particles[j]);
	if (dist < hi*1.0) {
	  rij.SetVector(Particles[i].X - Particles[j].X,Particles[i].Y - Particles[j].Y,Particles[i].Z - Particles[j].Z );
	  rij.Normalize();
	  rhs = rij*hi;

	  MkPoint pnt;
	  pnt.Set(rhs[0],rhs[1],rhs[2]);
	  pnt.Rotate(0,0,90);
	  Particles[j].X += pnt.X;
	  Particles[j].Y += pnt.Y;
	  Particles[j].Z += pnt.Z;
	  
	  X = Particles[j].X;
	  Y = Particles[j].Y;
	  I = int((X+1)*15/2); if (I>=15) I = 14; else if(I<0) I = 0;
	  J = int((Y+1)*15/2); if (J>=15) J = 14; else if(J<0) J = 0;
	  k = I+J*15;

	  if (k!=Particles[j].GetGridNum()){
	    SPHGrids[Particles[j].GetGridNum()].Unregister(j);
	    SPHGrids[k].Register(j);
	    Particles[j].SetGridNum(k);
	  }
	  //printf("i = %d, j = %d\n",i,j);
	}
      }
    }
  }
  */
for (i=0;i<Particles.GetSize();i++){  
  for (j=0;j<3;j++) Xi[j] = 0;
  for (j=0;j<grids.getSzX();j++) grids[j]=-1;
  
  gridnum = Particles[i].GetGridNum();
  I = gridnum%GI;
  J = gridnum/GJ;

  II = I-1; JJ = J-1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[0] = k;

  II = I; JJ = J-1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[1] = k;
  
  II = I+1; JJ = J-1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[2] = k;

  II = I-1; JJ = J;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[3] = k;

  II = I; JJ = J;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[4] = k;
  
  II = I+1; JJ = J;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[5] = k;

  II = I-1; JJ = J+1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[6] = k;

  II = I; JJ = J+1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[7] = k;
  
  II = I+1; JJ = J+1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[8] = k;

  //  printf("%d %d %d %d %d %d %d %d %d\n", grids[0], grids[1], grids[2], grids[3], grids[4], grids[5], grids[6], grids[7], grids[8]);
 
  for (l=0;l<grids.getSzX();l++) {
    if (grids[l] == -1) continue;
    //    printf("l = %d, grids[l] = %d, SPHGrids[grids[l]].GetNumOfParticle = %d\n",l,grids[l],SPHGrids[grids[l]].GetNumOfParticle());
    for (k=0;k<SPHGrids[grids[l]].GetNumOfParticle();k++) {
      j = SPHGrids[grids[l]][k];
      if (j==-1) continue;
      if (i==j) continue;
      //      printf("i = %d, j = %d\n",i,j);
      dist = CalDist(Particles[i], Particles[j]);
      hij = (Particles[i].GetRadius()+Particles[j].GetRadius())/2.0;
      if (dist > 2*hij) continue;

      rij.SetVector(Particles[j].X - Particles[i].X,Particles[j].Y - Particles[i].Y,Particles[j].Z - Particles[i].Z );
      rhs = rij*(mu*hij*(pow(hij/(dist+eps),2)+constant));
      Xi += rhs;
    }
  }

  Particles[i].X -= Xi[0];
  Particles[i].Y -= Xi[1];
  Particles[i].Z -= Xi[2];

  dx += fabs(Xi[0]);
  dy += fabs(Xi[1]);
  if (mu > 5.5 && Particles[i].X<-0.99 || Particles[i].X>0.99 || Particles[i].Y<-0.99 || Particles[i].Y>0.99) {
    //Particles[i].X = rand()%10000/5000.0-1.00;
    //Particles[i].Y = rand()%10000/5000.0-1.00;
  }

  X = Particles[i].X;
  Y = Particles[i].Y;
  I = int((X+1)*GI/2); if (I>=GI) I = GI-1; else if(I<0) I = 0;
  J = int((Y+1)*GJ/2); if (J>=GJ) J = GJ-1; else if(J<0) J = 0;
  k = I+J*GI;

  if (k!=Particles[i].GetGridNum()){
    SPHGrids[Particles[i].GetGridNum()].Unregister(i);
    SPHGrids[k].Register(i);
    Particles[i].SetGridNum(k);
  }
}
  i++;
  //  if (i>=Particles.GetSize()) {
  //    i=0;
    if (mu>2.5) mu = mu/1.002;
    else if (mu<2.5) mu = mu/1.001;
    //    if (mu<0.5 && flag) {mu = 1.0;flag =false;} 
    printf("mu =  %f, dx = %f, dy = %f\n",mu,dx,dy);
    //    getch();

    //  }
    return (dx+dy < 0.001) ? true:false;
}

float repos(float Mu)
{
  static int i=0,iter = 0;
  static float mu=20,eps=0.01,constant=0.01,hi,hij;
  int j,k,l;
  float X,Y;
  float dx=0,dy=0;
  float unval,dist;
  float rho_j=1, m_j=1;
  static bool flag=true;
  MkVector Xi(3), rij(3), rhs(3);
  int gridnum;
  MkInt grids(9);
  int I,J,II,JJ;
  /*
  if ((18 < mu && mu < 18.2) || (15 < mu && mu < 15.2) || (12 < mu && mu < 12.2) || (9 < mu && mu < 9.1)|| (6 < mu && mu < 6.1) || (4 < mu && mu < 4.1) || (2 < mu && mu < 2.1)) {
    for (i=0;i<Particles.GetSize();i++) {
      for (j=0;j<Particles.GetSize();j++) {
	if (i==j) continue;
	dist = CalDist(Particles[i],Particles[j]);
	if (dist < hi*1.0) {
	  rij.SetVector(Particles[i].X - Particles[j].X,Particles[i].Y - Particles[j].Y,Particles[i].Z - Particles[j].Z );
	  rij.Normalize();
	  rhs = rij*hi;

	  MkPoint pnt;
	  pnt.Set(rhs[0],rhs[1],rhs[2]);
	  pnt.Rotate(0,0,90);
	  Particles[j].X += pnt.X;
	  Particles[j].Y += pnt.Y;
	  Particles[j].Z += pnt.Z;
	  
	  X = Particles[j].X;
	  Y = Particles[j].Y;
	  I = int((X+1)*15/2); if (I>=15) I = 14; else if(I<0) I = 0;
	  J = int((Y+1)*15/2); if (J>=15) J = 14; else if(J<0) J = 0;
	  k = I+J*15;

	  if (k!=Particles[j].GetGridNum()){
	    SPHGrids[Particles[j].GetGridNum()].Unregister(j);
	    SPHGrids[k].Register(j);
	    Particles[j].SetGridNum(k);
	  }
	  //printf("i = %d, j = %d\n",i,j);
	}
      }
    }
  }
  */
  mu = Mu;
for (i=0;i<Particles.GetSize();i++){  
  if(Particles[i].GetFixity()) continue;
  for (j=0;j<3;j++) Xi[j] = 0;
  for (j=0;j<grids.getSzX();j++) grids[j]=-1;
  
  gridnum = Particles[i].GetGridNum();
  I = gridnum%GI;
  J = gridnum/GJ;

  II = I-1; JJ = J-1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[0] = k;

  II = I; JJ = J-1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[1] = k;
  
  II = I+1; JJ = J-1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[2] = k;

  II = I-1; JJ = J;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[3] = k;

  II = I; JJ = J;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[4] = k;
  
  II = I+1; JJ = J;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[5] = k;

  II = I-1; JJ = J+1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[6] = k;

  II = I; JJ = J+1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[7] = k;
  
  II = I+1; JJ = J+1;
  k = II+JJ*GI;
  if(0<=k && k<GI*GJ) grids[8] = k;

  //  printf("%d %d %d %d %d %d %d %d %d\n", grids[0], grids[1], grids[2], grids[3], grids[4], grids[5], grids[6], grids[7], grids[8]);
 
  for (l=0;l<grids.getSzX();l++) {
    if (grids[l] == -1) continue;
    //    printf("l = %d, grids[l] = %d, SPHGrids[grids[l]].GetNumOfParticle = %d\n",l,grids[l],SPHGrids[grids[l]].GetNumOfParticle());
    for (k=0;k<SPHGrids[grids[l]].GetNumOfParticle();k++) {
      j = SPHGrids[grids[l]][k];
      if (j==-1) continue;
      if (i==j) continue;
      //      printf("i = %d, j = %d\n",i,j);
      dist = CalDist(Particles[i], Particles[j]);
      hij = (Particles[i].GetRadius()+Particles[j].GetRadius())/2.0;
      if (dist > 2*hij) continue;

      rij.SetVector(Particles[j].X - Particles[i].X,Particles[j].Y - Particles[i].Y,Particles[j].Z - Particles[i].Z );
      rhs = rij*(mu*hij*(pow(hij/(dist+eps),2)+constant));
      Xi += rhs;
    }
  }

  Particles[i].X -= Xi[0];
  Particles[i].Y -= Xi[1];
  Particles[i].Z -= Xi[2];

  dx += fabs(Xi[0]);
  dy += fabs(Xi[1]);
  if (mu > 5.5 && Particles[i].X<-0.99 || Particles[i].X>0.99 || Particles[i].Y<-0.99 || Particles[i].Y>0.99) {
    //Particles[i].X = rand()%10000/5000.0-1.00;
    //Particles[i].Y = rand()%10000/5000.0-1.00;
  }

  X = Particles[i].X;
  Y = Particles[i].Y;
  I = int((X+1)*GI/2); if (I>=GI) I = GI-1; else if(I<0) I = 0;
  J = int((Y+1)*GJ/2); if (J>=GJ) J = GJ-1; else if(J<0) J = 0;
  k = I+J*GI;

  if (k!=Particles[i].GetGridNum()){
    SPHGrids[Particles[i].GetGridNum()].Unregister(i);
    SPHGrids[k].Register(i);
    Particles[i].SetGridNum(k);
  }
}
  i++;
  //  if (i>=Particles.GetSize()) {
  //    i=0;
    if (mu>2.5) mu = mu/1.002;
    else if (mu<2.5) mu = mu/1.001;
    //    if (mu<0.5 && flag) {mu = 1.0;flag =false;} 
    printf("mu =  %f, dx = %f, dy = %f\n",mu,dx,dy);
    //    getch();

    //  }
    return (dx+dy);
}

bool boxing(void)
{
  static float y_min=10, x_min=10, x_max = -10, _min;
  static float is_first = true;
  static bool is_settled=true;
  static bool is_boxed=false;
  static float Mu=20;
  float ac_dist;
  int i;
  if (is_first) {
    for (i=0;i<Particles.GetSize();i++) {
      if (Particles[i].X < x_min) x_min = Particles[i].X;
      if (Particles[i].X > x_max) x_max = Particles[i].X;
      if (Particles[i].Y < y_min) y_min = Particles[i].Y;
    }
    is_first = false;
    _min = (fabs(x_min) > fabs(x_max))? fabs(x_min) : fabs(x_max);
    _min = (fabs(_min) > fabs(y_min)) ?fabs(_min) : fabs(y_min);  
  }

  is_boxed = _min <= 1.001 ? true : false;
  
  if (is_boxed) return is_boxed;

  if (Mu>2.5) Mu = Mu/1.002;
  else if (Mu<2.5) Mu = Mu/1.001;

  ac_dist = repos(Mu); 
  if(ac_dist < 0.01) is_settled = true;

  printf("         y_min = %f    ac_dist = %f        ",_min, ac_dist);

  if (is_settled) {
    _min -= 0.01;
    for (i=0;i<Particles.GetSize();i++) {
      if (Particles[i].Y < -_min) {
	Particles[i].Y = -_min;
	Particles[i].SetFixity(true);
      }
      if (Particles[i].X > _min) {
	Particles[i].X = _min;
	Particles[i].SetFixity(true);
      }
      if (Particles[i].X < -_min) {
	Particles[i].X = -_min;
	Particles[i].SetFixity(true);
      }

      Mu = 20;
      is_settled = false;
    }
  }
  return is_boxed;
}

void temp_init(void)
{
  int i,j;
  float t=0;
  for (i=0;i<Particles.GetSize();i++) {
    if (Particles[i].X > 0.99) {
      Particles[i].SetTemp(100);
      Particles[i].SetTemp1(100);
    }
    else if(Particles[i].X < -0.99) {
      Particles[i].SetTemp(0);
      Particles[i].SetTemp1(0);
    }
    Particles[i].SetupTempColor();
  }
}


void temp_init2(void)
{
  int i,j,k;
  int I,J;
  float X,Y;
  float t=0;
  for (i=0;i<Particles.GetSize();i++) {
    if (Particles[i].X > -0.03 && Particles[i].X < 0.03 && Particles[i].Y > -0.03 && Particles[i].Y < 0.03) {
      Particles[i].SetTemp(100);
      Particles[i].SetTemp1(100);
    }
    Particles[i].SetupTempColor();
  }

  for (i=0;i<TI*TJ;i++) {
    I = i%TI;
    J = i/TJ;
    TGrids[i].Setup(-1+2.0/TI*I,-1+2.0/TJ*J,0,-1+2.0/TI*(I+1),-1+2.0/TJ*(J+1),0);
    TGrids[i].Setup(15);
  }
  for (i=0;i<Particles.GetSize();i++){
    X = Particles[i].X;
    Y = Particles[i].Y;
    I = int((X+1)*TI/2); if (I==TI) I = TI-1; if (I<0) I = 0;
    J = int((Y+1)*TJ/2); if (J==TJ) J = TJ-1; if (J<0) J = 0;
    k = I+J*TI;
    TGrids[k].Register(i);
    Particles[i].SetGridNum(k);
  }

}

bool temp(void)
{
  static int i=0;
  static float t=0;
  static bool is_ini=false;

  static int iter = 0;
  int j,l;
  float X,Y;
  float unval,dist;
  static bool flag=true;
  int gridnum;
  MkInt grids(9);
  int I,J,II,JJ;

  float rho_i=1, K = 1, rho_j=1, k=1, C=1, m_j=1, sum=0,temp,dt=0.00000003,Q=0;

  if (!is_ini) {temp_init2(); is_ini = true;}
  for (i=0;i<Particles.GetSize();i++) {

    if ((Particles[i].X > -0.03 && Particles[i].X < 0.03 && Particles[i].Y > -0.03 && Particles[i].Y < 0.03)) continue;
    //
    //    if ((Particles[i].X > -0.2 && Particles[i].X < 0.2 && Particles[i].Y > -0.2 && Particles[i].Y < 0.2) ||(Particles[i].X > 0.99 ) ) continue;
    //if (Particles[i].X > 0.99 || Particles[i].X < -0.99) continue;
    //if (Particles[i].X > 0.99) continue; // || Particles[i].X < -0.99) continue;

    for (j=0;j<grids.getSzX();j++) grids[j]=-1;
  
    gridnum = Particles[i].GetGridNum();
    I = gridnum%TI;
    J = gridnum/TJ;

    II = I-1; JJ = J-1;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[0] = k;
    
    II = I; JJ = J-1;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[1] = k;
    
    II = I+1; JJ = J-1;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[2] = k;
    
    II = I-1; JJ = J;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[3] = k;
    
    II = I; JJ = J;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[4] = k;
    
    II = I+1; JJ = J;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[5] = k;
    
    II = I-1; JJ = J+1;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[6] = k;
    
    II = I; JJ = J+1;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[7] = k;
    
    II = I+1; JJ = J+1;
    k = II+JJ*TI;
    if(0<=k && k<TI*TJ) grids[8] = k;

    temp = Particles[i].GetTemp();
    sum = 0;
    
    for (l=0;l<grids.getSzX();l++) {
      if (grids[l] == -1) continue;
      for (k=0;k<TGrids[grids[l]].GetNumOfParticle();k++) {
	j = TGrids[grids[l]][k];
	if (j==-1) continue;
	//for (j=0;j<Particles.GetSize();j++) {
	if (i==j) continue;
		
	dist = CalDist(Particles[i], Particles[j]);
	sum += K*(Particles[i].GetTemp()-Particles[j].GetTemp())/dist*SmoothFuncTemp.dWdR(dist)*m_j/rho_j/rho_i/C;
		
	temp += dt*(sum+Q/rho_i/C);
	Particles[i].SetTemp1(temp);
	Particles[i].SetupTempColor();
	//  printf("%d-th particles is updated to temperature %f\n",i,temp);
	
	//  getch();
      }
    }
  }
  unval = 0;
  for (j=0;j<Particles.GetSize();j++) {
    unval += Particles[j].RefreshTemp();
  }
  t += dt;
  printf("time = %f, unval =  %f\n",t,unval); 
  
  return unval < 0.01 ? true:false;

}

void temp_back(void)
{
  static int i=0;
  static float t=0;
  static float unval=0;
  int j;
  float rho_i=1, rho_j=1, k=1, C=1, m_j=1, sum=0,temp,dt=0.001,Q=0;
  static bool is_ini=false;

  if (!is_ini) {temp_init(); is_ini = true;}
  //  printf("temp %d\n", i);
  if (Particles[i].X > 0.99 || Particles[i].X < -0.99) {
    i++;  
    if (i>=Particles.GetSize()){ 
      i=0;
      for (j=0;j<Particles.GetSize();j++) {
	Particles[j].RefreshTemp();
      }
    }
    return;
  }

  temp = Particles[i].GetTemp();
  sum = 0;
  for (j=0;j<Particles.GetSize();j++){
    if (i==j) continue;
    float dist = CalDist(Particles[i], Particles[j]);
    //    sum += k*Particles[j].GetTemp()*Wijab(Particles[i], Particles[j])*m_j/rho_j/rho_i/C;
    sum += k*(Particles[i].GetTemp()-Particles[j].GetTemp())/dist*SmoothFunc.dWdR(dist)*m_j/rho_j/rho_i/C;
  }
  
  temp += dt*(sum+Q/rho_i/C);
  Particles[i].SetTemp1(temp);
  Particles[i].SetupTempColor();
  printf("%d-th particles is updated to temperature %f\n",i,temp);

  //  getch();

  i++;
  if (i>=Particles.GetSize()) {
    i=0;
    unval = 0;
    for (j=0;j<Particles.GetSize();j++) {
      unval += Particles[j].RefreshTemp();
    }
    t += dt;
    printf("time = %f, unval =  %f\n",t,unval);
  }
  
}

void solve(void)
{
  static int i,j;
  float sum;
  float dist;
  int cnt;
  for (i=0;i<Particles.GetSize();i++) {
    if (Particles[i].X > 0.45 || Particles[i].X < -0.45) continue;
    sum = 0;cnt = 0;
    for (j=0;j<Particles.GetSize();j++) {
      if (Particles[j].X < 0.45 && Particles[j].X > -0.45) {
	dist = CalDist(Particles[i],Particles[j]);
	if (dist < 0.1) {
	  cnt++;
	  sum += Particles[j].GetTemp();
	}
      }
    }
    if (cnt>0) Particles[i].SetTemp(sum/cnt);
    Particles[i].SetupTempColor();
  }
}

void scatter(void)
{
  static int i=0,j,k;
  int _near[3];
  float dist[3];
  float sum;
  float d;
  int cnt;

  for (k=0;k<2;k++){
    dist[k] = 10;
    _near[k] = 0;
  }  
  if (i>Particles.GetSize()) i=0;
    if (Particles[i].X < -0.45 || Particles[i].X > 0.45) {i++;printf("%d-th particle skipped\n",i);return;};
  //  if (Particles[i].Y < -0.45 || Particles[i].Y > 0.45) {i++;printf("%d-th particle skipped\n",i);return;};
  //  if (Particles[i].Z < -0.45 || Particles[i].Z > 0.45) {i++;printf("%d-th particle skipped\n",i);return;};

  sum = 0;cnt = 0;
  for (j=0;j<Particles.GetSize();j++) {
    d = CalDist(Particles[i],Particles[j]);
    if (d < fabs(dist[0])){
      dist[1] = dist[0];
      _near[1] = _near[0];
      dist[0] = d*(Particles[i].X-Particles[j].X)/fabs(Particles[i].X-Particles[j].X);
      _near[0] = j;
    }    
  }
  float dx;
  dx = fabs(dist[0]*dist[1])*(1/dist[0]+1/dist[1])/2;
  Particles[i].X += dx;
  i++;
  printf("%d-th particle calculated\n",i);
  Particles[i].SetR(1.0);
  Particles[i].SetG(1.0);
  Particles[i].SetB(1.0);
  
}

void
drawBox(void)
{
  int i;

  for (i = 0; i < 6; i++) {
    glBegin(GL_QUADS);
    glNormal3fv(&n[i][0]);
    glVertex3fv(&v[faces[i][0]][0]);
    glVertex3fv(&v[faces[i][1]][0]);
    glVertex3fv(&v[faces[i][2]][0]);
    glVertex3fv(&v[faces[i][3]][0]);
    glEnd();
  }
}

void
display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //  drawBox();
  glPushMatrix();
  //glRotatef(angleY,1.0,0.0,0.0);
  //glRotatef(angleY,0.0,1.0,0.0);
  //glRotatef(angleY,0.0,0.0,1.0);
  glRotatef(180, 0.,0.,1.);

  Particles.Draw();
  glPopMatrix();
  // this is the new line
  // previously it was: angle++;
  angleY+=deltaAngle;

  glutSwapBuffers();
  //  printf("angleX = %f, angleY = %f\n",angleX, angleY);
}

void
init(void)
{
  /* Setup cube vertex data. */
  v[0][0] = v[1][0] = v[2][0] = v[3][0] = -1;
  v[4][0] = v[5][0] = v[6][0] = v[7][0] = 1;
  v[0][1] = v[1][1] = v[4][1] = v[5][1] = -1;
  v[2][1] = v[3][1] = v[6][1] = v[7][1] = 1;
  v[0][2] = v[3][2] = v[4][2] = v[7][2] = 1;
  v[1][2] = v[2][2] = v[5][2] = v[6][2] = -1;

  /* Enable a single OpenGL light. */
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);		       // Enable color


  /* Use depth buffering for hidden surface elimination. */
  glEnable(GL_DEPTH_TEST);

  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective( /* field of view in degree */ 40.0,
    /* aspect ratio */ 1.0,
    /* Z near */ 1.0, /* Z far */ 10.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0.0, 0.0, -2.5,  /* eye is at (0,0,5) */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    0.0, 1.0, 0.);      /* up is in positive Y direction */

  /* Adjust cube position to be asthetic angle. */
  glTranslatef(0.0, 0.0, 1.0);
  //  glRotatef(60, 1.0, 0.0, 0.0);
  //  glRotatef(-20, 0.0, 0.0, 1.0);
}

void processMouse(int button, int state, int x, int y) {


	specialKey = glutGetModifiers();
	// if both a mouse button, and the ALT key, are pressed  then
	if ((state == GLUT_DOWN) && 
			(specialKey == GLUT_ACTIVE_ALT)) {

		// set the color to pure red for the left button
		if (button == GLUT_LEFT_BUTTON) {
			red = 1.0; green = 0.0; blue = 0.0;
		}
		// set the color to pure green for the middle button
		else if (button == GLUT_MIDDLE_BUTTON) {
			red = 0.0; green = 1.0; blue = 0.0;
		}
		// set the color to pure blue for the right button
		else {
			red = 0.0; green = 0.0; blue = 1.0;
		}
	}
	//printf("mouse button pressed\n");
	temp();
}

void processMouseActiveMotion(int x, int y) {

  float width = glutGet(GLUT_WINDOW_X);
  float height = glutGet(GLUT_WINDOW_Y);

	// the ALT key was used in the previous function
	if (specialKey != GLUT_ACTIVE_ALT) {
		// setting red to be relative to the mouse 
		// position inside the window
		if (x < 0)
			red = 0.0;
		else if (x > width)
			red = 1.0;
		else
			red = ((float) x)/height;
		// setting green to be relative to the mouse 
		// position inside the window
		if (y < 0)
			green = 0.0;
		else if (y > width)
			green = 1.0;
		else
			green = ((float) y)/height;
		// removing the blue component.
		blue = 0.0;
	}
}

void processMouseEntry(int state) {
	if (state == GLUT_LEFT)
		deltaAngle = 0.0;
	else
		deltaAngle = 0.1;
	//	printf("mouse entry\n");
}

void processMousePassiveMotion(int x, int y) {

  float width = glutGet(GLUT_WINDOW_X);
  float height = glutGet(GLUT_WINDOW_Y);

	// User must press the SHIFT key to change the 
	// rotation in the X axis
	if (specialKey != GLUT_ACTIVE_SHIFT) {

		// setting the angle to be relative to the mouse 
		// position inside the window
		if (x < 0)
			angleX = 0.0;
		else if (x > width)
		  			angleX = 0.0;
		else
		  			angleX = 0.0;
	}
	//	printf("mouse passive motion\n");
}

void update(void) {
  static bool is_scattered = false;
  static bool is_converged = true;
  static bool is_boxed = false;
  //    angleX += 0.2f;
    if (angleX > 360) {
        angleX -= 360;
    }

    if (!is_scattered) {
      is_scattered = pos();
    }
    if (is_scattered && !is_boxed) is_boxed = boxing();
    //if (is_scattered && is_boxed && !is_converged) is_converged = temp();
    if (is_scattered && is_boxed && is_converged) exit(-1);
    
    glutPostRedisplay(); //Tell GLUT that the scene has changed
    
    //Tell GLUT to call update again in 25 milliseconds
    //    printf("update called\n");
    //    glutTimerFunc(5, update, 0);
}

int main(int argc, char **argv)
{
  particles();
   
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_ALPHA);
  glutCreateWindow("red 3D lighted cube");
  glutDisplayFunc(display);

  glutMouseFunc(processMouse);
  glutMotionFunc(processMouseActiveMotion);
  glutPassiveMotionFunc(processMousePassiveMotion);
  glutEntryFunc(processMouseEntry);

  init();
    
  //  glutTimerFunc(5, update, 0);
  glutIdleFunc(update);

  glutMainLoop();

  return 0;
}
#pragma package(smart_init)
