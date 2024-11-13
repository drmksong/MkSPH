//---------------------------------------------------------------------------
#pragma hdrstop
#include "MkTestSPH.h"

#define GI 20
#define GJ 20
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
MkParticles Particles(1100);
MkLines Lines(5000);
MkSmoothFunc SmoothFunc, SmoothFuncTemp;
MkSPHGrid SPHGrids[GI*GJ];
MkSPH mySPH;

float angleX=0.5,angleY = 0.5;
float deltaAngle = 0.5;
int specialKey, lineCount;
float red, green, blue;
/*
int particles()
{
  int i,j,k,n_bp;
  float R, w,dw,d2w;
  int I,J;
  float X,Y,Z;

  SmoothFunc.SetupSmoothFunc(smtGaussian, 2, 0.04, 1);
  SmoothFuncTemp.SetupSmoothFunc(smtGaussian, 2, 0.1, 1);

  srand(time(NULL));

  for (i=0;i<Particles.GetSize();i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X =0.55* (rand()%10000/5000.0-1)+0.2;// (i%10)*0.2222222-1.0;
    Particles[i].YY = Particles[i].Y =0.55* (rand()%10000/5000.0-1)+0.2;// ((i/10)%10)*0.2222222-1.0;
    Particles[i].ZZ = Particles[i].Z = 0;// rand()%10000/5000.0-1.0;//((i/100)%10)*0.2222222-1.0;
    Particles[i].SetTemp(0);
    Particles[i].SetTemp1(0);
    Particles[i].SetRadius(fabs(Normal(0.034,0.000)));//0.01+0.04*(rand()%10000/10000.0));
    Particles[i].SetupTempColor();
  }

  for (i=0;i<4;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = (i == 0 || i == 1) ? -1 : 1;
    Particles[i].YY = Particles[i].Y = (i == 0 || i == 3) ? -1 : 1;
    Particles[i].ZZ = Particles[i].Z = 0;// rand()%10000/5000.0-1.0;//((i/100)%10)*0.2222222-1.0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.034,0.000)));//0.01+0.04*(rand()%10000/10000.0));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
  }

  n_bp = int(2.0/0.034/2+0.5);

  for (i=4;i<4+n_bp-2;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = 2.0*(float(i-4+1)/float(n_bp-1))-1.0;
    Particles[i].YY = Particles[i].Y = -1;
    Particles[i].ZZ = Particles[i].Z = 0;// rand()%10000/5000.0-1.0;//((i/100)%10)*0.2222222-1.0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.034,0.000)));//0.01+0.04*(rand()%10000/10000.0));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
  }

  for (i=4+n_bp-2;i<4+(n_bp-2)*2;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = 2.0*(float(i-4-(n_bp-2)+1)/float(n_bp-1))-1.0;
    Particles[i].YY = Particles[i].Y = 1;
    Particles[i].ZZ = Particles[i].Z = 0;// rand()%10000/5000.0-1.0;//((i/100)%10)*0.2222222-1.0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.034,0.000)));//0.01+0.04*(rand()%10000/10000.0));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
  }

  for (i=4+(n_bp-2)*2;i<4+(n_bp-2)*3;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = -1;
    Particles[i].YY = Particles[i].Y = 2.0*(float(i-4-(n_bp-2)*2+1)/float(n_bp-1))-1.0;
    Particles[i].ZZ = Particles[i].Z = 0;// rand()%10000/5000.0-1.0;//((i/100)%10)*0.2222222-1.0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.034,0.000)));//0.01+0.04*(rand()%10000/10000.0));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
  }

  for (i=4+(n_bp-2)*3;i<4+(n_bp-2)*4;i++) {
    MkColor col=(MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = 1;
    Particles[i].YY = Particles[i].Y = 2.0*(float(i-4-(n_bp-2)*3+1)/float(n_bp-1))-1.0;
    Particles[i].ZZ = Particles[i].Z = 0;// rand()%10000/5000.0-1.0;//((i/100)%10)*0.2222222-1.0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.034,0.000)));//0.01+0.04*(rand()%10000/10000.0));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
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
}

bool pos(void)
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
    if (mu>2.5) mu = mu/1.003;
    else if (mu<2.5) mu = mu/1.002;
    printf("mu =  %f, dx = %f, dy = %f\n",mu,dx,dy);

    return (dx+dy < 0.01) ? true:false;
}
*/
float shake()
{
  for (int i=0;i<Particles.GetSize();i++){  
    if(Particles[i].GetFixity()) continue;
    Particles[i].X += rand()%10000/500000.0-0.01;
    Particles[i].Y += rand()%10000/500000.0-0.01;
  }
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
  MkVector<float> Xi(3), rij(3), rhs(3);
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

bool mesh(void)
{
  int i,j,lineCount = 0;;
  float dist, hij;
  for(i=0;i<Particles.GetSize();i++) {
    for(j=i+1;j<Particles.GetSize();j++) {
      dist = CalDist(Particles[i],Particles[j]);
      hij = (Particles[i].GetRadius()+Particles[j].GetRadius())/2.0;
      if (dist < 2*hij*1.4) {
	Lines[lineCount].SetLine(Particles[i].X, Particles[i].Y, Particles[j].X, Particles[j].Y); 
	lineCount++;
      }
    }
  }
  return true;
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

void drawBox(void)
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

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //  drawBox();
  glPushMatrix();
  //glRotatef(angleY,1.0,0.0,0.0);
  //glRotatef(angleY,0.0,1.0,0.0);
  //glRotatef(angleY,0.0,0.0,1.0);
  glRotatef(180, 0.,0.,1.);

  mySPH.GetParticles().Draw();
  // Lines.Draw();  disabled temorarily due to error of MkContainer does not have Draw member function
  glPopMatrix();
  // this is the new line
  // previously it was: angle++;
  angleY+=deltaAngle;

  glutSwapBuffers();
  //  printf("angleX = %f, angleY = %f\n",angleX, angleY);
}

void init(void)
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
	shake();
	//printf("mouse button pressed\n");
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
/*
void update(void) {
  static bool is_scattered = false;
  static bool is_meshed = false;

  if (!is_scattered) {
    is_scattered = pos();
  }
  if (is_scattered) is_meshed = mesh();
  //  if (is_scattered && is_meshed) exit(-1);
  
  glutPostRedisplay(); //Tell GLUT that the scene has changed
  
  //Tell GLUT to call update again in 25 milliseconds
  //    printf("update called\n");
  //    glutTimerFunc(5, update, 0);
}
*/
void update2(void) 
{
  mySPH.Run();
  glutPostRedisplay(); //Tell GLUT that the scene has changed
}

void InitSPH(void)
{
  mySPH.Initialize(2000,20,20);
  mySPH.SetBound(-1,-1,1,1);
  mySPH.SetRefDensity(0.5);
  mySPH.SetGammaPolytropic(1.0);
  mySPH.SetSoundSpeed(340);
  mySPH.SetGravity(0,-10.0);
  mySPH.SetDeltaTime(0.000001);

  mySPH.GenParticle();
  mySPH.InitPosition();
  mySPH.InitSPHGrid();

}

int main(int argc, char **argv)
{
  //  particles();
  InitSPH();
   
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
  //  glutIdleFunc(update);
  glutIdleFunc(update2);

  glutMainLoop();

  return 0;
}
#pragma package(smart_init)
