//---------------------------------------------------------------------------
#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#include "MkSPH.hpp"

MkSPH::MkSPH()
{
  NumOfParticles = 0;
  GI = 0;
  GJ = 0;

  Particles.Clear();
  SmoothFunc.SetupSmoothFunc(smtGaussian, 2, 1, 1);
  SmoothFuncTemp.SetupSmoothFunc(smtGaussian, 2, 1, 1);

  SPHGrids.Clear();
}

MkSPH::~MkSPH()
{
  Particles.Clear();
  SPHGrids.Clear();
}

bool MkSPH::Initialize()
{
  NumOfParticles = 0;
  GI = 0;
  GJ = 0;

  Particles.Clear();

  SPHGrids.Clear();
  return true;
}

bool MkSPH::Initialize(int np, int gi, int gj)
{
  char s[256];
  NumOfParticles = np;
  GI = gi;
  GJ = gj;

  Particles.Initialize(NumOfParticles);
  //sprintf(s,"Particle num %d\n",NumOfParticles);
  //MkDebug(s);
  //sprintf(s,"Size of Particle  %d\n",Particles.GetSize());
  //MkDebug(s);

  //getch();

  SPHGrids.Initialize(GI * GJ);
}

bool MkSPH::InitAccel()
{
  int i;
  for (i = 0; i < Particles.GetSize(); i++)
  {
    MkParticle &par = Particles[i];
    par.GetXAccel(2) = 0;
    par.GetYAccel(2) = 0;
    par.GetXAccel(1) = 0;
    par.GetYAccel(1) = 0;
    par.GetXAccel(0) = 0;
    par.GetYAccel(0) = 0;
  }
  return true;
}

bool MkSPH::InitVelocity()
{
  int i;
  for (i = 0; i < Particles.GetSize(); i++)
  {
    MkParticle &par = Particles[i];
    par.GetXVelocity(2) = 0;
    par.GetYVelocity(2) = 0;
    par.GetXVelocity(1) = 0;
    par.GetYVelocity(1) = 0;
    par.GetXVelocity(0) = 0;
    par.GetYVelocity(0) = 0;
  }
  return true;
}

bool MkSPH::InitPosition()
{
  int i;
  for (i = 0; i < Particles.GetSize(); i++)
  {
    MkParticle &par = Particles[i];
    par.GetXPosition(2) = par.X;
    par.GetYPosition(2) = par.Y;
    par.GetXPosition(1) = par.X;
    par.GetYPosition(1) = par.Y;
    par.GetXPosition(0) = par.X;
    par.GetYPosition(0) = par.Y;
  }
  return true;
}

bool MkSPH::GenParticle2()
{
  int i, j, n_bp_x, n_bp_y;
  float R, w, dw, d2w;
  int I, J;
  float X, Y, Z;
  float rad;
  char s[256];

  sprintf(s, "GenParticles2 01 ::Size of Particle  %d\n", Particles.GetSize());
  MkDebug(s);

  for (i = 0; i < 5; i++)
  {
    for (j = 0; j < 5; j++)
    {
      I = j * 5 + i;
      MkColor col = (MkColor)0xFF00FF;
      Particles[I].XX = Particles[I].X = (UpX - LowX) * i / 4.0 + LowX;
      Particles[I].YY = Particles[I].Y = (UpY - LowY) * j / 4.0 + LowY;
      Particles[I].ZZ = Particles[I].Z = 0;
      Particles[I].SetTemp(I * 4);
      Particles[I].SetTemp1(I * 4);
      Particles[I].SetRadius(fabs(Normal(0.22, 0.000)));
      Particles[I].SetupTempColor();
      Particles[I].SetFixity(true);
      Particles[i].SetViscosity(0.01);
      rad = Particles[I].GetRadius();
      Particles[I].SetMass(0.4);
      //sprintf(s,"particle gen no %d i %d j %d \n",Particles.GetSize(),i,j);
      //MkDebug(s);
      //getch();
    }
  }

  sprintf(s, "GenParticles2 02 ::Size of Particle  %d\n", Particles.GetSize());
  MkDebug(s);

  //Particles[12].XX = Particles[12].X = -0.25;
  //Particles[12].YY = Particles[12].Y -= 0.1;
  //Particles[12].YY = Particles[7].Y -= 0.1;
  //Particles[12].SetFixity(false);
  //Particles[7].SetFixity(false);

  for (i = 0; i < Particles.GetSize(); i++)
  {
    if ((Particles[i].X > -0.6) && (Particles[i].X < 0.6) && (Particles[i].Y > -0.2) && (Particles[i].Y < 0.2))
    {
      Particles[i].SetFixity(false);
      Particles[i].Y -= 0.15;
      Particles[i].YY = Particles[i].Y;
    }
    else
    {
      Particles[i].SetTemp(0);
      Particles[i].SetTemp1(0);
      Particles[i].SetupTempColor();
    }
    //sprintf(s,"bnd cond  i %d  \n",i);
    //MkDebug(s);
    //getch();
  }
  /*
  MkDebug("!!!");
  for (i=0;i<Particles.GetSize();i++) {
      MkColor col=(MkColor)0xFF00FF;
      Particles[i].XX = Particles[i].X = 0.7*(UpX-LowX)*(rand() % 1000)/1000.0+LowX/2;
      Particles[i].YY = Particles[i].Y = 0.7*(UpY-LowY)*(rand() % 1000)/1000.0+LowY/2;
      Particles[i].ZZ = Particles[i].Z = 0;
      Particles[i].SetTemp(50);
      Particles[i].SetTemp1(50);
      Particles[i].SetRadius(fabs(Normal(0.55,0.000)));
      Particles[i].SetupTempColor();
      Particles[i].SetFixity(false);
      Particles[i].SetViscosity(0.2);
      rad = Particles[i].GetRadius();
      Particles[i].SetMass(1);
      xi = Particles[i].GetXPosition(0);
  }

  for (i=0;i<7;i++) {
    Particles[i].SetFixity(true);
    Particles[i].XX = Particles[i].X = (UpX-LowX)*i/6.0+LowX;
    Particles[i].YY = Particles[i].Y = LowY;
      Particles[i].SetTemp(0);
      Particles[i].SetTemp1(0);
      Particles[i].SetupTempColor();
  }

  for (i=7;i<14;i++){
    Particles[i].SetFixity(true);
    Particles[i].XX = Particles[i].X = (UpX-LowX)*(i-7)/6.0+LowX;
    Particles[i].YY = Particles[i].Y = UpY;
      Particles[i].SetTemp(30);
      Particles[i].SetTemp1(30);
      Particles[i].SetupTempColor();
  }

  for (i=14;i<21;i++) {
    Particles[i].SetFixity(true);
    Particles[i].XX = Particles[i].X = LowX;
    Particles[i].YY = Particles[i].Y = (UpY-LowY)*(i-14)/6.0+LowY;
      Particles[i].SetTemp(60);
      Particles[i].SetTemp1(60);
      Particles[i].SetupTempColor();
  }

  for (i=21;i<28;i++) {
    Particles[i].SetFixity(true);
    Particles[i].XX = Particles[i].X = UpX;
    Particles[i].YY = Particles[i].Y = (UpY-LowY)*(i-21)/6.0+LowY;
      Particles[i].SetTemp(90);
      Particles[i].SetTemp1(90);
      Particles[i].SetupTempColor();
  }
  */
}

bool MkSPH::GenParticle()
{
  int i, j, k, n_bp_x, n_bp_y;
  float R, w, dw, d2w;
  int I, J;
  float X, Y, Z;
  float rad;

  char s[256];
  sprintf(s, "LX %f LY %f UX %f UY %f \n", LowX, LowY, UpX, UpY);
  MkDebug(s);

  for (i = 0; i < Particles.GetSize(); i++)
  {
    MkColor col = (MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = 0.90 * (UpX - LowX) * (rand() % 10000 / 10000.0 - 1.0) - LowX - (UpX - LowX) * 0.05;
    Particles[i].YY = Particles[i].Y = 0.90 * (UpY - LowY) * (rand() % 10000 / 10000.0 - 1.0) - LowY - (UpY - LowY) * 0.05;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(0);
    Particles[i].SetTemp1(0);
    Particles[i].SetRadius(fabs(Normal(0.55, 0.000)));
    Particles[i].SetupTempColor();
    rad = Particles[i].GetRadius();
    Particles[i].SetMass(3.141592 * rad * rad);
  }

  for (i = 0; i < 4; i++)
  {
    MkColor col = (MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = (i == 0 || i == 1) ? LowX : UpX;
    Particles[i].YY = Particles[i].Y = (i == 0 || i == 3) ? LowY : UpY;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.55, 0.000)));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
    rad = Particles[i].GetRadius();
    Particles[i].SetMass(3.141592 * rad * rad);
  }

  n_bp_x = int((UpX - LowX) / 0.05 / 2.0 + 0.5);
  n_bp_y = int((UpY - LowY) / 0.05 / 2.0 + 0.5);

  for (i = 4; i < 4 + n_bp_x - 2; i++)
  {
    MkColor col = (MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = (UpX - LowX) * (float(i - 4 + 1) / float(n_bp_x - 1)) + LowX;
    Particles[i].YY = Particles[i].Y = LowY;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.55, 0.000)));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
    rad = Particles[i].GetRadius();
    Particles[i].SetMass(3.141592 * rad * rad);
  }

  for (i = 4 + n_bp_x - 2; i < 4 + (n_bp_x - 2) * 2; i++)
  {
    MkColor col = (MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = (UpX - LowX) * (float(i - 4 - (n_bp_x - 2) + 1) / float(n_bp_x - 1)) + LowX;
    Particles[i].YY = Particles[i].Y = UpY;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.55, 0.000)));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
    rad = Particles[i].GetRadius();
    Particles[i].SetMass(3.141592 * rad * rad);
  }

  for (i = 4 + (n_bp_x - 2) * 2; i < 4 + (n_bp_x - 2) * 2 + (n_bp_y - 2); i++)
  {
    MkColor col = (MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = LowX;
    Particles[i].YY = Particles[i].Y = (UpY - LowY) * (float(i - 4 - (n_bp_x - 2) * 2 + 1) / float(n_bp_y - 1)) + LowY;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.55, 0.000)));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
    rad = Particles[i].GetRadius();
    Particles[i].SetMass(3.141592 * rad * rad);
  }

  for (i = 4 + (n_bp_x - 2) * 2 + (n_bp_y - 2); i < 4 + (n_bp_x - 2) * 2 + (n_bp_y - 2) * 2; i++)
  {
    MkColor col = (MkColor)0xFF00FF;
    Particles[i].XX = Particles[i].X = UpX;
    Particles[i].YY = Particles[i].Y = (UpY - LowY) * (float(i - 4 - (n_bp_x - 2) * 2 - (n_bp_y - 2) * 1 + 1) / float(n_bp_y - 1)) + LowY;
    Particles[i].ZZ = Particles[i].Z = 0;
    Particles[i].SetTemp(50);
    Particles[i].SetTemp1(50);
    Particles[i].SetRadius(fabs(Normal(0.55, 0.000)));
    Particles[i].SetupTempColor();
    Particles[i].SetFixity(true);
    rad = Particles[i].GetRadius();
    Particles[i].SetMass(3.141592 * rad * rad);
  }
}

bool MkSPH::InitSPHGrid()
{
  int i, k;
  int I, J;
  float X, Y;
  for (i = 0; i < GI * GJ; i++)
  {
    I = i % GI;
    J = i / GJ;
    SPHGrids[i].Setup(LowX + (UpX - LowX) / GI * I, LowY + (UpY - LowY) / GJ * J, 0, LowX + (UpX - LowX) / GI * (I + 1), LowY + (UpY - LowY) / GJ * (J + 1), 0);
    SPHGrids[i].Setup(15);
  }
  for (i = 0; i < Particles.GetSize(); i++)
  {
    X = Particles[i].X;
    Y = Particles[i].Y;
    I = int((X - LowX) * GI / (UpX - LowX));
    if (I == GI)
      I = GI - 1;
    J = int((Y - LowY) * GJ / (UpY - LowY));
    if (J == GJ)
      J = GJ - 1;
    k = I + J * GI;
    SPHGrids[k].Register(i);
    Particles[i].SetGridNum(k);
  }
  for (i = 0; i < GI * GJ; i++)
  {
    printf("%d : ", i);
    SPHGrids[i].Out();
  }
}

bool MkSPH::UpdateSPHGrid()
{
  static int i = 0, iter = 0;
  static float mu = 20, eps = 0.01, constant = 0.01, hi, hij;
  int j, k, l;
  float X, Y;
  float dx = 0, dy = 0;
  float unval, dist;
  float rho_j = 1, m_j = 1;
  static bool flag = true;
  MkVector<float> Xi(3), rij(3), rhs(3);
  int gridnum;
  MkInt grids(9);
  int I, J, II, JJ;

  for (i = 0; i < Particles.GetSize(); i++)
  {
    X = Particles[i].X;
    Y = Particles[i].Y;
    I = int((X - LowX) * GI / (UpX - LowX));
    if (I >= GI)
      I = GI - 1;
    else if (I < 0)
      I = 0;
    J = int((Y - LowY) * GJ / (UpY - LowY));
    if (J >= GJ)
      J = GJ - 1;
    else if (J < 0)
      J = 0;
    k = I + J * GI;

    if (k != Particles[i].GetGridNum())
    {
      SPHGrids[Particles[i].GetGridNum()].Unregister(i);
      SPHGrids[k].Register(i);
      Particles[i].SetGridNum(k);
    }
  }
}

bool MkSPH::ComputeDensity()
{
  static int i = 0, iter = 0;
  static float mu = 20, eps = 0.01, constant = 0.01, hi, hij;
  char s[256];

  int j, k, l, cnt = 0;
  float X, Y;
  float dx = 0, dy = 0;
  float unval, dist;
  float rho_j = 1, m_j = 1;
  static bool flag = true;
  MkVector<float> Xi(3), rij(3), rhs(3);
  int gridnum;
  MkInt grids(9);
  int I, J, II, JJ;
  float density;
  static float maxd = -1000, mind = 1000;

  sprintf(s, "Number of Particles %d \n", Particles.GetSize());
  MkDebug(s);

  for (i = 0; i < Particles.GetSize(); i++)
  {
    density = 0;
    cnt = 0;
    //    if(Particles[i].GetFixity()) continue;

    for (j = 0; j < 3; j++)
      Xi[j] = 0;
    for (j = 0; j < grids.getSzX(); j++)
      grids[j] = -1;

    gridnum = Particles[i].GetGridNum();
    I = gridnum % GI;
    J = gridnum / GJ;

    II = I - 1;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[0] = k;

    II = I;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[1] = k;

    II = I + 1;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[2] = k;

    II = I - 1;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[3] = k;

    II = I;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[4] = k;

    II = I + 1;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[5] = k;

    II = I - 1;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[6] = k;

    II = I;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[7] = k;

    II = I + 1;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[8] = k;

    cnt = 0;
    for (l = 0; l < grids.getSzX(); l++)
    {
      if (grids[l] == -1)
        continue;

      for (k = 0; k < SPHGrids[grids[l]].GetNumOfParticle(); k++)
      {
        j = SPHGrids[grids[l]][k];
        if (j == -1)
          continue;
        if (i == j)
          continue;

        dist = CalDist(Particles[i], Particles[j]);
        hij = (Particles[i].GetRadius() + Particles[j].GetRadius()) / 2.0;
        //if (dist > 2*hij) continue;

        density += Particles[j].GetMass() * SmoothFunc.W(dist);
        cnt++;
      }
    }
    Particles[i].SetDensity(density);
    sprintf(s, "grid 0: %d 1: %d 2: %d 3: %d 4: %d 5: %d 6: %d 7: %d 8: %d ", grids[0], grids[1], grids[2], grids[3], grids[4], grids[5], grids[6], grids[7], grids[8]);
    MkDebug(s);
    sprintf(s, "%d-th cnt %d mass %f density is %f\n", i, cnt, Particles[i].GetMass(), density);
    MkDebug(s);

    /*
    maxd = max(maxd,density);
    mind = min(mind,density);
    sprintf(s,"maxd %f mind %f\n",maxd, mind);
    //MkDebug(s);
    if (maxd > 1) {
      Particles[i].SetTemp1(100*density/maxd);
      Particles[i].SetupTempColor();
      sprintf(s,"density %f maxd %f density color %f\n",density,maxd,50*density/maxd);
      MkDebug(s);
    }
    */
  }
}

bool MkSPH::ComputeDensityRate()
{
  return true;
}

bool MkSPH::ComputePressure()
{
  int i;
  float pressure, den;
  char s[256];
  static float maxd = -1000, mind = 1000;
  for (i = 0; i < Particles.GetSize(); i++)
  {
    den = Particles[i].GetDensity();
    if (fabs(Gamma_P) < 0.001)
      return false;
    if (fabs(RefDensity) < 0.001)
      return false;

    pressure = SoundSpeed * SoundSpeed / Gamma_P * RefDensity * (pow(den / RefDensity, Gamma_P) - 1);
    Particles[i].SetPressure(pressure);
    char s[256];
    sprintf(s, "%d-th mass %f density is %f pressure is %f\n", i, Particles[i].GetMass(), Particles[i].GetDensity(), pressure);
    MkDebug(s);

    /*
    maxd = max(maxd,pressure);
    mind = min(mind,pressure);
    sprintf(s,"maxd %f mind %f\n",maxd, mind);
    //MkDebug(s);
    if (maxd > 1) {
      Particles[i].SetTemp1(100*pressure/maxd);
      Particles[i].SetupTempColor();
      sprintf(s,"pressure %f maxd %f pressure color %f\n",pressure,maxd,100*pressure/maxd);
      //      MkDebug(s);
    }
    */
  }

  return true;
}

bool MkSPH::ComputeInternalForce2()
{
  int i, j;
  float xf, yf, dist, angle, force;
  float mi, mj, di, dj, pi, pj, smfc;
  MkLine line;
  char s[256];

  for (i = 0; i < Particles.GetSize(); i++)
  {
    xf = yf = 0;
    for (j = 0; j < Particles.GetSize(); j++)
    {
      if (i == j)
        continue;

      mi = Particles[i].GetMass();
      mj = Particles[j].GetMass();
      di = Particles[i].GetDensity();
      dj = Particles[j].GetDensity();
      pi = Particles[i].GetPressure();
      pj = Particles[j].GetPressure();

      line.SetLine(Particles[i], Particles[j]);
      dist = line.GetLength();
      angle = line.GetTheta();

      smfc = SmoothFunc.dWdR(dist);

      if ((abs(di) < 0.001) || (abs(dj) < 0.001))
        force = 0;

      //if ((abs(di) > 0.001) && (abs(dj) > 0.001)) {
      force = -mi * mj * (pi / di / di + pj / dj / dj) * SmoothFunc.dWdR(dist);
      //force = -mi*mj/(di*dj)*(pi+pj)/2*SmoothFunc.dWdR(dist);
      sprintf(s, "i %d j %d dist %f angle %f mi * mj %f pi/di/di %f pj/dj/dj  %f smfc %f force %f \n ", i, j, dist, angle, mi * mj, pi / di / di, pj / dj / dj, smfc, force);
      //MkDebug(s);
      //}
      xf += force * cos(angle * 3.14159 / 180.0);
      yf += force * sin(angle * 3.14159 / 180.0);
    }
    xf = xf - XGravity;
    yf = yf - YGravity;
    Particles[i].SetXForceInt(xf);
    Particles[i].SetYForceInt(yf);

    sprintf(s, "%d-th mass %f density is %f pressure %f xf %f yf %f\n", i, Particles[i].GetMass(), Particles[i].GetDensity(), Particles[i].GetPressure(), xf, yf); //Particles[i].GetXForceInt(),Particles[i].GetYForceInt());
    MkDebug(s);
  }
}

bool MkSPH::ComputeInternalForce() // Mohaghan 1992
{

  static int i = 0, iter = 0;
  static float mu = 20, eps = 0.01, constant = 0.01, hi, hij;
  int j, k, l, cnt;
  float X, Y;
  float dx = 0, dy = 0;
  float unval;
  float rho_j = 1, m_j = 1;
  static bool flag = true;
  MkVector<float> Xi(3), rij(3), rhs(3);
  int gridnum;
  MkInt grids(9);
  int I, J, II, JJ;
  float density;

  static float xf = 0, yf = 0;
  float force, pressure, dist, angle;
  float mi, mj, di, dj, pi, pj;
  MkLine line;
  static float maxd = -1000000, mind = 1000000;
  char s[256];

  for (i = 0; i < Particles.GetSize(); i++)
  {
    xf = 0;
    yf = 0;
    //if(Particles[i].GetFixity()) continue;

    for (j = 0; j < 3; j++)
      Xi[j] = 0;
    for (j = 0; j < grids.getSzX(); j++)
      grids[j] = -1;

    gridnum = Particles[i].GetGridNum();
    I = gridnum % GI;
    J = gridnum / GJ;

    II = I - 1;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[0] = k;

    II = I;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[1] = k;

    II = I + 1;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[2] = k;

    II = I - 1;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[3] = k;

    II = I;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[4] = k;

    II = I + 1;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[5] = k;

    II = I - 1;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[6] = k;

    II = I;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[7] = k;

    II = I + 1;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ && 0 <= II && II < GI && 0 <= JJ && JJ < GJ)
      grids[8] = k;
    cnt = 0;
    for (l = 0; l < grids.getSzX(); l++)
    {
      if (grids[l] == -1)
        continue;

      for (k = 0; k < SPHGrids[grids[l]].GetNumOfParticle(); k++)
      {
        j = SPHGrids[grids[l]][k];
        if (j == -1)
          continue;
        if (i == j)
          continue;
        cnt++;
        mi = Particles[i].GetMass();
        mj = Particles[j].GetMass();
        di = Particles[i].GetDensity();
        dj = Particles[j].GetDensity();
        pi = Particles[i].GetPressure();
        pj = Particles[j].GetPressure();

        line.SetLine(Particles[i], Particles[j]);
        dist = line.GetLength();
        angle = line.GetTheta();

        static float smfc;
        smfc = SmoothFunc.dWdR(dist);

        if (abs(di) < 0.001 || abs(dj) < 0.001)
          force = 0;
        if (abs(di) > 0.001 && abs(dj) > 0.001)
        {
          force = -mi * mj * (pi / di / di + pj / dj / dj) * smfc;
          //kbhit();
          //getch();
          //	sprintf(s,"dist %f mi %f mj %f di %f dj %f pi %f pj  %f smfc %f force %f \n",dist, mi,mj,di,dj,pi,pj,smfc,force);
          sprintf(s, "mi * mj %f pi/di/di %f pj/dj/dj  %f smfc %f force %f ", mi * mj, pi / di / di, pj / dj / dj, smfc, force);
          MkDebug(s);

          if (abs(force) < 0.001)
          {
            MkDebug(s);
            sprintf(s, "forc %f \n", force);
            MkDebug(s);
            //kbhit;
          }
        }
        MkDebug(s);
        xf += force * cos(angle);
        yf += force * sin(angle);
      }
    }
    sprintf(s, "%d-th mass %f density is %f pressure %f xf %f yf %f\n", i, Particles[i].GetMass(), Particles[i].GetDensity(), Particles[i].GetPressure(), Particles[i].GetXForceInt(), Particles[i].GetYForceInt());
    MkDebug(s);

    xf = xf - XGravity;
    yf = yf - YGravity;
    Particles[i].SetXForceInt(xf);
    Particles[i].SetYForceInt(yf);

    maxd = max(maxd, xf);
    mind = min(mind, xf);
    sprintf(s, "xf %f maxd %f mind %f\n", xf, maxd, mind);
    //    MkDebug(s);

    Particles[i].SetTemp1(100 * (xf - mind) / (maxd - mind));
    Particles[i].SetupTempColor();
    sprintf(s, "xf %f maxd %f mind %f xf color %f\n", xf, maxd, mind, 100 * (xf - mind) / (maxd - mind));
    //MkDebug(s);
    printf("cnt = %d\n", cnt);
  }
  //  char ch;
  //if(kbhit()) cin.get(ch);
  return true;
}

bool MkSPH::ComputeViscousForce2()
{
  int i, j;
  float xf, yf, dist, angle, force;
  float mi, mj, di, dj, pi, pj, vi, vj, vxi, vyi, vxj, vyj, xi, yi, xj, yj, smfc;
  MkLine line;
  char s[256];

  for (i = 0; i < Particles.GetSize(); i++)
  {
    xf = yf = 0;
    for (j = 0; j < Particles.GetSize(); j++)
    {
      if (i == j)
        continue;

      mi = Particles[i].GetMass();
      mj = Particles[j].GetMass();
      di = Particles[i].GetDensity();
      dj = Particles[j].GetDensity();
      pi = Particles[i].GetPressure();
      pj = Particles[j].GetPressure();
      vi = Particles[i].GetViscosity();
      vj = Particles[j].GetViscosity();
      vxi = Particles[i].GetXVelocity(0);
      vxj = Particles[j].GetXVelocity(0);
      vyi = Particles[i].GetYVelocity(0);
      vyj = Particles[j].GetYVelocity(0);
      xi = Particles[i].GetXPosition(0);
      xj = Particles[j].GetXPosition(0);
      yi = Particles[i].GetYPosition(0);
      yj = Particles[j].GetYPosition(0);

      line.SetLine(Particles[i], Particles[j]);
      dist = line.GetLength();
      angle = line.GetTheta();

      smfc = SmoothFunc.dWdR(dist);

      if ((abs(di) < 0.001) || (abs(dj) < 0.001))
        force = 0;

      //xf += mi*mj/(di*dj)*(vi+vj)/2*(vxj-vxi)*SmoothFunc.d2WdR2(dist);//*cos(angle*3.14159/180.0);
      //yf += mi*mj/(di*dj)*(vi+vj)/2*(vyj-vyi)*SmoothFunc.d2WdR2(dist);//*sin(angle*3.14159/180.0);
      //xf -= mi*mj/(di*dj)*(vi+vj)/2*(vxj-vxi)*SmoothFunc.dWdR(dist);//*cos(angle*3.14159/180.0);
      //yf -= mi*mj/(di*dj)*(vi+vj)/2*(vyj-vyi)*SmoothFunc.dWdR(dist);//*sin(angle*3.14159/180.0);

      xf += mi * mj / (di * dj) * (vi + vj) * (vxj - vxi) * (xj - xi) / (dist * dist + 0.01 * SmoothLen * SmoothLen) * SmoothFunc.dWdR(dist);
      yf += mi * mj / (di * dj) * (vi + vj) * (vyj - vyi) * (yj - yi) / (dist * dist + 0.01 * SmoothLen * SmoothLen) * SmoothFunc.dWdR(dist);
    }
    //    xf = xf - XGravity;
    //    yf = yf - YGravity;
    Particles[i].SetXForceVis(xf);
    Particles[i].SetYForceVis(yf);

    //sprintf(s,"%d-th mass %f density is %f pressure %f viscous xf %f viscous yf %f\n",i,Particles[i].GetMass(), Particles[i].GetDensity(),Particles[i].GetPressure(),xf,yf);//Particles[i].GetXForceInt(),Particles[i].GetYForceInt());
    sprintf(s, "%d-th vyi %f vyj is %f viscous xf %f viscous yf %f\n", i, vyi, vyj, xf, yf); //Particles[i].GetXForceInt(),Particles[i].GetYForceInt());
    MkDebug(s);
  }
  //getch();
}

bool MkSPH::ComputeArtificialViscousForce2()
{
  int i, j;
  float xf, yf, dist, angle, force;
  float mi, mj, di, dj, pi, pj, vi, vj, vxi, vyi, vxj, vyj, xi, yi, xj, yj, smfc;
  float ci, cj, hi, hj, z1, z2, phi_ij_x, phi_ij_y;
  MkLine line;
  char s[256];

  ci = cj = SoundSpeed;
  hi = hj = SmoothLen;
  z1 = 1.0;
  z2 = 2 * z1;

  for (i = 0; i < Particles.GetSize(); i++)
  {
    xf = yf = 0;
    for (j = 0; j < Particles.GetSize(); j++)
    {
      if (i == j)
        continue;

      mi = Particles[i].GetMass();
      mj = Particles[j].GetMass();
      di = Particles[i].GetDensity();
      dj = Particles[j].GetDensity();
      pi = Particles[i].GetPressure();
      pj = Particles[j].GetPressure();
      vi = Particles[i].GetViscosity();
      vj = Particles[j].GetViscosity();
      vxi = Particles[i].GetXVelocity(0);
      vxj = Particles[j].GetXVelocity(0);
      vyi = Particles[i].GetYVelocity(0);
      vyj = Particles[j].GetYVelocity(0);
      xi = Particles[i].GetXPosition(0);
      xj = Particles[j].GetXPosition(0);
      yi = Particles[i].GetYPosition(0);
      yj = Particles[j].GetYPosition(0);

      line.SetLine(Particles[i], Particles[j]);
      dist = line.GetLength();
      angle = line.GetTheta();

      smfc = SmoothFunc.dWdR(dist);
      if ((vxj - vxi) * (xj - xi) < 0)
        phi_ij_x = (vxj - vxi) * (xj - xi) / (dist * dist + 0.01 * hi * hj);
      else
        phi_ij_x = 0;
      if ((vyj - vyi) * (yj - yi) < 0)
        phi_ij_y = (vyj - vyi) * (yj - yi) / (dist * dist + 0.01 * hi * hj);
      else
        phi_ij_y = 0;

      if ((abs(di) < 0.001) || (abs(dj) < 0.001))
        force = 0;

      //xf += mi*mj/(di*dj)*(vi+vj)/2*(vxj-vxi)*SmoothFunc.d2WdR2(dist);//*cos(angle*3.14159/180.0);
      //yf += mi*mj/(di*dj)*(vi+vj)/2*(vyj-vyi)*SmoothFunc.d2WdR2(dist);//*sin(angle*3.14159/180.0);
      //xf -= mi*mj/(di*dj)*(vi+vj)/2*(vxj-vxi)*SmoothFunc.dWdR(dist);//*cos(angle*3.14159/180.0);
      //yf -= mi*mj/(di*dj)*(vi+vj)/2*(vyj-vyi)*SmoothFunc.dWdR(dist);//*sin(angle*3.14159/180.0);

      xf -= mi * mj * (-z1 * (ci + cj) * (hi + hj) * phi_ij_x / 2 / (di + dj) + z2 * pow((hi + hj) * phi_ij_x, 2) / 2 / (di + dj)) * SmoothFunc.dWdR(dist);
      yf -= mi * mj * (-z1 * (ci + cj) * (hi + hj) * phi_ij_y / 2 / (di + dj) + z2 * pow((hi + hj) * phi_ij_y, 2) / 2 / (di + dj)) * SmoothFunc.dWdR(dist);
    }
    //    xf = xf - XGravity;
    //    yf = yf - YGravity;
    Particles[i].SetXForceArtVis(xf);
    Particles[i].SetYForceArtVis(yf);

    //sprintf(s,"%d-th mass %f density is %f pressure %f a viscous xf %f a viscous yf %f\n",i,Particles[i].GetMass(), Particles[i].GetDensity(),Particles[i].GetPressure(),xf,yf);//Particles[i].GetXForceInt(),Particles[i].GetYForceInt());
    sprintf(s, "%d-th vyi %f vyj is %f  a viscous xf %f a viscous yf %f\n", i, vyi, vyj, xf, yf); //Particles[i].GetXForceInt(),Particles[i].GetYForceInt());
    MkDebug(s);
  }
  //getch();
}

bool MkSPH::ComputeViscousForce()
{

  static int i = 0, iter = 0;
  static float mu = 20, eps = 0.01, constant = 0.01, hi, hij;
  int j, k, l;
  float X, Y;
  float dx = 0, dy = 0;
  float unval;
  float rho_j = 1, m_j = 1;
  static bool flag = true;
  MkVector<float> Xi(3), rij(3), rhs(3);
  int gridnum;
  MkInt grids(9);
  int I, J, II, JJ;
  float density;

  float xf = 0, yf = 0, xforce, yforce, pressure, dist, angle;
  float mi, mj, di, dj, pi, pj, vi, vj, vxi, vyi, vxj, vyj;
  MkLine line;

  for (i = 0; i < Particles.GetSize(); i++)
  {
    xf = 0;
    yf = 0;

    if (Particles[i].GetFixity())
      continue;

    for (j = 0; j < 3; j++)
      Xi[j] = 0;
    for (j = 0; j < grids.getSzX(); j++)
      grids[j] = -1;

    gridnum = Particles[i].GetGridNum();
    I = gridnum % GI;
    J = gridnum / GJ;

    II = I - 1;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[0] = k;

    II = I;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[1] = k;

    II = I + 1;
    JJ = J - 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[2] = k;

    II = I - 1;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[3] = k;

    II = I;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[4] = k;

    II = I + 1;
    JJ = J;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[5] = k;

    II = I - 1;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[6] = k;

    II = I;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[7] = k;

    II = I + 1;
    JJ = J + 1;
    k = II + JJ * GI;
    if (0 <= k && k < GI * GJ)
      grids[8] = k;

    for (l = 0; l < grids.getSzX(); l++)
    {
      if (grids[l] == -1)
        continue;

      for (k = 0; k < SPHGrids[grids[l]].GetNumOfParticle(); k++)
      {
        j = SPHGrids[grids[l]][k];
        if (j == -1)
          continue;
        if (i == j)
          continue;

        mi = Particles[i].GetMass();
        mj = Particles[j].GetMass();
        di = Particles[i].GetDensity();
        dj = Particles[j].GetDensity();
        pi = Particles[i].GetPressure();
        pj = Particles[j].GetPressure();
        vi = Particles[i].GetViscosity();
        vj = Particles[j].GetViscosity();
        vxi = Particles[i].GetXVelocity(0);
        vxj = Particles[j].GetXVelocity(0);
        vyi = Particles[i].GetYVelocity(0);
        vyj = Particles[j].GetYVelocity(0);

        line.SetLine(Particles[i], Particles[j]);
        dist = line.GetLength();
        angle = line.GetTheta();

        xforce += mi * mj / (di * dj) * (vi + vj) / 2 * (vxj - vxi) * SmoothFunc.d2WdR2(dist);
        yforce += mi * mj / (di * dj) * (vi + vj) / 2 * (vyj - vyi) * SmoothFunc.d2WdR2(dist);
      }
    }
    Particles[i].SetXForceVis(xforce);
    Particles[i].SetYForceVis(yforce);
  }
}

bool MkSPH::ComputeAccel()
{
  int i;
  char s[256];
  for (i = 0; i < Particles.GetSize(); i++)
  {
    if (Particles[i].GetFixity())
      continue;
    MkParticle &par = Particles[i];
    Particles[i].GetXAccel(2) = Particles[i].GetXAccel(1);
    Particles[i].GetYAccel(2) = Particles[i].GetYAccel(1);
    Particles[i].GetXAccel(1) = Particles[i].GetXAccel(0);
    Particles[i].GetYAccel(1) = Particles[i].GetYAccel(0);
    Particles[i].GetXAccel(0) = (Particles[i].GetXForceInt() + Particles[i].GetXForceVis() + Particles[i].GetXForceArtVis()) / Particles[i].GetMass();
    Particles[i].GetYAccel(0) = (Particles[i].GetYForceInt() + Particles[i].GetYForceVis() + Particles[i].GetYForceArtVis()) / Particles[i].GetMass();
    sprintf(s, "i %d old x acc %f old y acc %f new x acc %f new y acc %f \n", i, par.GetXAccel(1), par.GetYAccel(1), par.GetXAccel(0), par.GetYAccel(0));
    MkDebug(s);
  }
  //getch();
}

bool MkSPH::ComputeVelocity()
{
  int i;
  char s[256];
  //  float beta = 0.25, gamma = 0.5; // trapezoidal
  //  float beta = 0.0, gamma = 0.5; // Verlet
  float beta = 0.0, gamma = 0.25; // central difference
  for (i = 0; i < Particles.GetSize(); i++)
  {
    if (Particles[i].GetFixity())
      continue;
    MkParticle &par = Particles[i];
    par.GetXVelocity(2) = par.GetXVelocity(1);
    par.GetYVelocity(2) = par.GetYVelocity(1);
    par.GetXVelocity(1) = par.GetXVelocity(0);
    par.GetYVelocity(1) = par.GetYVelocity(0);
    //par.GetXVelocity(0) = par.GetXVelocity(0) + DeltaTime*((1-gamma)*par.GetXAccel(1)+gamma*par.GetXAccel(0));
    //par.GetYVelocity(0) = par.GetYVelocity(0) + DeltaTime*((1-gamma)*par.GetYAccel(1)+gamma*par.GetYAccel(0));

    par.GetXVelocity(0) = par.GetXVelocity(1) + DeltaTime * par.GetXAccel(1);
    par.GetYVelocity(0) = par.GetYVelocity(1) + DeltaTime * par.GetYAccel(1);
    sprintf(s, "i %d old x vel %f old y vel %f new x vel %f new y vel %f \n", i, par.GetXVelocity(1), par.GetYVelocity(1), par.GetXVelocity(0), par.GetYVelocity(0));
    MkDebug(s);
  }
  //getch();
}

bool MkSPH::ComputePosition()
{
  int i;
  char s[256];
  //  float beta = 0.25, gamma = 0.5; // trapezoidal
  //  float beta = 0.0, gamma = 0.5; // Verlet
  float beta = 0.0, gamma = 0.25; // central difference
  for (i = 0; i < Particles.GetSize(); i++)
  {
    if (Particles[i].GetFixity())
      continue;
    MkParticle &par = Particles[i];
    par.GetXPosition(2) = par.GetXPosition(1);
    par.GetYPosition(2) = par.GetYPosition(1);
    par.GetXPosition(1) = par.GetXPosition(0);
    par.GetYPosition(1) = par.GetYPosition(0);
    //par.GetXPosition(0) = par.GetXPosition(1) - DeltaTime*par.GetXVelocity(0)+ DeltaTime*DeltaTime/2*((1-2*beta)*par.GetXAccel(1)+2*beta*par.GetXAccel(0));
    //par.GetYPosition(0) = par.GetYPosition(1) - DeltaTime*par.GetYVelocity(0)+ DeltaTime*DeltaTime/2*((1-2*beta)*par.GetYAccel(1)+2*beta*par.GetYAccel(0));
    par.GetXPosition(0) = par.GetXPosition(1) - DeltaTime * par.GetXVelocity(0);
    par.GetYPosition(0) = par.GetYPosition(1) - DeltaTime * par.GetYVelocity(0);
    sprintf(s, "i %d old x pos %f old y pos %f new x pos %f new y pos %f \n", i, par.GetXPosition(1), par.GetYPosition(1), par.GetXPosition(0), par.GetYPosition(0));
    MkDebug(s);
  }
  //getch();
}

bool MkSPH::UpdatePosition()
{
  int i;
  for (i = 0; i < Particles.GetSize(); i++)
  {
    if (Particles[i].GetFixity())
      continue;
    MkParticle &par = Particles[i];
    par.SetPoint(par.GetXPosition(0), par.GetYPosition(0));
  }
}

bool MkSPH::Run() // single step
{
  static bool isInit = false;

  if (!isInit)
  {
    if (NumOfParticles <= 0)
    {
      MkDebug("Please set Particle No. \n");
      return false;
    }
    if (GI * GJ <= 0)
    {
      MkDebug("Please set SPH Grid No. GI and GJ \n");
      return false;
    }
    if ((UpX - LowX) <= 0 && (UpY - LowY) <= 0)
    {
      MkDebug("Please check model domain \n");
      return false;
    }
    if (RefDensity <= 0)
    {
      MkDebug("Plase set Reference Density or Initial Density \n");
      return false;
    }
    if (Gamma_P <= 0)
    {
      MkDebug("Please set Polytropic Constant (Gamma, 1 - 7 )\n");
      return false;
    }
    if (SoundSpeed <= 0)
    {
      MkDebug("Please set speed of sound (could be lower than actual value but 10 times higer than max speed \n");
      return false;
    }
    if (XGravity == 0 && YGravity == 0)
    {
      MkDebug("PLease check that gravity set zero\n");
    }
    isInit = true;
  }
  ComputeDensity();
  //  ComputeDensityRate();
  ComputePressure();
  ComputeInternalForce2();
  //  ComputeViscousForce2();
  ComputeArtificialViscousForce2();
  ComputeAccel();
  ComputeVelocity();
  ComputePosition();
  UpdatePosition();
  MkDebug(".");
  Time += DeltaTime;
  //getch();
}
//---------------------------------------------------------------------------
#pragma package(smart_init)
