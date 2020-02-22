//---------------------------------------------------------------------------
#pragma hdrstop
#include "MkLiuTest.hpp"

//---------------------------------------------------------------------------
GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};  /* Gray diffuse light. */
GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0}; /* Infinite light location. */

GLfloat v[8][3]; /* Will be filled in with X,Y,Z vertexes. */

MkLiuParticles Particles;
//MkLiuGrid LiuGrids;
MkLiuSPH mySPH;

bool shear_cavity = true;
bool shock_tube = false;

void initGL(void)
{
  /* Enable a single OpenGL light. */
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL); // Enable color

  /* Use depth buffering for hidden surface elimination. */
  glEnable(GL_DEPTH_TEST);

  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective(/* field of view in degree */ 40.0,
                 /* aspect ratio */ 1.0,
                 /* Z near */ 1.0, /* Z far */ 30.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(-0.0, 0.5, 5.0, /* eye is at (0,0,5) */
            -0.0, 0.5, 0.0, /* center is at (0,0,0) */
            0.0, 1.0, 0.);  /* up is in positive Y direction */

  glTranslatef(0.0, 0.0, 1.0);
}

void InitLiuSPH(void)
{
  int ntotal = 0, nvirt = 0, gi = 100, gj = 100, gk = 1;
  int m, n, mp, np;
  double rho_ref = 1000;

  if (mySPH.GetLiuParam().shocktube)
  {
    ntotal = 400;
    nvirt = 100;

    mySPH.Initialize(ntotal, nvirt, gi, gj, gk);
    mySPH.SetDt(0.005);

    mySPH.Shock_Tube();
  }
  else if (mySPH.GetLiuParam().shearcavity)
  {
    m = 21;
    n = 21;
    mp = m - 1;
    np = n - 1;
    ntotal = mp * np;
    nvirt = 2 * (2 * mp + 1) + 2 * (2 * np - 1);

    MkDebug("nvirt = %d\n", nvirt);
    mySPH.Initialize(ntotal, nvirt, gi, gj, gk);

    mySPH.SetDt(5.e-3);
    mySPH.SetDim(3);
    mySPH.SetMaxTimeStep(3000);
    MkDebug("ntotal = %d\n", ntotal);

    mySPH.SetRhoRef(rho_ref);

    mySPH.Shear_Cavity();

    MkDebug("ntotal = %d\n", ntotal);
    mySPH.Virt_Part();
    mySPH.LiuParam.moni_particle = 0;
    mySPH.LiuParticles.SetColorType(ctRho);

    MkLiuParticles &par = mySPH.GetLiuParticles();
    int moni = 0;

    MkDebug("Moni_particle\n");
    MkDebug("x:%f   y:%f  \n", par[moni].X, par[moni].Y);
  }
  else if (mySPH.GetLiuParam().lowdensity)
  {
    MkDebug("MkLiuTest::InitSPH()\n");
    m = 50;
    n = 50;
    mp = m - 1;
    np = n - 1;
    ntotal = mp * np;
    nvirt = 0; // 2*(2*mp+1) + 2*(2*np-1);
    rho_ref = 1000.0;

    MkDebug("nvirt = %d\n", nvirt);
    mySPH.Initialize(ntotal, nvirt, gi, gj, gk);
    MkLiuBoundarys &bnd = mySPH.LiuBoundarys;
    bnd.Initialize(4);
    bnd[0].SetPoints(MkPoint(-1.7, 0.0, 0.0), MkPoint(-1.7, 0.0, 1.0), MkPoint(-1.7, 1.0, 0.0));
    bnd[1].SetPoints(MkPoint(0.0, -0.0, 0.0), MkPoint(1.0, -0.0, 0.0), MkPoint(0.0, -0.0, 1.0));
    bnd[2].SetPoints(MkPoint(1.0, 0.0, 0.0), MkPoint(1.0, 1.0, 0.0), MkPoint(1.0, 0.0, 1.0));
    bnd[3].SetPoints(MkPoint(0.0, 2.0, 0.0), MkPoint(0.0, 2.0, 1.0), MkPoint(1.0, 2.0, 0.0));
    // MkLiuGrids &grids = mySPH.LiuGrids;
    // grids.SetXMin(bnd[0].GetCenter().X);
    // grids.SetYMin(bnd[1].GetCenter().Y);
    // grids.SetZMin(-0.01);
    // grids.SetXMax(bnd[2].GetCenter().X);
    // grids.SetYMax(bnd[3].GetCenter().Y);
    // grids.SetZMax(0.01);

    mySPH.SetRhoRef(rho_ref);

    mySPH.SetDt(5.0e-3);
    mySPH.SetDim(3);
    mySPH.SetMaxTimeStep(50000);
    MkDebug("ntotal = %d\n", ntotal);
    mySPH.LowDensity();
    mySPH.Shake();
    MkDebug("ntotal = %d\n", ntotal);
    //mySPH.Virt_LD();
    mySPH.LiuParam.moni_particle = 0;
    mySPH.LiuParticles.SetColorType(ctRho);
    mySPH.SetRhoRef(rho_ref);

    MkLiuParticles &par = mySPH.GetLiuParticles();
    int moni = mySPH.LiuParam.moni_particle;

    MkDebug("Moni_particle\n");
    MkDebug("x:%f   y:%f  \n", par[moni].X, par[moni].Y);
    MkDebug("MkLiuTest::~InitSPH()\n");
  }
}

void twoparticles()
{
  int i, j, k, ntotal = 0, nvirt = 0, gi = 100, gj = 100, gk = 1;
  int m, n, mp, np;
  double rho_ref;

  m = 21;
  mp = m - 1;
  n = 21;
  np = n - 1;
  ntotal = mp * np;
  nvirt = 0;
  rho_ref = 1200.0;

  MkDebug("nvirt = %d\n", nvirt);
  mySPH.Initialize(ntotal, nvirt, gi, gj, gk);

  // MkLiuBoundarys &bnd = mySPH.LiuBoundarys;
  // bnd.Initialize(1);
  // bnd[0].SetPoints(MkPoint(-1.0, -10.0, 0.0), MkPoint(1., -10.0, 0.0), MkPoint(0.0, -10.0, 1.0));

  mySPH.SetDt(1.e-2);
  mySPH.SetDim(3);
  mySPH.SetMaxTimeStep(9999999);
  MkDebug("ntotal = %d\n", ntotal);

  mySPH.LiuParam.moni_particle = 0;
  mySPH.LiuParticles.SetColorType(ctRho);
  mySPH.SetRhoRef(rho_ref);

  MkLiuParticles &par = mySPH.GetLiuParticles();
  int moni = mySPH.LiuParam.moni_particle;
  double dx = 1.0 / (mp - 1);

  for (i = 0; i < mp; i++)
  {
    for (j = 0; j < np; j++)
    {
      k = j + i * np;
      par[k].X = 0.0 + dx * 1 * i;
      par[k].Y = 0.0 + dx * 1 * j;
    }
  }

  for (i = 0; i < ntotal; i++)
  {
    par[i].XVel = 0.;
    par[i].YVel = 0.;
    par[i].Rho = rho_ref;
    par[i].Mass = dx * dx * par[i].Rho;
    par[i].Eta = 10000000000;
    par[i].Press = 0.0;
    par[i].Energy = 357.1;
    par[i].ParticleType = 2;
    par[i].SmoothLen = dx;
    par[i].Radius = dx / 2.0;
  }

  MkDebug("Moni_particle %d\n", moni);
  MkDebug("x:%f   y:%f  \n", par[moni].X, par[moni].Y);
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();
  glRotatef(180, 0., 0., 1.);

  mySPH.Draw();

  glPopMatrix();
  glutSwapBuffers();
}

void update(void)
{
  /*
  static int cnt=0;
  double vel=0.0001/8.e-3;
  MkLiuParticles &par = mySPH.GetLiuParticles();
  //  par[0].X += 0.001;
  //  par[0].Y += 0.001;
  for (int i=0;i<9;i++) {
    par[i].X = 0.07+0.0001*cnt;
    par[i].Y = (i+0.25)*(1/9.0);
    //    par[i].Y += cnt*0.0002;
    //    par[i].YVel[0] = vel;
  }

  for (int i=9;i<18;i++) {
    par[i].X = 1.0-0.0001*cnt;;
    par[i].Y = (i+0.5-9)*(1/9.0);
    //    par[i].Y += -cnt*0.0002;
    //    par[i].YVel[0] = -vel;
  }
    cnt++;
  */
  /*
  static int cnt=1;
  double vel=0.001/mySPH.Dt;
  MkLiuParticles &par = mySPH.GetLiuParticles();
  if (cnt > 50) {
    par[0].X += 0.001;
    par[0].Y += 0.001;
    par[0].XVel = vel;
    par[0].YVel = vel;
  }*/
  /*
  if (cnt > 100) {
    par[0].X = cos((cnt%1000)/500.0*3.14159)/2.0+0.5;
    par[0].Y = sin((cnt%1000)/500.0*3.14159)/2.0+0.5;
    //    par[0].XVel[0] = vel;
    //par[0].YVel[0] = vel;
  }
  */
  //  cnt++;

  //  mySPH.NTotal = (cnt<19*19) ? cnt:19*19;
  mySPH.Run();

  glutPostRedisplay(); //Tell GLUT that the scene has changed
}

void update2(void)
{
  static bool first = true;
  static int n = 1;

  srand(2);
  //  if(first)
  if (n > 100)
    exit(-1);
  if (n % 3 == 0)
  {
    for (int i = 0; i < mySPH.GetLiuParticles().GetSize(); i++)
      mySPH.GetLiuParticles()[i].SetRho(1000);
    mySPH.Direct_Find();
    n = 0;
    first = false;
  }
  mySPH.Sum_Density();
  mySPH.Int_Force();
  //usleep(999999);usleep(999999);usleep(999999);usleep(999999);usleep(999999);
  glutPostRedisplay();
  n++;
}

/*
int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_ALPHA);
    glutCreateWindow("Liu SPH Test");
    glutDisplayFunc(display);

    initGL();
    InitLiuSPH();
    
    glutIdleFunc(update);
    glutMainLoop();

  return 0;
}
*/

void kernelcheck()
{
  MkLiuKernel k1;
  float sml = 0.5, sum = 0;
  int dim = 1;

  k1.SetupLiuKernel(knlGaussian, dim, sml);
  //k1.SetupLiuKernel(knlQuintic, dim, sml);
  //k1.SetupLiuKernel(knlCubicSpline, dim, sml);

  for (int i = -200; i < 200; i++)
  {
    float div = 100;
    float x = sml * i / div;
    float r = fabs(x / sml);
    float dr = 1 / div;
    float dx = fabs(0.3 * r), dy = fabs(0.3 * r), dz = sqrt(r * r - dx * dx - dy * dy);
    sum += ((dim == 3) ? 2 : 1) * ((dim == 1) ? 1 : 3.14159) * pow(r, dim - 1) * dr * k1.W(r);
    printf("%f %f %f %f %f %f %f\n", x, r, k1.W(r), k1.dWdR(r), k1.dWdX(r, dx, dy, dz), k1.dWdY(r, dx, dy, dz), k1.dWdZ(r, dx, dy, dz));
  }
  printf("sum is %f\n", sum);
}

void densitycheck(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_ALPHA);
  glutCreateWindow("Liu SPH Test");
  glutDisplayFunc(display);

  initGL();
  InitLiuSPH();

  glutIdleFunc(update);
  glutMainLoop();
}

void tpcheck(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_ALPHA);
  glutCreateWindow("Liu SPH Test");
  glutDisplayFunc(display);

  initGL();
  //InitLiuSPH();
  twoparticles();

  glutIdleFunc(update);
  glutMainLoop();
}

void bndcheck(int argc, char **argv)
{
  double rho_norm;

  InitLiuSPH();

  rho_norm = mySPH.Calc_Bnd_Norm();

  MkDebug("Rho_Norm is %f\n", rho_norm);
}

int main(int argc, char **argv)
{
  //densitycheck(argc, argv);
  //bndcheck(argc,argv);
  // boost::shared_ptr<int> int_ptr = boost::make_shared<int>(5);

  // MkInt iarr(5);
  // MkDouble darr(5);
  // MkArcs arc(5);
  // MkCubes cube(5);
  // MkCircles cir(5);
  // MkCylinders cyl(5);
  // MkLines lines(5);
  // MkPlanes planes(5);
  // MkJointPlanes jplanes(5);
  // MkPennyJoints pjs(5);
  // MkPolygons polys(5);
  // MkRects rects(5);
  // MkSpheres sphs(5);
  // MkTriangles tris(5);

  // try
  // {
  //   std::cout << iarr(5) << std::endl;
  // }
  // catch (MkInt::Range &r)
  // {
  //   MkDebug("MkInt::Range Error is captured %s, %d \n", r.what(), r.N);
  // }

  // std::cout
  //     << iarr.getSzX() << std::endl;

  // for (int i = 0; i < iarr.getSzX(); i++)
  // {
  //   iarr(i) = i;
  //   darr(i) = i;
  //   printf("i = %d d = %f\n", iarr(i), darr(i));
  //   std::cout << arc[i].ClassName() << std::endl;
  //   std::cout << cube[i].ClassName() << std::endl;
  //   std::cout << cir[i].ClassName() << std::endl;
  //   std::cout << cyl[i].ClassName() << std::endl;
  //   std::cout << lines[i].ClassName() << std::endl;
  //   std::cout << planes[i].ClassName() << std::endl;
  //   std::cout << jplanes[i].ClassName() << std::endl;
  //   std::cout << pjs[i].ClassName() << std::endl;
  //   std::cout << polys[i].ClassName() << std::endl;
  //   std::cout << rects[i].ClassName() << std::endl;
  //   std::cout << sphs[i].ClassName() << std::endl;
  //   std::cout << tris[i].ClassName() << std::endl;
  // }

  // std::cout << "\nMkLiuTest is performing\n"
  //           << *int_ptr
  //           << std::endl;

  // glutInit(&argc, argv);
  // MkDebug("1\n");
  // glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_ALPHA);
  // MkDebug("2\n");
  // glutCreateWindow("Liu SPH Test");
  // MkDebug("3\n");
  // glutDisplayFunc(display);
  // MkDebug("4\n");
  // initGL();
  // MkDebug("5\n");
  // mySPH.GetLiuParam().shocktube = false;
  // mySPH.GetLiuParam().lowdensity = true;
  // std::cout << "mySPH.GetLiuParam().shocktube " << mySPH.GetLiuParam().shocktube << std::endl;
  // InitLiuSPH();
  // MkDebug("6\n");
  // glutIdleFunc(update);
  // MkDebug("7\n");
  // glutMainLoop();
  // MkDebug("8\n");

  tpcheck(argc, argv);
  return 0;
}

#pragma package(smart_init)
