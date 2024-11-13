#include "MkLiuSPH_mk2.hpp"

MkLiuSPH::MkLiuSPH()
{
  Clear();
}

MkLiuSPH::~MkLiuSPH()
{
  Clear();
}

bool MkLiuSPH::Initialize()
{
  Clear();  
}

bool MkLiuSPH::Initialize(int np, int nvirt, int gx, int gy, int gz)
{
  NTotal = np;
  NVirt = nvirt;

  LiuParticles.Initialize(NTotal+NVirt);
  LiuKernel.SetKernelType(knlCubicSpline);
  LiuKernel.SetDim(3);
  LiuKernel.SetSmoothLen(0.001/40.0);
}

void MkLiuSPH::Clear()
{
//MkDebug("MkLiuSPH::Clear()\n");

  Dim = 0;
  NTotal=0;
  NVirt=0;
      
  MaxInteraction = 0;
  NIac=0;        

  MaxTimeStep=0; 
  CurrentTimeStep=0;   
  Dt=0;       
  Rho_Ref = 0;
  Rho_Norm = 0;

  LiuParticles.Clear();
  LiuBoundarys.Clear();
  LiuKernel.Clear();
  LiuPairs.Clear();    
  LiuBndPairs.Clear();
//MkDebug("MkLiuSPH::~Clear()\n");
}

void MkLiuSPH::Draw()
{
  LiuParticles.Draw();
  LiuBoundarys.Draw();
} 
 
void MkLiuSPH::Run()
{
  static int i, j, k, d;
  static double time, temp_rho, temp_u;
  static bool is_first = true;
  MkLiuParticles &par=LiuParticles;

  int maxn=NTotal+NVirt;
//MkDebug("MkLiuSPH::Run() Dim %d, maxn %d\n",Dim, maxn);

  MkDouble v_min(Dim, maxn), u_min(maxn), rho_min(maxn); 

//MkDebug("MkLiuSPH::Run()\n");
  if(is_first) {
    i = j = k = d = 0;
    time = temp_rho = temp_u = 0;

    for (i = 0;i< NTotal;i++) {
      par[i].XVel = 0;
      par[i].YVel = 0;
      par[i].ZVel = 0;
    }
  }
//MkDebug("MkLiuSPH::Run() #1 \n");
  if (CurrentTimeStep> MaxTimeStep) {is_first = false; exit(-1); return;}
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    printf("______________________________________________\n");
    printf("  current number of time step = %d     current time=%f\n",CurrentTimeStep, double(time+Dt));
    printf("______________________________________________\n");
  }
//MkDebug("MkLiuSPH::Run() #2 \n");

//   If not first time step, then update thermal energy, density and  
//   velocity half a time step   

  if (!is_first) {
    for (i = 0;i< NTotal;i++) {
      u_min(i) = par[i].Energy;
      temp_u=0.;
      if (Dim==1) temp_u=-LiuParam.nsym*par[i].Press*par[i].XVel/par[i].X/par[i].Rho;
      par[i].Energy = par[i].Energy + (Dt/2.)* (par[i].DUDt+temp_u);
      if(par[i].Energy < 0) par[i].Energy = 0.;
             
      if (!LiuParam.summation_density) {
	rho_min(i) = par[i].Rho;
	temp_rho=0.;
	if (Dim==1) temp_rho=-LiuParam.nsym*par[i].Rho*par[i].XVel/par[i].X;
	par[i].Rho = par[i].Rho +(Dt/2.)*( par[i].DRhoDt+ temp_rho);
      }

      v_min(0, i) = par[i].XVel;
      par[i].XVel = par[i].XVel + (Dt/2.)*par[i].DVXDt;
      v_min(1, i) = par[i].YVel;
      par[i].YVel = par[i].YVel + (Dt/2.)*par[i].DVYDt;
      v_min(2, i) = par[i].ZVel;
      par[i].ZVel = par[i].ZVel + (Dt/2.)*par[i].DVZDt;
    }
  }
//---  Definition of variables out of the function vector:     

  //MkDebug("MkLiuSPH::Run() #3 \n");
  Single_Step();
  //MkDebug("MkLiuSPH::Run() #4 \n");  
  if (is_first) {
   
    for (i=0;i<NTotal;i++) {
      temp_u=0.;
      if (Dim==1) temp_u=-LiuParam.nsym*par[i].Press*par[i].XVel/par[i].X/par[i].Rho;
      par[i].Energy = par[i].Energy + (Dt/2.)*(par[i].DUDt + temp_u);
      if(par[i].Energy<0)  par[i].Energy = 0.;
          
      if (!LiuParam.summation_density ){
	temp_rho=0.;
	if (Dim==1) temp_rho=-LiuParam.nsym*par[i].Rho*par[i].XVel/par[i].X;
	par[i].Rho = par[i].Rho + (Dt/2.)* (par[i].DRhoDt+temp_rho);
      }
          
      par[i].XVel = par[i].XVel + (Dt/2.) * par[i].DVXDt + par[i].XAVel;
      par[i].X = par[i].X + Dt * par[i].XVel;
      par[i].YVel = par[i].YVel + (Dt/2.) * par[i].DVYDt + par[i].YAVel;
      par[i].Y = par[i].Y + Dt * par[i].YVel;
      par[i].ZVel = par[i].ZVel + (Dt/2.) * par[i].DVZDt + par[i].ZAVel;
      par[i].Z = par[i].Z + Dt * par[i].ZVel;

    }
  }               
  else {
    for (i=0;i<NTotal;i++) {
      temp_u=0.;
      if (Dim==1) temp_u=-LiuParam.nsym*par[i].Press*par[i].XVel/par[i].X/par[i].Rho;
      par[i].Energy = u_min(i) + Dt*(par[i].DUDt+temp_u);
      if(par[i].Energy<0)  par[i].Energy = 0.;
             
      if (!LiuParam.summation_density ) {
	temp_rho=0.;
	if (Dim==1) temp_rho=-LiuParam.nsym*par[i].Rho*par[i].XVel/par[i].X;
	par[i].Rho = rho_min(i) + Dt*(par[i].DRhoDt+temp_rho);
      }
                 
      par[i].XVel = v_min(0, i) + Dt * par[i].DVXDt + par[i].XAVel;
      par[i].X = par[i].X + Dt * par[i].XVel;
      par[i].YVel = v_min(1, i) + Dt * par[i].DVYDt + par[i].YAVel;
      par[i].Y = par[i].Y + Dt * par[i].YVel;
      par[i].ZVel = v_min(2, i) + Dt * par[i].DVZDt + par[i].ZAVel;
      par[i].Z = par[i].Z + Dt * par[i].ZVel;
    }
         
  }
//MkDebug("MkLiuSPH::Run() #5 \n"); 
  if ((CurrentTimeStep%LiuParam.save_step)==0)  Output();
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    printf("\n");
    //      123456789ABC 123456789ABC 123456789ABC 123456789ABC 123456789ABC 123456789ABC 123456789ABC 123456789ABC 123456789ABC  
    printf("      x           y            Rho         Press         xvel        yvel           dvx         dvy          dvz\n");
    for (i=0;i<NTotal;i++) 
    printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
	   par[i].X, 
	   par[i].Y,
	   par[i].Rho,
	   par[i].Press,
	   par[i].XVel,
	   par[i].YVel,
	   par[i].DVXDt,
	   par[i].DVYDt,
	   par[i].DVZDt);
    }

  time = time + Dt;
  CurrentTimeStep++;
  is_first = false;
//MkDebug("~MkLiuSPH::Run()\n");
  return;
}

//---------------------------------------------------------------------- 
// Subroutine to determine the right hand side of a differential  
// equation in a single step for performing time integration  
 
// In this routine and its subroutines the SPH algorithms are performed. 
//   CurrentTimeStep: Current timestep number                            [in] 
//   Dt       : Timestep                                           [in] 
//   NTotal   :  Number of particles                               [in] 
//   Hsml     :  Smoothing Length                                  [in] 
//   Mass     :  Particle masses                                   [in] 
//   X        :  Particle position                                 [in] 
//   VX       :  Particle velocity                                 [in] 
//   U        :  Particle internal energy                          [in] 
//   Rho      :  Density                                       [in/out] 
//   p        :  Pressure                                         [out] 
//   t        :  Temperature                                   [in/out] 
//   TDSDt    :  Production of viscous entropy t*ds/dt            [out] 
//   dx       :  dx = VX = dx/dt                                  [out] 
//   dvx      :  dvx = dvx/dt, force per unit mass                [out] 
//   DUDt       :  du  = du/dt                                      [out] 
//   DSDt       :  ds  = ds/dt                                      [out]      ???
//   DRhoDt   :  drhodt =  drho/dt                                  [out] 
//   ParticleType    :  Type of particle                                 [in] 
//   av       :  Monaghan average velocity                        [out] 
 
void MkLiuSPH::Single_Step()  
{
  int i, d;
  int maxn = NTotal+NVirt;
  MkLiuParticles &par = LiuParticles;
 
//MkDebug("\n  single_step()\n");
//MkDebug("  S1 ");
  for ( i=0;i<NTotal;i++) {
    par[i].AVDUDt = 0.; 
    par[i].AHDUDt = 0.; 
    par[i].ADDUDt = 0.; 
    par[i].APDUDt = 0.; 
    par[i].INDVXDt = 0.;
    par[i].ARDVXDt = 0.; 
    par[i].ADDVXDt = 0.; 
    par[i].APDVXDt = 0.; 
    par[i].EXDVXDt = 0.;
    par[i].INDVYDt = 0.;
    par[i].ARDVYDt = 0.; 
    par[i].ADDVYDt = 0.; 
    par[i].APDVYDt = 0.; 
    par[i].EXDVYDt = 0.;
    par[i].INDVZDt = 0.;
    par[i].ARDVZDt = 0.; 
    par[i].ADDVZDt = 0.; 
    par[i].APDVZDt = 0.; 
    par[i].EXDVZDt = 0.;
  }
//MkDebug("  S2 ");  
//---  Positions of virtual (boundary) particles:  
 
  if (LiuParam.virtual_part) {
    //Virt_Part(); 
    //Virt_LD(); 
  }
//MkDebug("  S3 ");      
//---  Interaction parameters, calculating neighboring particles 
//   and optimzing smoothing length 
 
  if (LiuParam.nnps==1) { 
    Direct_Find();
    Direct_Bnd_Find();
  }

  //MkDebug("  S4 ");
//---  Density approximation or change rate 

  if (LiuParam.summation_density) Norm_Bnd_Density();//Uni_Density(); //MK_Density();//Sum_Density();
  else Con_Density();
 
//---  Dynamic viscosity: 
 
//MkDebug("  S5 ");
  if (LiuParam.visc) Viscosity(); 
        
//---  Internal forces: 
  
  Int_Force(); 
  //Int_Bnd_Force();
               
//MkDebug("  S6 ");    

//---  Artificial viscosity: 
 
  //Art_Drag(); // artificial dragging similar to artificial viscocity 

  Art_Bnd_Repel(); // artificial repelling to give minimum spacing of particles

  if (LiuParam.visc_artificial) 
    Art_Visc(); //not working...

//MkDebug("  S7 ");
       
//---  External forces: 
 
//MkDebug("  S8 ");

  if (LiuParam.ex_force) {Ext_Force();}
 
//MkDebug("  S9 ");

//   Calculating the neighboring particles and undating DATA.HSML 
       
  if (LiuParam.sle!=0) H_Upgrade(); 
 
  if (LiuParam.heat_artificial) Art_Heat();
      
//   Calculating average velocity of each partile for avoiding penetration 
 
  if (LiuParam.average_velocity) Av_Vel();
 
//---  Convert velocity, force, and energy to f and dfdt   
 
  for (i=0;i<NTotal;i++) {
    par[i].DVXDt=par[i].INDVXDt+par[i].EXDVXDt+par[i].ARDVXDt+par[i].ADDVXDt+par[i].APDVXDt; 
    par[i].DVYDt=par[i].INDVYDt+par[i].EXDVYDt+par[i].ARDVYDt+par[i].ADDVYDt+par[i].APDVYDt; 
    par[i].DVZDt=par[i].INDVZDt+par[i].EXDVZDt+par[i].ARDVZDt+par[i].ADDVZDt+par[i].APDVZDt;
    par[i].DUDt =par[i].DUDt+par[i].INDUDt+par[i].AVDUDt+par[i].AHDUDt+par[i].APDUDt;
  }
  // //MkDebug("  S10 ");
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    //      123456789abc 123456789abc 123456789abc 123456789abc 123456789abc
    printf("\n") ;
    printf("**** Information for particle **** %d\n",LiuParam.moni_particle);
    printf("    xvel      internal a   artifical a     art drag   art repel    external a     total a \n");
    printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
	   par[LiuParam.moni_particle].XVel,
	   par[LiuParam.moni_particle].INDVXDt,
	   par[LiuParam.moni_particle].ARDVXDt,
	   par[LiuParam.moni_particle].ADDVXDt,
	   par[LiuParam.moni_particle].APDVXDt,
	   par[LiuParam.moni_particle].EXDVXDt,
	   par[LiuParam.moni_particle].DVXDt);
  }

//MkDebug("  ~single_step()\n");
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the artificial heat(Fulk, 1994, p, a-17)  
//  See Equ.(4.74)   
 
//  NTotal : Number of particles                                  [in] 
//  Hsml   : Smoothing Length                                     [in] 
//  Mass   : Particle masses                                      [in] 
//  X      : Coordinates of all particles                         [in] 
//  VX     : Velocities of all particles                          [in] 
//  Rho    : Density                                              [in] 
//  U      : specific internal energy                             [in] 
//  C      : Sound veolcity                                       [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  W      : Kernel for all interaction pairs                     [in] 
//  DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//  AHDUDt   : produced artificial heat, adding to energy Eq.      [out] 

void MkLiuSPH::Art_Heat()
{ 
  int i,j,k,d ;
  double dx, vr, rr, h, mc, mrho, mhsml, hvcc, mui, muj, muij, rdwdx, g1,g2;
  MkDouble vcc(NTotal);
  double dvx[Dim];
  MkLiuParticles &par = LiuParticles;

  i=j=k=d =0;
  dx= vr= rr= h= mc= mrho= mhsml= hvcc= mui= muj= muij= rdwdx= g1=g2=0;
       
  //---  Parameter for the artificial heat conduction:
      
  g1=0.1; 
  g2=1.0;
  for (i=0;i<NTotal+NVirt;i++) { 
    vcc(i) = 0.e0;
    par[i].AHDUDt = 0.e0; 
  }
     
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 

    dvx[0] = par[j].XVel - par[i].XVel;
    dvx[1] = par[j].YVel - par[i].YVel;
    dvx[2] = par[j].ZVel - par[i].ZVel;

    hvcc = dvx[0]*LiuPairs[k].DWDX(1-1);

    for (d=1;d<Dim;d++) { 
      hvcc = hvcc + dvx[d]*LiuPairs[k].DWDX(d) ;
    }    
    vcc(i) = vcc(i) + par[j].Mass*hvcc/par[j].Rho; 
    vcc(j) = vcc(j) + par[i].Mass*hvcc/par[i].Rho; 
  }
    
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    mhsml= (par[i].SmoothLen+par[j].SmoothLen)/2.; 
    LiuKernel.SetSmoothLen(mhsml);

    mrho = 0.5e0*(par[i].Rho + par[j].Rho);
    rr = 0.e0; 
    rdwdx = 0.e0; 

    dx = par[i].X -  par[j].X;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*LiuPairs[k].DWDX(d-1);
    dx = par[i].Y -  par[j].Y;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*LiuPairs[k].DWDX(d-1);
    dx = par[i].Z -  par[j].Z;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*LiuPairs[k].DWDX(d-1);


    mui=g1*par[i].SmoothLen*par[i].SoundSpeed + 
        g2*par[i].SmoothLen*par[i].SmoothLen*(fabs(vcc(i))-vcc(i));
    muj=g1*par[j].SmoothLen*par[j].SoundSpeed + 
        g2*par[j].SmoothLen*par[j].SmoothLen*(fabs(vcc(j))-vcc(j));
    muij= 0.5*(mui+muj);
    h = muij/(mrho*(rr+0.01*mhsml*mhsml))*rdwdx ;
    par[i].AHDUDt = par[i].AHDUDt + par[j].Mass*h*(par[i].Energy-par[j].Energy) ;
    par[j].AHDUDt = par[j].AHDUDt + par[i].Mass*h*(par[j].Energy-par[i].Energy); 
  }
 
  for (i=0;i<NTotal+NVirt;i++) { 
    par[i].AHDUDt = 2.0e0*par[i].AHDUDt;           
  //MkDebug("%d, %12.6f\n",i,par[i].AHDUDt);
  }
}

 
//---------------------------------------------------------------------- 
//  Subroutine to calculate the artificial viscosity (Monaghan, 1992)  
//  See Equ.(4.66) Equ.(4.62) 
 
//  NTotal : Number of particles (including virtual particles)    [in] 
//  Hsml   : Smoothing Length                                     [in] 
//  Mass   : Particle masses                                      [in] 
//  X      : Coordinates of all particles                         [in] 
//  VX     : Velocities of all particles                          [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Rho    : Density                                              [in] 
//  C      : Temperature                                     [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  W      : Kernel for all interaction pairs                     [in] 
//  DWDX   : Derivative of kernel( with respect to x, y and z     [in] 
//  ARDVXDt  : Acceleration with respect to x, y and z            [out]  
//  AVDUDt   : Change of specific internal energy                 [out] 

void MkLiuSPH::Art_Visc()
{
  int i,j,k,d ;
  double dx, alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml;
  double dvx[Dim];
  MkLiuParticles &par = LiuParticles;

  //MkDebug("art_visc()\n");
   i=j=k=d =0;
   dx= alpha= beta= etq= piv= muv= vr= rr= h= mc= mrho= mhsml=0;

//  Parameter for the artificial viscosity: 
//  Shear viscosity 
  alpha = 1.e-4; 
      
//  Bulk viscosity 
  beta  = 1.e-4;  
       
//  Parameter to avoid singularities 
  etq   = 0.1e0;
            
  for (i=0;i<NTotal+NVirt;i++){ 
    par[i].ARDVXDt = 0.e0;
    par[i].ARDVYDt = 0.e0;
    par[i].ARDVZDt = 0.e0;
    par[i].AVDUDt = 0.e0;
  }
      
//  Calculate SPH sum for artificial viscosity 
       
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    mhsml= (par[i].SmoothLen+par[j].SmoothLen)/2.; 
    LiuKernel.SetSmoothLen(mhsml);
    vr = 0.e0; 
    rr = 0.e0; 

    dvx[0] = par[i].XVel - par[j].XVel; 
    dx     = par[i].X - par[j].X; 
    vr     = vr + dvx[0]*dx; 
    rr     = rr + dx*dx; 
    dvx[1] = par[i].YVel - par[j].YVel; 
    dx     = par[i].Y - par[j].Y; 
    vr     = vr + dvx[1]*dx; 
    rr     = rr + dx*dx; 
    dvx[2] = par[i].ZVel - par[j].ZVel; 
    dx     = par[i].Z - par[j].Z; 
    vr     = vr + dvx[2]*dx; 
    rr     = rr + dx*dx; 

//  Artificial viscous force only if v_ij * r_ij  0 
 
    if (vr<0.e0) { 
 
//  Calculate muv_ij = data.Hsml v_ij * r_ij / ( r_ij^2 + data.Hsml^2 etq^2 ) 
             
      muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq);
 
//  Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij 
 
      mc   = 0.5e0*(par[i].SoundSpeed + par[j].SoundSpeed);
      mrho = 0.5e0*(par[i].Rho + par[j].Rho); 
      piv  = (beta*muv - alpha*mc)*muv/mrho;

//  Calculate SPH sum for artificial viscous force 
 
      h = -piv*LiuPairs[k].DWDX(1-1);

      par[i].ARDVXDt = par[i].ARDVXDt + par[j].Mass*h;
      par[j].ARDVXDt = par[j].ARDVXDt - par[i].Mass*h;
      par[i].AVDUDt = par[i].AVDUDt - par[j].Mass*dvx[0]*h;
      par[j].AVDUDt = par[j].AVDUDt - par[i].Mass*dvx[0]*h;
    
      //MkDebug("i %d muv %f mc %f mrho %f piv %f h %f mass %f ardvxdt %f\n",i,muv,mc,mrho,piv,h,par[i].Mass,par[i].ARDVXDt); 

      h = -piv*LiuPairs[k].DWDX(2-1);
      par[i].ARDVYDt = par[i].ARDVYDt + par[j].Mass*h;
      par[j].ARDVYDt = par[j].ARDVYDt - par[i].Mass*h;
      par[i].AVDUDt = par[i].AVDUDt - par[j].Mass*dvx[1]*h;
      par[j].AVDUDt = par[j].AVDUDt - par[i].Mass*dvx[1]*h;
      h = -piv*LiuPairs[k].DWDX(3-1);
      par[i].ARDVZDt = par[i].ARDVZDt + par[j].Mass*h;
      par[j].ARDVZDt = par[j].ARDVZDt - par[i].Mass*h;
      par[i].AVDUDt = par[i].AVDUDt - par[j].Mass*dvx[2]*h;
      par[j].AVDUDt = par[j].AVDUDt - par[i].Mass*dvx[2]*h;



    }
  }
 
//  Change of specific internal energy: 
 
  for (i=0;i<NTotal+NVirt;i++) { 
    par[i].AVDUDt = 0.5e0*par[i].AVDUDt;
  //MkDebug("i:%d, ARDVXDt:%12.6f",i,par[i].ARDVXDt );
  }
  
  //usleep(999999);usleep(999999);usleep(999999);
  
  //MkDebug("~art_visc()\n");
}

void MkLiuSPH::Art_Drag()
{
  int i;
  double radius=0.01;
  MkLiuParticles &par = LiuParticles;
           
  for (i=0;i<NTotal;i++){ 
    par[i].ADDVXDt = 0.e0;
    par[i].ADDVYDt = 0.e0;
    par[i].ADDVZDt = 0.e0;
  }

  for (i=0;i<NTotal;i++){ 
    par[i].ADDVXDt = -6*M_PI*par[i].Eta*radius*par[i].XVel;
    par[i].ADDVYDt = -6*M_PI*par[i].Eta*radius*par[i].YVel;
    par[i].ADDVZDt = -6*M_PI*par[i].Eta*radius*par[i].ZVel;
  }

  // require treatment of internal energy how???
}

void MkLiuSPH::Art_Repel()
{
  int i,j,k;
  double dist, dx,dy,dz,repel;
  MkLiuParticles &par = LiuParticles;

  for (i=0;i<NTotal+NVirt;i++){ 
    par[i].APDVXDt = 0.e0;
    par[i].APDVYDt = 0.e0;
    par[i].APDVZDt = 0.e0;
    par[i].AVDUDt = 0.e0;
  }

  for (k=0;k<NIac;k++) {
    i = LiuPairs[k].I;
    j = LiuPairs[k].J;

    dx = LiuParticles[i].X - LiuParticles[j].X;
    dy = LiuParticles[i].Y - LiuParticles[j].Y; 
    dz = LiuParticles[i].Z - LiuParticles[j].Z; 
    dist = sqrt(dx*dx+dy*dy+dz*dz);
    
    if (dist > par[i].Radius + par[j].Radius) continue;
    repel = 50*(exp(-dist) - exp(-(par[i].Radius + par[j].Radius)));

    //MkDebug("Repelling %d %d with repel %f\n",i,j,repel);
    par[i].APDVXDt = par[i].APDVXDt + repel*dx/dist;
    par[j].APDVXDt = par[j].APDVXDt - repel*dx/dist;
    par[i].APDVYDt = par[i].APDVYDt + repel*dy/dist;
    par[j].APDVYDt = par[j].APDVYDt - repel*dy/dist;
    par[i].APDVZDt = par[i].APDVZDt + repel*dz/dist;
    par[j].APDVZDt = par[j].APDVZDt - repel*dz/dist;
  }

  // require to handle the internal energy due to repel

}

void MkLiuSPH::Art_Bnd_Repel()
{
  int i,j,k;
  double dist, dx,dy,dz,repel;
  MkPoint pnt;
  MkLiuParticles &par = LiuParticles;
  MkLiuBoundarys &bnd = LiuBoundarys;

  for (i=0;i<NTotal+NVirt;i++){ 
    par[i].APDVXDt = 0.e0;
    par[i].APDVYDt = 0.e0;
    par[i].APDVZDt = 0.e0;
    par[i].AVDUDt = 0.e0;
  }
  /*
  for (k=0;k<NIac;k++) {
    i = LiuPairs[k].I;
    j = LiuPairs[k].J;

    dx = par[i].X - par[j].X;
    dy = par[i].Y - par[j].Y; 
    dz = par[i].Z - par[j].Z; 
    dist = sqrt(dx*dx+dy*dy+dz*dz);
    
    if (dist > par[i].Radius + par[j].Radius) continue;
    repel = 500*(exp(-dist) - exp(-(par[i].Radius + par[j].Radius)));

    //MkDebug("Repelling %d %d with repel %f\n",i,j,repel);
    par[i].APDVXDt = par[i].APDVXDt + repel*dx/dist;
    par[j].APDVXDt = par[j].APDVXDt - repel*dx/dist;
    par[i].APDVYDt = par[i].APDVYDt + repel*dy/dist;
    par[j].APDVYDt = par[j].APDVYDt - repel*dy/dist;
    par[i].APDVZDt = par[i].APDVZDt + repel*dz/dist;
    par[j].APDVZDt = par[j].APDVZDt - repel*dz/dist;
  }
  */
  for (k=0;k<NBIac;k++) {
    i = LiuBndPairs[k].I;  //particles
    j = LiuBndPairs[k].J;  //boundaries

    dist = bnd[j].GetDistance(par[i]);
    pnt =  bnd[j].GetNearestPoint(par[i]);

    dx = par[i].X - pnt.X; 
    dy = par[i].Y - pnt.Y; 
    dz = par[i].Z - pnt.Z; 
    
    if (dist > par[i].Radius + bnd[j].Radius) continue;
    repel = 500*(exp(-dist) - exp(-(par[i].Radius + bnd[j].Radius)));

    //MkDebug("Repelling %d %d with repel %f\n",i,j,repel);
    par[i].APDVXDt = par[i].APDVXDt + repel*dx/dist;
    par[i].APDVYDt = par[i].APDVYDt + repel*dy/dist;
    par[i].APDVZDt = par[i].APDVZDt + repel*dz/dist;
  }
  // require to handle the internal energy due to repel

}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the internal forces on the right hand side  
//  of the Navier-Stokes equations, i.e. the pressure gradient and the 
//  gradient of the viscous stress tensor, used by the time integration.  
//  Moreover the entropy production due to viscous dissipation, tds/dt,  
//  and the change of internal energy per mass, de/dt, are calculated.  

//    CurrentTimeStep: Current timestep number                            [in] 
//    Dt     :   Time step                                          [in] 
//    NTotal : Number of particles                                  [in] 
//    Hsml   : Smoothing Length                                     [in] 
//    Mass   : Particle masses                                      [in] 
//    VX     : Velocities of all particles                          [in] 
//    nia//  : Number of interaction pairs                          [in] 
//    Rho    : Density                                              [in] 
//    Eta    : Dynamic viscosity                                    [in] 
//    Pair_I : List of first partner of interaction pair            [in] 
//    Pair_J : List of second partner of interaction pair           [in] 
//    DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//    ParticleType  : Type of particle (material types)                    [in] 
//    U      : Particle internal energy                             [in] 
//    X      : Particle coordinates                                 [in] 
//    t      : Particle temperature                             [in/out] 
//    C      : Particle sound speed                                [out] 
//    P      : Particle pressure                                   [out] 
//    INDVXDt  : Acceleration with respect to x, y and z             [out]  
//    TDSDt  : Production of viscous entropy                       [out] 
//    INDUDt   : Change of specific internal energy                  [out] 

void MkLiuSPH::Int_Force() 
{
  int i, j, k, d;
  double hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij;
  int maxn = NTotal + NVirt;
  MkDouble dvx(Dim), txx(maxn), tyy(maxn),tzz(maxn), txy(maxn), txz(maxn), tyz(maxn), vcc(maxn);
  MkLiuParticles &par=LiuParticles;
//   Initialization of shear tensor, velocity divergence,  
//   viscous energy, internal energy, acceleration  

   i= j= k= d=0;
   hxx= hyy= hzz= hxy= hxz= hyz= h= hvcc= he= rhoij=0;

  for (i=0;i<NTotal+NVirt;i++) {
    txx(i) = 0.e0; 
    tyy(i) = 0.e0; 
    tzz(i) = 0.e0; 
    txy(i) = 0.e0; 
    txz(i) = 0.e0; 
    tyz(i) = 0.e0; 
    vcc(i) = 0.e0; 
    par[i].TDSDt = 0.e0; 
    par[i].INDUDt = 0.e0; 
    par[i].INDVXDt = 0.e0;
    par[i].INDVYDt = 0.e0;
    par[i].INDVZDt = 0.e0;
  }
 
//   Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c 
 
  if (LiuParam.visc) {
    for (k=0;k<NIac;k++) { 
      i = LiuPairs[k].I; 
      j = LiuPairs[k].J; 
      dvx(0) = par[j].XVel - par[i].XVel;
      dvx(1) = par[j].YVel - par[i].YVel;
      dvx(2) = par[j].ZVel - par[i].ZVel;

      if (Dim==1) {
	hxx = 2.e0*dvx(1)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==2) {
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) -  dvx(1)*LiuPairs[k].DWDX(2-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - dvx(0)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==3){
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) - 
	           dvx(1)*LiuPairs[k].DWDX(2-1) - 
	           dvx(2)*LiuPairs[k].DWDX(3-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hxz = dvx(0)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - 
	           dvx(0)*LiuPairs[k].DWDX(1-1) - 
                   dvx(2)*LiuPairs[k].DWDX(3-1);
	hyz = dvx(1)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(2-1);
	hzz = 2.e0*dvx(2)*LiuPairs[k].DWDX(3-1) - 
	           dvx(0)*LiuPairs[k].DWDX(1-1) - 
                   dvx(1)*LiuPairs[k].DWDX(2-1);
      }
      hxx = 2.e0/3.e0*hxx;
      hyy = 2.e0/3.e0*hyy;
      hzz = 2.e0/3.e0*hzz;
      if (Dim==1) {
	txx(i) = txx(i) + par[j].Mass*hxx/par[j].Rho;
	txx(j) = txx(j) + par[i].Mass*hxx/par[i].Rho;
      }
      else if (Dim==2) {
	txx(i) = txx(i) + par[j].Mass*hxx/par[j].Rho;
	txx(j) = txx(j) + par[i].Mass*hxx/par[i].Rho;
	txy(i) = txy(i) + par[j].Mass*hxy/par[j].Rho;
	txy(j) = txy(j) + par[i].Mass*hxy/par[i].Rho;             
	tyy(i) = tyy(i) + par[j].Mass*hyy/par[j].Rho; 
	tyy(j) = tyy(j) + par[i].Mass*hyy/par[i].Rho;
      }
      else if (Dim==3) {
	txx(i) = txx(i) + par[j].Mass*hxx/par[j].Rho;
	txx(j) = txx(j) + par[i].Mass*hxx/par[i].Rho;    
	txy(i) = txy(i) + par[j].Mass*hxy/par[j].Rho;
	txy(j) = txy(j) + par[i].Mass*hxy/par[i].Rho;
	txz(i) = txz(i) + par[j].Mass*hxz/par[j].Rho;
	txz(j) = txz(j) + par[i].Mass*hxz/par[i].Rho;
	tyy(i) = tyy(i) + par[j].Mass*hyy/par[j].Rho;
	tyy(j) = tyy(j) + par[i].Mass*hyy/par[i].Rho;
	tyz(i) = tyz(i) + par[j].Mass*hyz/par[j].Rho;
	tyz(j) = tyz(j) + par[i].Mass*hyz/par[i].Rho;
	tzz(i) = tzz(i) + par[j].Mass*hzz/par[j].Rho;
	tzz(j) = tzz(j) + par[i].Mass*hzz/par[i].Rho;
      }
 
//   Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz: 
 
      hvcc = 0.;
      for (d=0;d<Dim;d++) {
	hvcc = hvcc + dvx(d)*LiuPairs[k].DWDX(d);
      }
      vcc(i) = vcc(i) + par[j].Mass*hvcc/par[j].Rho;
      vcc(j) = vcc(j) + par[i].Mass*hvcc/par[i].Rho;
    }
  }
 
  for (i=0;i<NTotal+NVirt;i++) {
//   Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab 
    if (LiuParam.visc) {
      if (Dim==1) {
	par[i].TDSDt = txx(i)*txx(i);
      }
      else if (Dim==2) {
	par[i].TDSDt = txx(i)*txx(i) + 2.e0*txy(i)*txy(i)+ tyy(i)*tyy(i);
      }
      else if (Dim==3) {
	par[i].TDSDt = txx(i)*txx(i)+ tyy(i)*tyy(i)+ tzz(i)*tzz(i)
                      + 2.e0*txy(i)*txy(i)+ 2.e0*txz(i)*txz(i)+ 2.e0*tyz(i)*tyz(i)  ;
      }
      par[i].TDSDt = 0.5e0*par[i].Eta/par[i].Rho*par[i].TDSDt;
    }
 
//   Pressure from equation of state 
 
    if (abs(par[i].ParticleType)==1) {
      par[i].Press = P_ideal_gas(par[i].Rho, par[i].Energy);
      par[i].SoundSpeed = C_ideal_gas(par[i].Rho, par[i].Energy);
    }
    else if (abs(par[i].ParticleType)==2) {
      par[i].Press = P_art_water(par[i].Rho);
      par[i].SoundSpeed = C_art_water(par[i].Rho);
    }

  //MkDebug("%d X:%f Y:%f Rho %f Pressure %f\n",i, par[i].X, par[i].Y, par[i].Rho, par[i].Press);
  } 
//    Calculate SPH sum for pressure force -p,a/rho 
//    and viscous force (eta Tab),b/rho 
//    and the internal energy change de/dt due to -p/rho vc,c 
 
  for (k=0;k<NIac;k++) {
    i = LiuPairs[k].I;
    j = LiuPairs[k].J;
    he = 0.e0;
       
//   For SPH algorithm 1 
 
    rhoij = 1.e0/(par[i].Rho*par[j].Rho);
    if(LiuParam.pa_sph==1){
      
      for (d=0;d<Dim ;d++) {
         
//   Pressure part 
                     
	h = -(par[i].Press + par[j].Press)*LiuPairs[k].DWDX(d);
	if (d==0) he = he + (par[j].XVel - par[i].XVel)*h;
	if (d==1) he = he + (par[j].YVel - par[i].YVel)*h;
	if (d==2) he = he + (par[j].ZVel - par[i].ZVel)*h;
 
//   Viscous force 
 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h = h + (par[i].Eta*txx(i) + par[j].Eta*txx(j))*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2) {
	      h = h + (par[i].Eta*txy(i) + par[j].Eta*txy(j))*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h = h + (par[i].Eta*txz(i) + par[j].Eta*txz(j))*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (par[i].Eta*txy(i) + par[j].Eta*txy(j))*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyy(i) + par[j].Eta*tyy(j))*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (par[i].Eta*tyz(i) + par[j].Eta*tyz(j))*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (par[i].Eta*txz(i) + par[j].Eta*txz(j))*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyz(i) + par[j].Eta*tyz(j))*LiuPairs[k].DWDX(2-1)
	          + (par[i].Eta*tzz(i) + par[j].Eta*tzz(j))*LiuPairs[k].DWDX(3-1);
	  }
	}
	h = h*rhoij;
	//MkDebug("Int_Force i:%d, j:%d, h:%f, rhoij:%f\n",i,j,h,rhoij);

	if (d==0) { 
	  par[i].INDVXDt = par[i].INDVXDt + par[j].Mass*h;
	  par[j].INDVXDt = par[j].INDVXDt - par[i].Mass*h;
	}
	if (d==1) { 
	  par[i].INDVYDt = par[i].INDVYDt + par[j].Mass*h;
	  par[j].INDVYDt = par[j].INDVYDt - par[i].Mass*h;
	}
	if (d==2) { 
	  par[i].INDVZDt = par[i].INDVZDt + par[j].Mass*h;
	  par[j].INDVZDt = par[j].INDVZDt - par[i].Mass*h;
	}
      }

      he = he*rhoij;
      par[i].INDUDt = par[i].INDUDt + par[j].Mass*he;
      par[j].INDUDt = par[j].INDUDt + par[i].Mass*he;
    }
//   For SPH algorithm 2 
           
    else if (LiuParam.pa_sph==2){
      double r2i, r2j;
      r2i = pow(par[i].Rho,2);
      r2j = pow(par[j].Rho,2);
      for (d=0;d<Dim;d++) {
	h = -(par[i].Press/r2i + par[j].Press/r2j)*LiuPairs[k].DWDX(d);

	if (d==1) he = he + (par[j].XVel - par[i].XVel)*h;
	if (d==2) he = he + (par[j].YVel - par[i].YVel)*h;
	if (d==3) he = he + (par[j].ZVel - par[i].ZVel)*h;

 //   Viscous force 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h += (par[i].Eta*txx(i)/r2i + par[j].Eta*txx(j)/r2j)*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2){
	      h += (par[i].Eta*txy(i)/r2i+par[j].Eta*txy(j)/r2j)*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h += (par[i].Eta*txz(i)/r2i+par[j].Eta*txz(j)/r2j)*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }             
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (par[i].Eta*txy(i)/r2i
		   + par[j].Eta*txy(j)/r2j)*LiuPairs[k].DWDX(1-1)
		  + (par[i].Eta*tyy(i)/r2i
		   + par[j].Eta*tyy(j)/r2j)*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (par[i].Eta*tyz(i)/r2i
		     + par[j].Eta*tyz(j)/r2j)*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (par[i].Eta*txz(i)/r2i 
		   + par[j].Eta*txz(j)/r2j)*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyz(i)/r2i 
	           + par[j].Eta*tyz(j)/r2j)*LiuPairs[k].DWDX(2-1)
	          + (par[i].Eta*tzz(i)/r2i 
	           + par[j].Eta*tzz(j)/r2j)*LiuPairs[k].DWDX(3-1);
	  }
	}
	//MkDebug("Int_Force i:%d, j:%d, h:%f, Mass[i]:%f,Mass[i]*h:%f\n",i,j,h,par[i].Mass,par[i].Mass*h);
	if(d==0) {
	  par[i].INDVXDt = par[i].INDVXDt + par[j].Mass*h;
	  par[j].INDVXDt = par[j].INDVXDt - par[i].Mass*h;
	}
	if(d==1) {
	  par[i].INDVYDt = par[i].INDVYDt + par[j].Mass*h; 
	  par[j].INDVYDt = par[j].INDVYDt - par[i].Mass*h; 
	}
	if(d==2) {
	  par[i].INDVZDt = par[i].INDVZDt + par[j].Mass*h;
	  par[j].INDVZDt = par[j].INDVZDt - par[i].Mass*h;
	}
      }
      par[i].INDUDt = par[i].INDUDt + par[j].Mass*he;
      par[j].INDUDt = par[j].INDUDt + par[i].Mass*he;
    }
  }
 
//   Change of specific internal energy de/dt = T ds/dt - p/rho vc,c: 
 
  for (i=0;i<NTotal+NVirt;i++) {
    par[i].INDUDt = par[i].TDSDt + 0.5e0*par[i].INDUDt;
    //    if(fabs(par[i].INDVXDt)>0.000001) MkDebug("i:%d, INDVXDt:%12.6f\n",i,par[i].INDVXDt);
  }
}

void MkLiuSPH::Int_Bnd_Force() 
{
  int i, j, k, d,I,J;
  double hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij;
  int maxn = NTotal + NVirt;
  int bndn = LiuBoundarys.GetSize()+NTotal;
  MkDouble dvx(Dim);
  MkDouble txx(maxn), tyy(maxn),tzz(maxn), txy(maxn), txz(maxn), tyz(maxn), vcc(maxn);
  MkDouble b_txx(bndn), b_tyy(bndn),b_tzz(bndn), b_txy(bndn), b_txz(bndn), b_tyz(bndn), b_vcc(bndn);  //0..bnd.GetSize()..bnd.GetSize()+NTotal  bnd first and then particles
  MkLiuParticles &par=LiuParticles;
  MkLiuBoundarys &bnd=LiuBoundarys;

//   Initialization of shear tensor, velocity divergence,  
//   viscous energy, internal energy, acceleration  

  i= j= k= d=0;
  hxx= hyy= hzz= hxy= hxz= hyz= h= hvcc= he= rhoij=0;

  for (i=0;i<NTotal+NVirt;i++) {
    txx(i) = 0.e0; 
    tyy(i) = 0.e0; 
    tzz(i) = 0.e0; 
    txy(i) = 0.e0; 
    txz(i) = 0.e0; 
    tyz(i) = 0.e0; 
    vcc(i) = 0.e0; 
    par[i].TDSDt = 0.e0; 
    par[i].INDUDt = 0.e0; 
    par[i].INDVXDt = 0.e0;
    par[i].INDVYDt = 0.e0;
    par[i].INDVZDt = 0.e0;
  }

  for (i=0;i<bnd.GetSize();i++) {
    b_txx(i) = 0.e0; 
    b_tyy(i) = 0.e0; 
    b_tzz(i) = 0.e0; 
    b_txy(i) = 0.e0; 
    b_txz(i) = 0.e0; 
    b_tyz(i) = 0.e0; 
    b_vcc(i) = 0.e0; 
    bnd[i].TDSDt = 0.e0; 
    bnd[i].INDUDt = 0.e0; 
    bnd[i].INDVXDt = 0.e0;
    bnd[i].INDVYDt = 0.e0;
    bnd[i].INDVZDt = 0.e0;
  }
 
//   Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c 
 
  if (LiuParam.visc) {
    for (k=0;k<NIac;k++) { 
      i = LiuPairs[k].I; 
      j = LiuPairs[k].J; 
      dvx(0) = par[j].XVel - par[i].XVel;
      dvx(1) = par[j].YVel - par[i].YVel;
      dvx(2) = par[j].ZVel - par[i].ZVel;

      if (Dim==1) {
	hxx = 2.e0*dvx(1)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==2) {
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) -  dvx(1)*LiuPairs[k].DWDX(2-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - dvx(0)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==3){
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) - 
	           dvx(1)*LiuPairs[k].DWDX(2-1) - 
	           dvx(2)*LiuPairs[k].DWDX(3-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hxz = dvx(0)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - 
	           dvx(0)*LiuPairs[k].DWDX(1-1) - 
                   dvx(2)*LiuPairs[k].DWDX(3-1);
	hyz = dvx(1)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(2-1);
	hzz = 2.e0*dvx(2)*LiuPairs[k].DWDX(3-1) - 
	           dvx(0)*LiuPairs[k].DWDX(1-1) - 
                   dvx(1)*LiuPairs[k].DWDX(2-1);
      }
      hxx = 2.e0/3.e0*hxx;
      hyy = 2.e0/3.e0*hyy;
      hzz = 2.e0/3.e0*hzz;
      if (Dim==1) {
	txx(i) = txx(i) + par[j].Mass*hxx/par[j].Rho;
	txx(j) = txx(j) + par[i].Mass*hxx/par[i].Rho;
      }
      else if (Dim==2) {
	txx(i) = txx(i) + par[j].Mass*hxx/par[j].Rho;
	txx(j) = txx(j) + par[i].Mass*hxx/par[i].Rho;
	txy(i) = txy(i) + par[j].Mass*hxy/par[j].Rho;
	txy(j) = txy(j) + par[i].Mass*hxy/par[i].Rho;             
	tyy(i) = tyy(i) + par[j].Mass*hyy/par[j].Rho; 
	tyy(j) = tyy(j) + par[i].Mass*hyy/par[i].Rho;
      }
      else if (Dim==3) {
	txx(i) = txx(i) + par[j].Mass*hxx/par[j].Rho;
	txx(j) = txx(j) + par[i].Mass*hxx/par[i].Rho;    
	txy(i) = txy(i) + par[j].Mass*hxy/par[j].Rho;
	txy(j) = txy(j) + par[i].Mass*hxy/par[i].Rho;
	txz(i) = txz(i) + par[j].Mass*hxz/par[j].Rho;
	txz(j) = txz(j) + par[i].Mass*hxz/par[i].Rho;
	tyy(i) = tyy(i) + par[j].Mass*hyy/par[j].Rho;
	tyy(j) = tyy(j) + par[i].Mass*hyy/par[i].Rho;
	tyz(i) = tyz(i) + par[j].Mass*hyz/par[j].Rho;
	tyz(j) = tyz(j) + par[i].Mass*hyz/par[i].Rho;
	tzz(i) = tzz(i) + par[j].Mass*hzz/par[j].Rho;
	tzz(j) = tzz(j) + par[i].Mass*hzz/par[i].Rho;
      }
 
//   Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz: 
 
      hvcc = 0.;
      for (d=0;d<Dim;d++) {
	hvcc = hvcc + dvx(d)*LiuPairs[k].DWDX(d);
      }
      vcc(i) = vcc(i) + par[j].Mass*hvcc/par[j].Rho;
      vcc(j) = vcc(j) + par[i].Mass*hvcc/par[i].Rho;
    }

    for (k=0;k<NBIac;k++) { 
      i = LiuBndPairs[k].I; // particle
      j = LiuBndPairs[k].J; // boundary
      dvx(0) = bnd[j].XVel - par[i].XVel;
      dvx(1) = bnd[j].YVel - par[i].YVel;
      dvx(2) = bnd[j].ZVel - par[i].ZVel;

      if (Dim==1) {
	hxx = 2.e0*dvx(1)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==2) {
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) -  dvx(1)*LiuPairs[k].DWDX(2-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - dvx(0)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==3){
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) - 
	           dvx(1)*LiuPairs[k].DWDX(2-1) - 
	           dvx(2)*LiuPairs[k].DWDX(3-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hxz = dvx(0)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - 
	           dvx(0)*LiuPairs[k].DWDX(1-1) - 
                   dvx(2)*LiuPairs[k].DWDX(3-1);
	hyz = dvx(1)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(2-1);
	hzz = 2.e0*dvx(2)*LiuPairs[k].DWDX(3-1) - 
	           dvx(0)*LiuPairs[k].DWDX(1-1) - 
                   dvx(1)*LiuPairs[k].DWDX(2-1);
      }
      hxx = 2.e0/3.e0*hxx;
      hyy = 2.e0/3.e0*hyy;
      hzz = 2.e0/3.e0*hzz;
      I = i; J = bnd.GetSize()+j; // particle I bnd J
      if (Dim==1) {
	b_txx(I) += bnd[j].Mass*hxx/bnd[j].Rho;
	b_txx(J) += par[i].Mass*hxx/par[i].Rho;
      }
      else if (Dim==2) {
	b_txx(I) += bnd[j].Mass*hxx/bnd[j].Rho;
	b_txx(J) += par[i].Mass*hxx/par[i].Rho;
	b_txy(I) += bnd[j].Mass*hxy/bnd[j].Rho;
	b_txy(J) += par[i].Mass*hxy/par[i].Rho;             
	b_tyy(I) += bnd[j].Mass*hyy/bnd[j].Rho; 
	b_tyy(J) += par[i].Mass*hyy/par[i].Rho;
      }
      else if (Dim==3) {
	b_txx(I) += bnd[j].Mass*hxx/bnd[j].Rho;
	b_txx(J) += par[i].Mass*hxx/par[i].Rho;    
	b_txy(I) += bnd[j].Mass*hxy/bnd[j].Rho;
	b_txy(J) += par[i].Mass*hxy/par[i].Rho;
	b_txz(I) += bnd[j].Mass*hxz/bnd[j].Rho;
	b_txz(J) += par[i].Mass*hxz/par[i].Rho;
	b_tyy(I) += bnd[j].Mass*hyy/bnd[j].Rho;
	b_tyy(J) += par[i].Mass*hyy/par[i].Rho;
	b_tyz(I) += bnd[j].Mass*hyz/bnd[j].Rho;
	b_tyz(J) += par[i].Mass*hyz/par[i].Rho;
	b_tzz(I) += bnd[j].Mass*hzz/bnd[j].Rho;
	b_tzz(J) += par[i].Mass*hzz/par[i].Rho;
      }
 
//   Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz: 
 
      hvcc = 0.;
      for (d=0;d<Dim;d++) {
	hvcc = hvcc + dvx(d)*LiuPairs[k].DWDX(d);
      }
      b_vcc(I) += bnd[j].Mass*hvcc/bnd[j].Rho;
      b_vcc(J) += par[i].Mass*hvcc/par[i].Rho;
    }
  }
 
  for (i=0;i<NTotal+NVirt;i++) {
//   Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab 
    if (LiuParam.visc) {
      if (Dim==1) {
	par[i].TDSDt = txx(i)*txx(i);
      }
      else if (Dim==2) {
	par[i].TDSDt = txx(i)*txx(i) + 2.e0*txy(i)*txy(i)+ tyy(i)*tyy(i);
      }
      else if (Dim==3) {
	par[i].TDSDt = txx(i)*txx(i)+ tyy(i)*tyy(i)+ tzz(i)*tzz(i)
                      + 2.e0*txy(i)*txy(i)+ 2.e0*txz(i)*txz(i)+ 2.e0*tyz(i)*tyz(i)  ;
      }
      par[i].TDSDt = 0.5e0*par[i].Eta/par[i].Rho*par[i].TDSDt;
    }
 
//   Pressure from equation of state 
 
    if (abs(par[i].ParticleType)==1) {
      par[i].Press = P_ideal_gas(par[i].Rho, par[i].Energy);
      par[i].SoundSpeed = C_ideal_gas(par[i].Rho, par[i].Energy);
    }
    else if (abs(par[i].ParticleType)==2) {
      par[i].Press = P_art_water(par[i].Rho);
      par[i].SoundSpeed = C_art_water(par[i].Rho);
    }

  //MkDebug("%d X:%f Y:%f Rho %f Pressure %f\n",i, par[i].X, par[i].Y, par[i].Rho, par[i].Press);
  } 

//    Calculate SPH sum for pressure force -p,a/rho 
//    and viscous force (eta Tab),b/rho 
//    and the internal energy change de/dt due to -p/rho vc,c 
 
  for (k=0;k<NIac;k++) {
    i = LiuPairs[k].I;
    j = LiuPairs[k].J;
    he = 0.e0;
       
//   For SPH algorithm 1 
 
    rhoij = 1.e0/(par[i].Rho*par[j].Rho);
    if(LiuParam.pa_sph==1){
      
      for (d=0;d<Dim ;d++) {
         
//   Pressure part 
                     
	h = -(par[i].Press + par[j].Press)*LiuPairs[k].DWDX(d);
	if (d==0) he = he + (par[j].XVel - par[i].XVel)*h;
	if (d==1) he = he + (par[j].YVel - par[i].YVel)*h;
	if (d==2) he = he + (par[j].ZVel - par[i].ZVel)*h;
 
//   Viscous force 
 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h = h + (par[i].Eta*txx(i) + par[j].Eta*txx(j))*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2) {
	      h = h + (par[i].Eta*txy(i) + par[j].Eta*txy(j))*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h = h + (par[i].Eta*txz(i) + par[j].Eta*txz(j))*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (par[i].Eta*txy(i) + par[j].Eta*txy(j))*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyy(i) + par[j].Eta*tyy(j))*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (par[i].Eta*tyz(i) + par[j].Eta*tyz(j))*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (par[i].Eta*txz(i) + par[j].Eta*txz(j))*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyz(i) + par[j].Eta*tyz(j))*LiuPairs[k].DWDX(2-1)
	          + (par[i].Eta*tzz(i) + par[j].Eta*tzz(j))*LiuPairs[k].DWDX(3-1);
	  }
	}
	h = h*rhoij;
	//MkDebug("Int_Force i:%d, j:%d, h:%f, rhoij:%f\n",i,j,h,rhoij);

	if (d==0) { 
	  par[i].INDVXDt = par[i].INDVXDt + par[j].Mass*h;
	  par[j].INDVXDt = par[j].INDVXDt - par[i].Mass*h;
	}
	if (d==1) { 
	  par[i].INDVYDt = par[i].INDVYDt + par[j].Mass*h;
	  par[j].INDVYDt = par[j].INDVYDt - par[i].Mass*h;
	}
	if (d==2) { 
	  par[i].INDVZDt = par[i].INDVZDt + par[j].Mass*h;
	  par[j].INDVZDt = par[j].INDVZDt - par[i].Mass*h;
	}
      }

      he = he*rhoij;
      par[i].INDUDt = par[i].INDUDt + par[j].Mass*he;
      par[j].INDUDt = par[j].INDUDt + par[i].Mass*he;
    }
//   For SPH algorithm 2 
           
    else if (LiuParam.pa_sph==2){
      double r2i, r2j;
      r2i = pow(par[i].Rho,2);
      r2j = pow(par[j].Rho,2);
      for (d=0;d<Dim;d++) {
	h = -(par[i].Press/r2i + par[j].Press/r2j)*LiuPairs[k].DWDX(d);

	if (d==1) he = he + (par[j].XVel - par[i].XVel)*h;
	if (d==2) he = he + (par[j].YVel - par[i].YVel)*h;
	if (d==3) he = he + (par[j].ZVel - par[i].ZVel)*h;

 //   Viscous force 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h += (par[i].Eta*txx(i)/r2i + par[j].Eta*txx(j)/r2j)*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2){
	      h += (par[i].Eta*txy(i)/r2i+par[j].Eta*txy(j)/r2j)*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h += (par[i].Eta*txz(i)/r2i+par[j].Eta*txz(j)/r2j)*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }             
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (par[i].Eta*txy(i)/r2i
		   + par[j].Eta*txy(j)/r2j)*LiuPairs[k].DWDX(1-1)
		  + (par[i].Eta*tyy(i)/r2i
		   + par[j].Eta*tyy(j)/r2j)*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (par[i].Eta*tyz(i)/r2i
		     + par[j].Eta*tyz(j)/r2j)*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (par[i].Eta*txz(i)/r2i 
		   + par[j].Eta*txz(j)/r2j)*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyz(i)/r2i 
	           + par[j].Eta*tyz(j)/r2j)*LiuPairs[k].DWDX(2-1)
	          + (par[i].Eta*tzz(i)/r2i 
	           + par[j].Eta*tzz(j)/r2j)*LiuPairs[k].DWDX(3-1);
	  }
	}
	//MkDebug("Int_Force i:%d, j:%d, h:%f, Mass[i]:%f,Mass[i]*h:%f\n",i,j,h,par[i].Mass,par[i].Mass*h);
	if(d==0) {
	  par[i].INDVXDt = par[i].INDVXDt + par[j].Mass*h;
	  par[j].INDVXDt = par[j].INDVXDt - par[i].Mass*h;
	}
	if(d==1) {
	  par[i].INDVYDt = par[i].INDVYDt + par[j].Mass*h; 
	  par[j].INDVYDt = par[j].INDVYDt - par[i].Mass*h; 
	}
	if(d==2) {
	  par[i].INDVZDt = par[i].INDVZDt + par[j].Mass*h;
	  par[j].INDVZDt = par[j].INDVZDt - par[i].Mass*h;
	}
      }
      par[i].INDUDt = par[i].INDUDt + par[j].Mass*he;
      par[j].INDUDt = par[j].INDUDt + par[i].Mass*he;
    }
  }

  for (i=0;i<bnd.GetSize();i++) {
//   Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab 
    if (LiuParam.visc) {
      if (Dim==1) {
	bnd[i].TDSDt = b_txx(i)*b_txx(i);
      }
      else if (Dim==2) {
	bnd[i].TDSDt = b_txx(i)*b_txx(i) + 2.e0*b_txy(i)*b_txy(i)+ b_tyy(i)*b_tyy(i);
      }
      else if (Dim==3) {
	bnd[i].TDSDt = b_txx(i)*b_txx(i)+ b_tyy(i)*b_tyy(i)+ b_tzz(i)*b_tzz(i)
                      + 2.e0*b_txy(i)*b_txy(i)+ 2.e0*b_txz(i)*b_txz(i)+ 2.e0*b_tyz(i)*b_tyz(i)  ;
      }
      bnd[i].TDSDt = 0.5e0*bnd[i].Eta/bnd[i].Rho*bnd[i].TDSDt;
    }
 
//   Pressure from equation of state 
 
    if (abs(bnd[i].BoundaryType)==1) {
      bnd[i].Press = P_ideal_gas(bnd[i].Rho, bnd[i].Energy);
      bnd[i].SoundSpeed = C_ideal_gas(bnd[i].Rho, bnd[i].Energy);
    }
    else if (abs(bnd[i].BoundaryType)==2) {
      bnd[i].Press = P_art_water(bnd[i].Rho);
      bnd[i].SoundSpeed = C_art_water(bnd[i].Rho);
    }

  //MkDebug("%d X:%f Y:%f Rho %f Pressure %f\n",i, par[i].X, par[i].Y, par[i].Rho, par[i].Press);
  } 

  for (k=0;k<NBIac;k++) {
    i = LiuBndPairs[k].I;
    j = LiuBndPairs[k].J;
    I = i;
    J = bnd.GetSize()+j;
    he = 0.e0;
       
//   For SPH algorithm 1 
 
    rhoij = 1.e0/(par[i].Rho*bnd[j].Rho);
    if(LiuParam.pa_sph==1){
      
      for (d=0;d<Dim ;d++) {
         
//   Pressure part 
                     
	h = -(par[i].Press + bnd[j].Press)*LiuPairs[k].DWDX(d);
	if (d==0) he = he + (bnd[j].XVel - par[i].XVel)*h;
	if (d==1) he = he + (bnd[j].YVel - par[i].YVel)*h;
	if (d==2) he = he + (bnd[j].ZVel - par[i].ZVel)*h;
 
//   Viscous force 
 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h = h + (par[i].Eta*b_txx(I) + bnd[j].Eta*b_txx(J))*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2) {
	      h = h + (par[i].Eta*txy(I) + bnd[j].Eta*txy(J))*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h = h + (par[i].Eta*txz(I) + bnd[j].Eta*txz(J))*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (par[i].Eta*txy(I) + bnd[j].Eta*txy(J))*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyy(I) + bnd[j].Eta*tyy(J))*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (par[i].Eta*tyz(I) + bnd[j].Eta*tyz(J))*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (par[i].Eta*txz(I) + bnd[j].Eta*txz(J))*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*tyz(I) + bnd[j].Eta*tyz(J))*LiuPairs[k].DWDX(2-1)
	          + (par[i].Eta*tzz(I) + bnd[j].Eta*tzz(J))*LiuPairs[k].DWDX(3-1);
	  }
	}
	h = h*rhoij;
	//MkDebug("Int_Force i:%d, j:%d, h:%f, rhoij:%f\n",i,j,h,rhoij);

	if (d==0) { 
	  par[i].INDVXDt = par[i].INDVXDt + bnd[j].Mass*h;
	  bnd[j].INDVXDt = bnd[j].INDVXDt - par[i].Mass*h;
	}
	if (d==1) { 
	  par[i].INDVYDt = par[i].INDVYDt + bnd[j].Mass*h;
	  bnd[j].INDVYDt = bnd[j].INDVYDt - par[i].Mass*h;
	}
	if (d==2) { 
	  par[i].INDVZDt = par[i].INDVZDt + bnd[j].Mass*h;
	  bnd[j].INDVZDt = bnd[j].INDVZDt - par[i].Mass*h;
	}
      }

      he = he*rhoij;
      par[i].INDUDt = par[i].INDUDt + bnd[j].Mass*he;
      bnd[j].INDUDt = bnd[j].INDUDt + par[i].Mass*he;
    }
//   For SPH algorithm 2 
           
    else if (LiuParam.pa_sph==2){
      double r2i, r2j;
      r2i = pow(par[i].Rho,2);
      r2j = pow(bnd[j].Rho,2);
      for (d=0;d<Dim;d++) {
	h = -(par[i].Press/r2i + bnd[j].Press/r2j)*LiuPairs[k].DWDX(d);

	if (d==1) he = he + (bnd[j].XVel - par[i].XVel)*h;
	if (d==2) he = he + (bnd[j].YVel - par[i].YVel)*h;
	if (d==3) he = he + (bnd[j].ZVel - par[i].ZVel)*h;

 //   Viscous force 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h += (par[i].Eta*b_txx(I)/r2i + bnd[j].Eta*b_txx(J)/r2j)*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2){
	      h += (par[i].Eta*b_txy(I)/r2i+bnd[j].Eta*b_txy(J)/r2j)*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h += (par[i].Eta*b_txz(I)/r2i+bnd[j].Eta*b_txz(J)/r2j)*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }             
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (par[i].Eta*b_txy(I)/r2i
		   + bnd[j].Eta*b_txy(J)/r2j)*LiuPairs[k].DWDX(1-1)
		  + (par[i].Eta*b_tyy(I)/r2i
		   + bnd[j].Eta*b_tyy(J)/r2j)*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (par[i].Eta*b_tyz(I)/r2i
		     + bnd[j].Eta*b_tyz(J)/r2j)*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (par[i].Eta*b_txz(I)/r2i 
		   + bnd[j].Eta*b_txz(J)/r2j)*LiuPairs[k].DWDX(1-1) 
	          + (par[i].Eta*b_tyz(I)/r2i 
	           + bnd[j].Eta*b_tyz(J)/r2j)*LiuPairs[k].DWDX(2-1)
	          + (par[i].Eta*b_tzz(I)/r2i 
	           + bnd[j].Eta*b_tzz(J)/r2j)*LiuPairs[k].DWDX(3-1);
	  }
	}
	//MkDebug("Int_Force i:%d, j:%d, h:%f, Mass[i]:%f,Mass[i]*h:%f\n",i,j,h,bnd[i].Mass,bnd[i].Mass*h);
	if(d==0) {
	  par[i].INDVXDt = par[i].INDVXDt + bnd[j].Mass*h;
	  bnd[j].INDVXDt = bnd[j].INDVXDt - par[i].Mass*h;
	}
	if(d==1) {
	  par[i].INDVYDt = par[i].INDVYDt + bnd[j].Mass*h; 
	  bnd[j].INDVYDt = bnd[j].INDVYDt - par[i].Mass*h; 
	}
	if(d==2) {
	  par[i].INDVZDt = par[i].INDVZDt + bnd[j].Mass*h;
	  bnd[j].INDVZDt = bnd[j].INDVZDt - par[i].Mass*h;
	}
      }
      par[i].INDUDt = par[i].INDUDt + bnd[j].Mass*he;
      bnd[j].INDUDt = bnd[j].INDUDt + par[i].Mass*he;
    }
  }
 
//   Change of specific internal energy de/dt = T ds/dt - p/rho vc,c: 
 
  for (i=0;i<NTotal+NVirt;i++) {
    par[i].INDUDt = par[i].TDSDt + 0.5e0*par[i].INDUDt;
    //    if(fabs(par[i].INDVXDt)>0.000001) MkDebug("i:%d, INDVXDt:%12.6f\n",i,par[i].INDVXDt);
  }
}

//-------------------------------------------------------------------------- 
//  Subroutine to calculate the external forces, e.g. gravitational forces.       
//  The forces from the interactions with boundary virtual particles  
//  are also calculated here as external forces. 
 
//  here as the external force.  
//  NTotal  : Number of particles                                 [in] 
//  Mass    : Particle masses                                     [in] 
//  X       : Coordinates of all particles                        [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  ParticleType   : type of particles                                   [in] 
//  Hsml   : Smoothing Length                                     [in] 
//  EXDVXDt   : Acceleration with respect to x, y and z            [out]  

void MkLiuSPH::Ext_Force() 
{ 

  int i, j, k, d ;
  double rr, f, rr0, dd, p1, p2;
  MkDouble dx(Dim);
  MkLiuParticles &par=LiuParticles;

   i= j= k= d =0;
   rr= f= rr0= dd= p1= p2=0;
   //MkDebug("ext_force() ");
  for (i = 0;i< NTotal+NVirt;i++) { 
    par[i].EXDVXDt = 0.;
    par[i].EXDVYDt = 0.;
    par[i].EXDVZDt = 0.;
  }
         
//  Consider self-gravity or not ? 
 
  if (LiuParam.self_gravity) { 
    for (i = 0; i<NTotal+NVirt;i++) { 
      if (Dim == 1) par[i].EXDVXDt = -9.8;
      if (Dim == 2) par[i].EXDVYDt = -9.8;
      if (Dim == 3) par[i].EXDVYDt = -0.98;//par[i].EXDVZDt = -9.8;
    }
  }  
 
//  Boundary particle force and penalty anti-penetration force.  
  rr0 = 1.25e-5; 
  dd = 1.e-2; 
  p1 = 12; 
  p2 = 4; 
       
  for ( k=0;k<NIac;k++){ 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J;   
    if(par[i].ParticleType>0&&par[j].ParticleType<0) {   
      rr = 0.;       

      dx(0) =  par[i].X - par[j].X; 
      rr = rr + dx(0)*dx(0); 
      dx(1) =  par[i].Y - par[j].Y; 
      rr = rr + dx(1)*dx(1); 
      dx(2) =  par[i].Z - par[j].Z; 
      rr = rr + dx(2)*dx(2); 

      rr = sqrt(rr); 
      if(rr<rr0) { 
	f = (pow(rr0/rr,p1)-pow(rr0/rr,p2))/rr*rr ;

	par[i].EXDVXDt = par[i].EXDVXDt + dd*dx(0)*f;
	par[i].EXDVYDt = par[i].EXDVYDt + dd*dx(1)*f;
	par[i].EXDVZDt = par[i].EXDVZDt + dd*dx(2)*f;
      } 
    }         
  }   
  //MkDebug("~ext_force()\n ");
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the average velocity to correct velocity 
//  for preventing.penetration (monaghan, 1992) 
 
//  NTotal : Number of particles                                  [in] 
//  Mass   : Particle masses                                      [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  W      : Kernel for all interaction pairs                     [in] 
//  VX     : Velocity of each particle                            [in] 
//  Rho    : Density of each particle                             [in] 
//  AveVel     : Average velocityof each particle                    [out] 

void MkLiuSPH::Av_Vel()
{ 
  int i,j,k,d;        
  double  vcc, epsilon; 
  MkDouble dvx(Dim);
  MkLiuParticles &par = LiuParticles;

   i=j=k=d=0;        
   vcc= epsilon=0; 
       
//  epsilon --- a small constants chosen by experience, may lead to instability. 
//  for example, for the 1 dimensional shock tube problem, the E = 0.3 
 
  epsilon = 0.3;
       
  for (i = 0;i< NTotal;i++ ) { 
    par[i].XAVel = 0.0;
    par[i].YAVel = 0.0;
    par[i].ZAVel = 0.0;
  }
      
  for (k=0;k<NIac;k++) {        
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J;

    dvx(0) = par[i].XVel - par[i].XVel;
    par[i].XAVel = par[i].XAVel - 2*par[j].Mass*dvx(0)/(par[i].Rho+par[j].Rho)*LiuPairs[k].W;
    par[j].XAVel = par[j].XAVel + 2*par[i].Mass*dvx(0)/(par[i].Rho+par[j].Rho)*LiuPairs[k].W;
    dvx(1) = par[i].YVel - par[i].YVel;
    par[i].YAVel = par[i].YAVel - 2*par[j].Mass*dvx(1)/(par[i].Rho+par[j].Rho)*LiuPairs[k].W;
    par[j].YAVel = par[j].YAVel + 2*par[i].Mass*dvx(1)/(par[i].Rho+par[j].Rho)*LiuPairs[k].W;
    dvx(2) = par[i].ZVel - par[i].ZVel;
    par[i].ZAVel = par[i].ZAVel - 2*par[j].Mass*dvx(2)/(par[i].Rho+par[j].Rho)*LiuPairs[k].W;
    par[j].ZAVel = par[j].ZAVel + 2*par[i].Mass*dvx(2)/(par[i].Rho+par[j].Rho)*LiuPairs[k].W;

  }         
  for (i = 0;i< NTotal;i++) {
    par[i].XAVel = epsilon * par[i].XAVel;
    par[i].YAVel = epsilon * par[i].YAVel;
    par[i].ZAVel = epsilon * par[i].ZAVel;
  } 
}


//----------------------------------------------------------------------- 
//  Subroutine to evolve smoothing length 
 
//  Dt     : time step                                            [in] 
//  NTotal : Number of particles                                  [in] 
//  Mass   : Particle masses                                      [in] 
//  VX     : Velocities of all particles                          [in] 
//  Rho    : Density                                              [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//  Hsml   : Smoothing Length                                 [in/out] 

void MkLiuSPH::H_Upgrade() 
{    
  int i,j,k,d; 
  double fac, hvcc;
  int maxn=NTotal+NVirt;
  MkDouble dvx(Dim), vcc(maxn), dhsml(maxn);
  MkLiuParticles &par = LiuParticles;

   i=j=k=d=0; 
   fac= hvcc=0;
 
  if (LiuParam.sle==0 ) {      
 //---  Keep smoothing length unchanged.  
     return; 
  }      
  else if (LiuParam.sle==2) {
       
//---  dh/dt = (-1/Dim)*(h/rho)*(drho/dt). 
 
    for (i=0;i<NTotal ;i++) {
      vcc(i) = 0.e0 ;
    }
       
    for (k=0;k<NIac;k++) {
      i = LiuPairs[k].I; 
      j = LiuPairs[k].J ;

      dvx(0) = par[j].XVel - par[i].XVel;
      dvx(1) = par[j].YVel - par[i].YVel;
      dvx(2) = par[j].ZVel - par[i].ZVel;

      hvcc = dvx(0)*LiuPairs[k].DWDX(1-1);
      for (d=1;d<Dim;d++) {
	hvcc = hvcc + dvx(d)*LiuPairs[k].DWDX(d-1);
      }
      vcc(i) = vcc(i) + par[j].Mass*hvcc/par[j].Rho; 
      vcc(j) = vcc(j) + par[i].Mass*hvcc/par[i].Rho;
    }
         
    for (i = 0;i< NTotal;i++) {
      dhsml(i) = (par[i].SmoothLen/Dim)*vcc(i);
      par[i].SmoothLen = par[i].SmoothLen + Dt*dhsml(i);
      if (par[i].SmoothLen<=0) par[i].SmoothLen = par[i].SmoothLen - Dt*dhsml(i);  
    }
  }
     
  else if(LiuParam.sle==1) {
    fac = 2.0;
    for (i = 0;i< NTotal;i++) {
      par[i].SmoothLen = fac * (par[i].Mass/pow(par[i].Rho,(1./Dim)));
    }
        
  }
}

//---------------------------------------------------------------------- 
//   Subroutine for loading or generating initial particle information 
 
//   X-- coordinates of particles                                 [out] 
//   VX-- velocities of particles                                 [out] 
//   Mass-- mass of particles                                     [out] 
//   Rho-- dnesities of particles                                 [out] 
//   p-- pressure  of particles                                   [out] 
//   U-- internal energy of particles                             [out] 
//   ParticleType-- types of particles                                   [out] 
//   Hsml-- smoothing lengths of particles                        [out] 
//   ntotal-- total particle number                               [out] 
 
void MkLiuSPH::Input() 
{
  if (LiuParam.shocktube) Shock_Tube();
  else if (LiuParam.shearcavity) Shear_Cavity();
  else if (LiuParam.lowdensity) LowDensity();
}        

//----------------------------------------------------------------------            
//   Subroutine for saving particle information to external disk file 
 
//   X-- coordinates of particles                                  [in] 
//   VX-- velocities of particles                                  [in] 
//   Mass-- mass of particles                                      [in] 
//   Rho-- dnesities of particles                                  [in] 
//   P-- pressure  of particles                                    [in] 
//   U-- internal energy of particles                              [in] 
//   C-- sound velocity of particles                               [in] 
//   ParticleType-- types of particles                                    [in] 
//   Hsml-- smoothing lengths of particles                         [in] 
//   NTotal-- total particle number                                [in] 
 
void MkLiuSPH::Output()  
{

}


//----------------------------------------------------------------------      
//   This subroutine is used to generate initial data for the  
//   1 d noh shock tube problem 
//   X-- coordinates of particles                                 [out] 
//   VX-- velocities of particles                                 [out] 
//   Mass-- mass of particles                                     [out] 
//   Rho-- dnesities of particles                                 [out] 
//   P-- pressure  of particles                                   [out] 
//   U-- internal energy of particles                             [out] 
//   ParticleType-- types of particles                                   [out] 
//        =1   ideal gas 
//   Hsml-- smoothing lengths of particles                        [out] 
//   ntotal-- total particle number                               [out] 
        
void MkLiuSPH::Shock_Tube() 
{
  int i, d; 
  double space_x;
  MkLiuParticles &par = LiuParticles;
 
  NTotal=400; 
  space_x=0.6/80.;       

  //MkDebug("  shock_tube() : \n");

  for (i=0;i<NTotal;i++) {
    par[i].Mass=0.75/400.; 
    par[i].SmoothLen=0.015; 
    par[i].ParticleType=1;

    par[i].X = par[i].Y = par[i].Z = 0;
    par[i].XVel = par[i].YVel = par[i].ZVel = 0;
  }

  for (i=0;i<320;i++) {
    par[i].X=-0.6+space_x/4.*(i-1); 
  }
       
  for (i=320;i<NTotal;i++) {
    par[i].X=0.+space_x*(i-320);
  }         

  for (i=0;i<NTotal;i++) {
    if (par[i].X<=1.e-8)  {
	par[i].Energy=2.5; 
	par[i].Rho=1. ;
	par[i].Press=1. ;
      }
    if (par[i].X>1.e-8){
      par[i].Energy=1.795; 
      par[i].Rho=0.25 ;
      par[i].Press=0.1795 ;
    }   
  }
  //MkDebug("  ~shock_tube(): \n ");
}    

 
//----------------------------------------------------------------------      
//   This subroutine is used to generate initial data for the  
//   2 d shear driven cavity probem with Re = 1 
//   X-- coordinates of particles                                 [out] 
//   VX-- velocities of particles                                 [out] 
//   Mass-- mass of particles                                     [out] 
//   Rho-- dnesities of particles                                 [out] 
//   P-- pressure  of particles                                   [out] 
//   U-- internal energy of particles                             [out] 
//   ParticleType-- types of particles                                   [out] 
//        =2   water 
//   h-- smoothing lengths of particles                           [out] 
//   NTotal-- total particle number                               [out] 
void MkLiuSPH::Shear_Cavity() 
{
  int i, j, d, m, n, mp, np, k; 
  double xl, yl, dx, dy, rho_ref=1000; 
  MkLiuParticles &par=LiuParticles;
 
//   Giving data.Mass and smoothing length as well as other data. 
 
//MkDebug("shear_cavity()\n");
  m = 21; 
  n = 21; 
  mp = m-1; 
  np = n-1; 
  NTotal = mp * np; 
  xl = 1.0; 
  yl = 1.0; 
  dx = xl/mp; 
  dy = yl/np; 
 
  for (i = 0;i < mp;i++) {
    for (j = 0;j < np;j++) {
      k = j + i*np; 
      par[k].X = (i)*dx + dx/2.0;
      par[k].Y = (j)*dy + dy/2.0;
    }
  }
 
  for (i = 0;i < mp*np ;i++) {
    par[i].XVel = 0.; 
    par[i].YVel = 0.;       
    par[i].Rho  = rho_ref;    
    par[i].Mass = dx*dy*par[i].Rho;   
    par[i].Eta = 100000000000;
    par[i].Press= 0.;    
    par[i].Energy=357.1; 
    par[i].ParticleType = 2; 
    par[i].SmoothLen = dx;  
    par[i].Radius = dx/2.0;
  }
//MkDebug("~shear_cavity()\n");
}

void MkLiuSPH::LowDensity()
{
  int i, j, d, m, n, mp, np, k; 
  double xl, yl, dx, dy, rho_ref=1000; 
 
//   Giving data.Mass and smoothing length as well as other data. 
 
//MkDebug("low_density()\n");
  m = 20; 
  n = 20; 
  mp = m-1; 
  np = n-1; 
  NTotal = mp * np; 
  xl = 1.e-0; 
  yl = 1.e-0; 
  dx = xl/mp; 
  dy = yl/np; 
 
  for (i = 0;i < mp;i++) {
    for (j = 0;j < np;j++) {
      k = j + i*np; 
      LiuParticles[k].X = (i)*dx/1.1 + dx/2.0; //rand()%10000/15000.0+0.25;
      LiuParticles[k].Y = (j)*dy/1.1 + dy/2.0; //rand()%10000/15000.0+0.25;
    }
  }
 
  for (i = 0;i < mp*np ;i++) {
    LiuParticles[i].XVel = 0.; 
    LiuParticles[i].YVel = 0.;       
    LiuParticles[i].Rho  = rho_ref;
    LiuParticles[i].Mass = dx*dy*LiuParticles[i].Rho;   
    LiuParticles[i].Press= 0.;    
    LiuParticles[i].Energy=357.1; 
    LiuParticles[i].ParticleType = 2; 
    LiuParticles[i].SmoothLen = dx; 
    LiuParticles[i].Radius = dx/1.9; 
  }

  LiuGrids.SetSmoothLen(dx);
  LiuGrids.Setup(LiuParticles);

//MkDebug("~low_density()\n");
}

void MkLiuSPH::Virt_LD()
{
  int i, j, d, im, mp,np,nvirt; 
  double xl, dx, v_inf, rho_virt = 300; 
  FILE *fp1,*fp2, *fp3;
  MkLiuParticles &par = LiuParticles;

   i= j= d= im= mp=0; 
   xl= dx= v_inf=0; 
 
//MkDebug("    virt_ld()\n");

  if (LiuParam.vp_input){
    //    LoadVPFromFile();
  }
  else  {
    nvirt = 0; 
    mp = 5;
    np = 5;
    xl = 1.0; 
    dx = xl / mp; 
    v_inf = 1.0; 

//   Monaghan type virtual particle on the Upper side 
  //MkDebug("    virt_ld()VP1 nvirt is %d ",nvirt); 
    for (i = 0;i< 2*mp+1;i++) {
      par[NTotal + nvirt].X = i*dx/2;
      par[NTotal + nvirt].Y = xl;   
      par[NTotal + nvirt].XVel = 0;//v_inf; 
      par[NTotal + nvirt].YVel = 0.; 
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_ld()VP2 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Lower side 
 
    for (i = 0;i< 2*mp+1;i++) {
      par[NTotal + nvirt].X = i*dx/2;  
      par[NTotal + nvirt].Y = 0.;   
      par[NTotal + nvirt].XVel = 0.; 
      par[NTotal + nvirt].YVel = 0.; 
      nvirt = nvirt + 1; 
    }
  //MkDebug("    virt_ld()VP3 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Left side 
 
    for (i = 0;i< 2*np-1;i++) {
      par[NTotal + nvirt].X = 0.;
      par[NTotal + nvirt].Y = (i+1)*dx/2;
      par[NTotal + nvirt].XVel = 0.; 
      par[NTotal + nvirt].YVel = 0.;
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_ld()VP4 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Right side 
 
    for (i = 0;i< 2*np-1;i++) {
      par[NTotal + nvirt].X = xl;
      par[NTotal + nvirt].Y = (i+1)*dx/2;
      par[NTotal + nvirt].XVel = 0.; 
      par[NTotal + nvirt].YVel = 0.; 
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_ld()VP5 nvirt is %d ",nvirt); 
    for (i = 0;i< nvirt;i++) {
      par[NTotal + i].Rho = rho_virt;
      par[NTotal + i].Mass = par[NTotal + i].Rho  * dx * dx; 
      par[NTotal + i].Press = 0.; 
      par[NTotal + i].Energy = 357.1; 
      par[NTotal + i].ParticleType = -2; 
      par[NTotal + i].SmoothLen = dx/2; 
    }
  } 
  if ((CurrentTimeStep%LiuParam.save_step)==0) {

  }
 //MkDebug("    virt_ld()VP6 nvirt is %d ",nvirt); 
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    if (LiuParam.int_stat) {
      printf(" >> Statistics: Virtual boundary particles:\n"); 
      printf("          Number of virtual particles:%d\n",NVirt);
    }
  }

//MkDebug("    ~virt_ld()\n");

}

void MkLiuSPH::Shake()
{
  for (int i = 0;i < NTotal ;i++) {
    LiuParticles[i].X += 5000*1e-6 - (rand()%10000)*1e-6; 
    LiuParticles[i].Y += 5000*1e-6 - (rand()%10000)*1e-6; 
  }
}

void MkLiuSPH::Droplet(double rad, double spacing, double rho)
{
  
}

double MkLiuSPH::Calc_Bnd_Norm() // return contribution from 1 particle without self density
{
  int i,j,k, np = 2;
  double hsml = LiuParticles[0].GetSmoothLen();
  double mass = LiuParticles[0].GetMass();
  double dx,dy,dz;
  MkLiuParticles par;
  MkLiuBoundarys &bnd=LiuBoundarys;
  MkLiuPairs _pair;

//MkDebug("\n      calc_bnd_norm()\n");

  par.Initialize(np);
  _pair.Initialize(np-1);

  LiuKernel.SetSmoothLen(hsml);

//MkDebug("CBN 01 ");
  par[0].X = par[0].Y = par[0].Z = 0;
  par[0].SmoothLen = hsml;
  par[0].Mass = 0; // no self density

//MkDebug("CBN 02 ");

  for (i=1;i<np;i++) {
    par[i].X = hsml*cos(i*M_PI/3);
    par[i].Y = hsml*sin(i*M_PI/3);
    par[i].Z = 0;
    par[i].SmoothLen = hsml;
    par[i].Mass = mass;
  }

//MkDebug("CBN 03 ");

  for (i=1;i<np;i++) {
    _pair[i-1].I = 0;
    _pair[i-1].J = i;

    dx = par[0].X - par[i].X; 
    dy = par[0].Y - par[i].Y; 
    dz = par[0].Z - par[i].Z; 

    _pair[i-1].Dist = sqrt(dx*dx+dy*dy+dz*dz);
    _pair[i-1].dX = dx;
    _pair[i-1].dY = dy;
    _pair[i-1].dZ = dz;

    _pair[i-1].W = LiuKernel.W(_pair[i-1].Dist);
  }

//MkDebug("CBN 04\n ");

  par[0].Rho = par[0].Mass*LiuKernel.W(0); 

//MkDebug("mass %f, hsml %f \n",mass, hsml);
//MkDebug("par[%d].Rho = %f\n",0,par[0].Rho);
  for (k=0;k<np-1;k++) { 
    i = _pair[k].I; 
    j = _pair[k].J; 
    par[i].Rho  += par[j].Mass*_pair[k].W; 
  //MkDebug("par[%d].Rho = %f, mass %f W %f\n",0,par[0].Rho,par[j].Mass,_pair[k].W );
  }

//MkDebug("\n      ~calc_bnd_norm()\n");
  return par[0].Rho;
 
}

void MkLiuSPH::Norm_Bnd_Density()
{
  static bool norm_bnd_calc=false;

  int i, j, k;
  double selfdens, r, dist, ang, shsml;
  int maxn = NTotal + NVirt;
  MkLiuParticles &par = LiuParticles;
  MkLiuBoundarys &bnd=LiuBoundarys;

  if(!norm_bnd_calc) {
    Rho_Norm = Calc_Norm();
    Rho_Bnd_Norm = Calc_Bnd_Norm();
    norm_bnd_calc = true;
    MkDebug("Rho_Norm %f Rho_Bnd_Norm %f\n",Rho_Norm, Rho_Bnd_Norm);
    //usleep(999999);usleep(999999);usleep(999999);
  }

   i= j= k= 0;
   selfdens= r=0;
 
//MkDebug("\n      norm_bnd_density()\n");

 //MkDebug("      NBD1 ");         

//MkDebug("NBD2 ");          
//     Self density of each particle: Wii (Kernel for distance 0) 
//     and take contribution of particle itself: 
 
  r=0.;
       
//  Secondly calculate the rho integration over the space 
//MkDebug("NBD4 ");           
  for (i=0;i<NTotal+NVirt;i++) { 
    LiuKernel.SetSmoothLen(par[i].GetSmoothLen());
    par[i].Rho = (Rho_Ref/Rho_Norm)*par[i].Mass*LiuKernel.W(r); 
  }
//MkDebug("NBD5 \n");           
//  Calculate SPH sum for rho: 
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    par[i].Rho += (Rho_Ref/Rho_Norm)*par[j].Mass*LiuPairs[k].W; 
    par[j].Rho += (Rho_Ref/Rho_Norm)*par[i].Mass*LiuPairs[k].W; 
  }

  for (k=0;k<NBIac;k++) {
    i = LiuBndPairs[k].I; // par
    j = LiuBndPairs[k].J; // bnd
    dist = LiuBndPairs[k].Dist;
    shsml = par[i].GetSmoothLen()+bnd[j].GetSmoothLen();
    assert(shsml > 0.0001);
    ang = min(acos(dist/shsml),2.1419/2);
    //MkDebug("before: i %d rho %f ang %f dist %f shsml %f \n",i, par[i].Rho,ang, dist, shsml);
    //par[i].Rho += (Rho_Ref/Rho_Norm)*(3*ang/M_PI)*Rho_Bnd_Norm;
    //MkDebug("after:  i %d rho %f ang %f dist %f shsml %f \n\n\n\n",i, par[i].Rho,ang, dist, shsml);
  }
  //usleep(999999);usleep(999999);usleep(999999);
  /*
  for (k=0 ;k<(int) min(10,NIac);k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
  //MkDebug("k-%d, Roh(%d):%f,Roh(%d):%f,W:%f \n",k,i,par[i].Rho,j,par[j].Rho,LiuPairs[k].W);
  }*/
//MkDebug("NBD6 \n");           

//  Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w) 
      
//MkDebug("\n      ~norm_bnd_density()\n");       


}

double MkLiuSPH::Calc_Norm()
{
  int i,j,k, np = 7;
  double hsml = LiuParticles[0].GetSmoothLen();
  double mass = LiuParticles[0].GetMass();
  double dx,dy,dz;
  MkLiuParticles par;
  MkLiuPairs _pair;

//MkDebug("\n      calc_norm()\n");

  par.Initialize(np);
  _pair.Initialize(np-1);

  LiuKernel.SetSmoothLen(hsml);

//MkDebug("CN 01 ");
  par[0].X = par[0].Y = par[0].Z = 0;
  par[0].SmoothLen = hsml;
  par[0].Mass = mass;

  MkDebug("hsml %f, mass %f\n",hsml, mass);
  //usleep(999999);usleep(999999);usleep(999999);
//MkDebug("CN 02 ");

  for (i=1;i<np;i++) {
    par[i].X = hsml*cos(i*M_PI/3);
    par[i].Y = hsml*sin(i*M_PI/3);
    par[i].Z = 0;
    par[i].SmoothLen = hsml;
    par[i].Mass = mass;
  }

//MkDebug("CN 03 ");

  for (i=1;i<np;i++) {
    _pair[i-1].I = 0;
    _pair[i-1].J = i;

    dx = par[0].X - par[i].X; 
    dy = par[0].Y - par[i].Y; 
    dz = par[0].Z - par[i].Z; 

    _pair[i-1].Dist = sqrt(dx*dx+dy*dy+dz*dz);
    _pair[i-1].dX = dx;
    _pair[i-1].dY = dy;
    _pair[i-1].dZ = dz;

    _pair[i-1].W = LiuKernel.W(_pair[i-1].Dist);
    _pair[i-1].dWdX = LiuKernel.dWdX(_pair[i-1].Dist,dx,dy,dz);
    _pair[i-1].dWdY = LiuKernel.dWdY(_pair[i-1].Dist,dx,dy,dz);
    _pair[i-1].dWdZ = LiuKernel.dWdZ(_pair[i-1].Dist,dx,dy,dz);
  }

//MkDebug("CN 04 ");

  par[0].Rho = par[0].Mass*LiuKernel.W(0); 

  for (k=0;k<np-1;k++) { 
    i = _pair[k].I; 
    j = _pair[k].J; 
    par[i].Rho  += par[j].Mass*_pair[k].W; 
  }

//MkDebug("\n      ~calc_norm()\n");

  return par[0].Rho;
 
}

void MkLiuSPH::Norm_Density()
{
  static bool norm_calc=false;

  int i, j, k, d;
  double selfdens, r;
  int maxn = NTotal + NVirt;
  MkLiuParticles &par = LiuParticles;

  if(!norm_calc) {Rho_Norm = Calc_Norm();norm_calc = true;}

   i= j= k= d=0;
   selfdens= r=0;
 
//MkDebug("\n      norm_density()\n");

 //MkDebug("      ND1 ");         

//MkDebug("ND2 ");          
//     Self density of each particle: Wii (Kernel for distance 0) 
//     and take contribution of particle itself: 
 
  r=0.;
       
//  Secondly calculate the rho integration over the space 
//MkDebug("ND4 ");           
  for (i=0;i<NTotal+NVirt;i++) { 
    LiuKernel.SetSmoothLen(par[i].GetSmoothLen());
    par[i].Rho = (Rho_Ref/Rho_Norm)*par[i].Mass*LiuKernel.W(r); 
  }
//MkDebug("ND5 \n");           
//  Calculate SPH sum for rho: 
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    par[i].Rho += (Rho_Ref/Rho_Norm)*par[j].Mass*LiuPairs[k].W; 
    par[j].Rho += (Rho_Ref/Rho_Norm)*par[i].Mass*LiuPairs[k].W; 
  }
  /*
  for (k=0 ;k<(int) min(10,NIac);k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
  //MkDebug("k-%d, Roh(%d):%f,Roh(%d):%f,W:%f \n",k,i,par[i].Rho,j,par[j].Rho,LiuPairs[k].W);
  }*/
//MkDebug("ND6 \n");           

//  Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w) 
      
//MkDebug("\n      ~norm_density()\n");       


}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the density with SPH summation algorithm. 
//  See Equ.(4.35) 
 
//    NTotal : Number of particles                                  [in] 
//    Hsml   : Smoothing Length                                     [in] 
//    Mass   : Particle masses                                      [in] 
//    NIac   : Number of interaction pairs                          [in] 
//    Pair_I : List of first partner of interaction pair            [in] 
//    Pair_J : List of second partner of interaction pair           [in] 
//    W      : Kernel for all interaction pairs                     [in] 
//    ParticleType   : type of particles                                   [in] 
//    X       : Coordinates of all particles                        [in] 
//    Rho    : Density                                             [out] 

void MkLiuSPH::Sum_Density() 
{
  int i, j, k, d;
  double selfdens, r;
  int maxn = NTotal + NVirt;
  MkDouble hv(Dim), wi(maxn);
  MkLiuParticles &par = LiuParticles;

   i= j= k= d=0;
   selfdens= r=0;
 
//MkDebug("\n      sum_density()\n");
//    wi(data.MaxN)---integration of the kernel itself 
 //MkDebug("      SD1 ");         
  for (d=0;d<Dim;d++) { 
    hv(d) = 0.e0; 
  }
//MkDebug("SD2 ");          
//     Self density of each particle: Wii (Kernel for distance 0) 
//     and take contribution of particle itself: 
 
  r=0.;
       
//     Firstly calculate the integration of the kernel over the space 
 
  for (i=0;i<NTotal+NVirt ;i++) {
    LiuKernel.SetSmoothLen(par[i].GetSmoothLen());
    selfdens = LiuKernel.W(r);
    wi(i)=selfdens*par[i].Mass/par[i].Rho; 
  }
//MkDebug("SD3 ");           
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    wi(i) = wi(i) + par[j].Mass/par[j].Rho*LiuPairs[k].W; 
    wi(j) = wi(j) + par[i].Mass/par[i].Rho*LiuPairs[k].W; 
  }
 
//  Secondly calculate the rho integration over the space 
//MkDebug("SD4 ");           
  for (i=0;i<NTotal+NVirt;i++) { 
    LiuKernel.SetSmoothLen(par[i].GetSmoothLen());
    selfdens = LiuKernel.W(r);
    par[i].Rho = selfdens*par[i].Mass; 
  }
//MkDebug("SD5 \n");           
//  Calculate SPH sum for rho: 
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    par[i].Rho = par[i].Rho + par[j].Mass*LiuPairs[k].W; 
    par[j].Rho = par[j].Rho + par[i].Mass*LiuPairs[k].W; 
  }
  for (k=0 ;k<(int) min(10,NIac);k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
  //MkDebug("k-%d, Roh(%d):%f,Roh(%d):%f,W:%f \n",k,i,par[i].Rho,j,par[j].Rho,LiuPairs[k].W);
  }
//MkDebug("SD6 \n");           

//  Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w) 
      
  if (LiuParam.nor_density) {  

    for (i=0;i< NTotal+NVirt;i++) { 
      par[i].Rho=par[i].Rho/wi(i);
    }

    for (k=NIac-10;k<NIac;k++) { 
      i = LiuPairs[k].I; 
      j = LiuPairs[k].J; 
    //MkDebug("k-%d, Roh(%d):%f,Roh(%d):%f,wi:%f,wj:%f \n",k,i,par[i].Rho,j,par[j].Rho,wi(i), wi(j));
  }


  }  
//MkDebug("\n      ~sum_density()\n");       
}
       
//---------------------------------------------------------------------- 
//  Subroutine to calculate the density with SPH continuity approach. 
//  See Equ.(4.34) 
 
//  NTotal : Number of particles                                  [in] 
//  Mass   : Particle masses                                      [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  DWDX   : derivation of Kernel for all interaction pairs       [in] 
//  VX     : Velocities of all particles                          [in] 
//  ParticleType   : type of particles                                   [in] 
//  X      : Coordinates of all particles                         [in] 
//  Rho    : Density                                              [in] 
//  DRhoDt : Density change rate of each particle                [out]    

void MkLiuSPH::Con_Density()
{ 
       
  int i,j,k,d;
  double vcc;
  MkDouble dvx(Dim);
  MkLiuParticles &par = LiuParticles;

   i=j=k=d=0;
   vcc=0;
       
  for (i = 0;i< NTotal+NVirt;i++) { 
    par[i].DRhoDt = 0.; 
  }
      
  for (k=0;k<NIac;k++) {       
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 

    dvx(0) = par[i].XVel - par[j].XVel;
    dvx(1) = par[i].YVel - par[j].YVel;
    dvx(2) = par[i].ZVel - par[j].ZVel;

    vcc = dvx(0)*LiuPairs[k].DWDX(1-1);
    for (d=1;d<Dim;d++) { 
      vcc = vcc + dvx(d)*LiuPairs[k].DWDX(d-1); 
    }
    par[i].DRhoDt = par[i].DRhoDt + par[j].Mass*vcc; 
    par[j].DRhoDt = par[j].DRhoDt + par[i].Mass*vcc;
  }
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the smoothing funciton for each particle and 
//  the interaction parameters used by the SPH algorithm. Interaction  
//  pairs are determined by directly comparing the particle distance  
//  with the corresponding smoothing length. 
//  See p.148 in Chapter 4 
 
//    CurrentTimeStep : Current time step                                 [in] 
//    NTotal    : Number of particles                               [in] 
//    Hsml      : Smoothing Length                                  [in] 
//    X         : Coordinates of all particles                      [in] 
//    NIac      : Number of interaction pairs                      [out] 
//    Pair_I    : List of first partner of interaction pair        [out] 
//    Pair_J    : List of second partner of interaction pair       [out] 
//    W         : Kernel for all interaction pairs                 [out] 
//    DWDX      : Derivative of kernel with respect to x, y and z  [out] 
//    CountIac  : Number of neighboring particles                  [out] 

int MkLiuSPH::Pair_Count() 
{
  int i, j, niac=0, scale_k;  
  double  dist, mhsml, dx,dy,dz;      
  MkLiuParticles &par = LiuParticles;

  if (LiuParam.skf==1) scale_k = 2;  
  else scale_k = 3; 

  for (i=0;i<NTotal+NVirt-1;i++) {      
    for (j = i+1;j< NTotal+NVirt;j++){ 
      float Xi,Xj,Yi,Yj;
      Xi = par[i].X;
      Xj = par[j].X;
      Yi = par[i].Y;
      Yj = par[j].Y;

      dx = par[i].X - par[j].X; 
      dy = par[i].Y - par[j].Y; 
      dz = par[i].Z - par[j].Z; 
      dist = sqrt(dx*dx+dy*dy+dz*dz);
      mhsml = (par[i].SmoothLen+par[j].SmoothLen)/2.; 
      if (dist < scale_k*mhsml) {
	//MkDebug("%d (%f,%f), %d(%f,%f) \n",i,Xi,Yi,j,Xj,Yj);
	niac = niac + 1;  
      }
    }
  }
//MkDebug("Number of pair count %d\n",niac);
  return niac;
}

void MkLiuSPH::Direct_Find() 
{ 
  int i, j, d,  sumiac, maxiac, miniac, niac, noiac, maxp, minp, scale_k;  
  double dx,dy,dz, dist, mhsml;    
  MkLiuParticles &par = LiuParticles;  

  i= j= d=  sumiac= maxiac= miniac= niac = noiac= maxp= minp= scale_k=0;  
  dist= mhsml=0;      

//    Smoothing kernel function  
//    skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//        = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//        = 3, Quintic kernel (Morris 1997) 

//MkDebug("    direct_find()\n");
//MkDebug("    D1 ");

  if (LiuParam.skf==1) scale_k = 2;
  else scale_k = 3 ; 

  NIac = Pair_Count();
  LiuPairs.Initialize(NIac+1);
      
  for (i=0;i<NTotal+NVirt;i++) par[i].CountIac = 0; 
       
  for (i=0;i<NTotal+NVirt-1;i++) {    
    for (j = i+1;j< NTotal+NVirt;j++){ 

      dx = par[i].X - par[j].X; 
      dy = par[i].Y - par[j].Y; 
      dz = par[i].Z - par[j].Z; 
      dist = sqrt(dx*dx+dy*dy+dz*dz);
      mhsml = (par[i].SmoothLen+par[j].SmoothLen)/2.; 
      LiuKernel.SetSmoothLen(mhsml);
      if (dist < scale_k*mhsml) { 
	par[i].CountIac = par[i].CountIac + 1; 
	par[j].CountIac = par[j].CountIac + 1; 

	LiuPairs[niac].I = i; 
	LiuPairs[niac].J = j; 

	LiuPairs[niac].Dist = dist;
	LiuPairs[niac].dX = dx;
	LiuPairs[niac].dY = dy;
	LiuPairs[niac].dZ = dz;

//    Kernel and derivations of kernel 
	LiuPairs[niac].W = LiuKernel.W(dist);
	LiuPairs[niac].dWdX = LiuKernel.dWdX(dist,dx,dy,dz);
	LiuPairs[niac].dWdY = LiuKernel.dWdY(dist,dx,dy,dz);
	LiuPairs[niac].dWdZ = LiuKernel.dWdZ(dist,dx,dy,dz);
	niac = niac + 1;
      }
    }
  }
//MkDebug("D3 "); 
//    Statistics for the interaction 
 
  sumiac = 0; 
  maxiac = 0; 
  miniac = 1000; 
  noiac  = 0; 
  for (i=0;i<NTotal+NVirt;i++) { 
    sumiac = sumiac + par[i].CountIac; 
    if (par[i].CountIac>maxiac){ 
      maxiac = par[i].CountIac; 
      maxp = i; 
    }
    if (par[i].CountIac<miniac) {
      miniac = par[i].CountIac; 
      minp = i; 
    }
    if (par[i].CountIac==0) noiac  = noiac + 1;
  }
  
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    if (LiuParam.int_stat) {
      printf(" >> Statistics: interactions per particle:\n"); 
      printf("**** Particle: %d  maximal interactions: %d\n",maxp,maxiac) ;
      printf("**** Particle: %d minimal interactions: %d\n",minp,miniac) ;
      printf("**** Average :%f\n",double(sumiac)/double(NTotal+NVirt)) ;
      printf("**** Total pairs : %d\n",niac) ;
      printf("**** Particles with no interactions:%d\n",noiac) ;
    }
  }
//MkDebug("D4 "); 
//MkDebug("    ~direct_find()\n");
}

int MkLiuSPH::Pair_Grid_Count() 
{
  int i, j, niac=0, scale_k;  
  int I,J,K,ii,jj;
  double  dist, mhsml, dx,dy,dz;      
  MkLiuParticles &par = LiuParticles;

  if (LiuParam.skf==1) scale_k = 2;  
  else scale_k = 3; 

  for (I=0;I<LiuGrids.GetNX();I++) {
    for (J=0;J<LiuGrids.GetNY();J++) {
      for (K=0;K<LiuGrids.GetNZ();K++) {
	MkLiuGrid &grid = LiuGrids(I,J,K);
	MkInt & paref = LiuGrids.GetParticles(I,J,K);
	for (ii=0;ii<grid.GetNumOfParticle();ii++) {
	  for (jj=0;jj<paref.getSzX();jj++) {
	    i = grid.GetParticleRef()[ii];
	    j = paref[jj];
	    if (i>=j) continue;
	    
	    float Xi,Xj,Yi,Yj;
	    Xi = par[i].X;
	    Xj = par[j].X;
	    Yi = par[i].Y;
	    Yj = par[j].Y;

	    dx = par[i].X - par[j].X; 
	    dy = par[i].Y - par[j].Y; 
	    dz = par[i].Z - par[j].Z; 
	    dist = sqrt(dx*dx+dy*dy+dz*dz);
	    mhsml = (par[i].SmoothLen+par[j].SmoothLen)/2.; 
	    if (dist < scale_k*mhsml) {
	      //MkDebug("%d (%f,%f), %d(%f,%f) \n",i,Xi,Yi,j,Xj,Yj);
	      niac = niac + 1;  
	    }
	  }
	}
      }
    }
  }

//MkDebug("Number of pair count %d\n",niac);
  return niac;
}

void MkLiuSPH::Direct_Grid_Find() 
{ 
  int i, j, d,  sumiac, maxiac, miniac, niac, noiac, maxp, minp, scale_k;  
  int I, J, K, ii, jj;
  double dx,dy,dz, dist, mhsml;    
  MkLiuParticles &par = LiuParticles;  

  i= j= d=  sumiac= maxiac= miniac= niac = noiac= maxp= minp= scale_k=0;  
  dist= mhsml=0;      

//    Smoothing kernel function  
//    skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//        = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//        = 3, Quintic kernel (Morris 1997) 

  //MkDebug("    direct_grid_find()\n");
  //MkDebug("    DG1 ");

  LiuGrids.Update(LiuParticles);

  //MkDebug("    DG2 ");

  if (LiuParam.skf==1) scale_k = 2;
  else scale_k = 3 ; 

  NIac = Pair_Count();
  LiuPairs.Initialize(NIac+1);
  for (i=0;i<NTotal+NVirt;i++) par[i].CountIac = 0; 

  //MkDebug("    DG3 ");
  //MkDebug("NX %d, NY  %d, NZ %d \n",LiuGrids.GetNX(),LiuGrids.GetNY(),LiuGrids.GetNZ());

  for (I=0;I<LiuGrids.GetNX();I++) { 
    for (J=0;J<LiuGrids.GetNY();J++) {
      for (K=0;K<LiuGrids.GetNZ();K++) {
	//MkDebug("I %d J %d  K %d \n",I,J,K);
	MkLiuGrid &grid = LiuGrids(I,J,K);

	if (LiuGrids.GetNumOfParticles(I,J,K) <= 0) continue; 

	MkInt & paref = LiuGrids.GetParticles(I,J,K);
	MkDebug("paref has particles %d\n",paref.getSzX());
	for (ii=0;ii<grid.GetNumOfParticle();ii++) {
	  for (jj=0;jj<paref.getSzX();jj++) {
	    i = grid.GetParticleRef()[ii];
	    j = paref[jj];
	    MkDebug("i %d, j %d\n",i,j);
	    if (i>=j) continue;
	    MkDebug("    DG3-2 ");
	    dx = par[i].X - par[j].X; 
	    dy = par[i].Y - par[j].Y; 
	    dz = par[i].Z - par[j].Z; 
	    dist = sqrt(dx*dx+dy*dy+dz*dz);
	    mhsml = (par[i].SmoothLen+par[j].SmoothLen)/2.; 
	    LiuKernel.SetSmoothLen(mhsml);
	    if (dist < scale_k*mhsml) { 
	      par[i].CountIac = par[i].CountIac + 1; 
	      par[j].CountIac = par[j].CountIac + 1; 

	      LiuPairs[niac].I = i; 
	      LiuPairs[niac].J = j; 

	      LiuPairs[niac].Dist = dist;
	      LiuPairs[niac].dX = dx;
	      LiuPairs[niac].dY = dy;
	      LiuPairs[niac].dZ = dz;

//    Kernel and derivations of kernel 
	      LiuPairs[niac].W = LiuKernel.W(dist);
	      LiuPairs[niac].dWdX = LiuKernel.dWdX(dist,dx,dy,dz);
	      LiuPairs[niac].dWdY = LiuKernel.dWdY(dist,dx,dy,dz);
	      LiuPairs[niac].dWdZ = LiuKernel.dWdZ(dist,dx,dy,dz);
	      niac = niac + 1;
	    }
	  }
	}
      }
    }
  }
      
  //MkDebug("D4 "); 
//    Statistics for the interaction 
 
  sumiac = 0; 
  maxiac = 0; 
  miniac = 1000; 
  noiac  = 0; 
  for (i=0;i<NTotal+NVirt;i++) { 
    sumiac = sumiac + par[i].CountIac; 
    if (par[i].CountIac>maxiac){ 
      maxiac = par[i].CountIac; 
      maxp = i; 
    }
    if (par[i].CountIac<miniac) {
      miniac = par[i].CountIac; 
      minp = i; 
    }
    if (par[i].CountIac==0) noiac  = noiac + 1;
  }
  
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    if (LiuParam.int_stat) {
      printf(" >> Statistics: interactions per particle:\n"); 
      printf("**** Particle: %d  maximal interactions: %d\n",maxp,maxiac) ;
      printf("**** Particle: %d minimal interactions: %d\n",minp,miniac) ;
      printf("**** Average :%f\n",double(sumiac)/double(NTotal+NVirt)) ;
      printf("**** Total pairs : %d\n",niac) ;
      printf("**** Particles with no interactions:%d\n",noiac) ;
    }
  }
  //MkDebug("D5 "); 
  //MkDebug("    ~direct_grid_find()\n");
}

int MkLiuSPH::Pair_Bnd_Count() 
{
  int i, j, nbiac=0, scale_k;  
  double dist, mhsml;      
  MkLiuParticles &par = LiuParticles;
  MkLiuBoundarys &bnd = LiuBoundarys;

  if (LiuParam.skf==1) scale_k = 2;  
  else scale_k = 3; 

  for (i=0;i<NTotal;i++) {      
    for (j = 0;j< bnd.GetSize();j++){ 
      dist = bnd[j].GetDistance(par[i]);
      mhsml = (par[i].SmoothLen+bnd[j].SmoothLen)/2.; 
      if (dist < scale_k*mhsml) {
	nbiac = nbiac + 1;  
      }
    }
  }
//MkDebug("Number of pair count %d\n",nbiac);//usleep(999999);usleep(999999);usleep(999999);
  return nbiac;
}

void MkLiuSPH::Direct_Bnd_Find() 
{ 
  int i, j, d,  sumiac, maxiac, miniac, nbiac, noiac, maxp, minp, scale_k;  
  double dist, mhsml,dx,dy,dz;    
  MkLiuParticles &par = LiuParticles;  
  MkLiuBoundarys &bnd = LiuBoundarys;
  MkPoint pnt;

  i= j= d=  sumiac= maxiac= miniac= nbiac = noiac= maxp= minp= scale_k=0;  
  dist= mhsml=0;      

//    Smoothing kernel function  
//    skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//        = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//        = 3, Quintic kernel (Morris 1997) 

  //MkDebug("    direct_bnd_find()\n");
  //MkDebug("    DB1 ");

  if (LiuParam.skf==1) scale_k = 2;
  else scale_k = 3 ; 

  NBIac = Pair_Bnd_Count();
  LiuBndPairs.Initialize(NBIac+1);
      
  for (i=0;i<NTotal;i++) {    
    for (j = 0;j<bnd.GetSize();j++){ 

      dist = bnd[j].GetDistance(par[i]);
      pnt =  bnd[j].GetNearestPoint(par[i]);

      dx = par[i].X - pnt.X; 
      dy = par[i].Y - pnt.Y; 
      dz = par[i].Z - pnt.Z; 

      mhsml = (bnd[j].SmoothLen+par[i].SmoothLen)/2.; 
      LiuKernel.SetSmoothLen(mhsml);
      if (dist < scale_k*mhsml) { 
	//MkDebug("%d, X %f Y %f Z %f dist %f\n",i,par[i].X, par[i].Y, par[i].Z,dist);
	par[i].CountIac++; 

	LiuBndPairs[nbiac].I = i; //particle
	LiuBndPairs[nbiac].J = j; //boundary

	LiuBndPairs[nbiac].Dist = dist;
	LiuBndPairs[nbiac].dX = dx;
	LiuBndPairs[nbiac].dY = dy;
	LiuBndPairs[nbiac].dZ = dz;

//    Kernel and derivations of kernel 
	LiuBndPairs[nbiac].W = LiuKernel.W(dist);
	LiuBndPairs[nbiac].dWdX = LiuKernel.dWdX(dist,dx,dy,dz);
	LiuBndPairs[nbiac].dWdY = LiuKernel.dWdY(dist,dx,dy,dz);
	LiuBndPairs[nbiac].dWdZ = LiuKernel.dWdZ(dist,dx,dy,dz);
	nbiac = nbiac + 1;
      }
    }
  }
//MkDebug("DB3 "); 
//    Statistics for the interaction 
 
  sumiac = 0; 
  maxiac = 0; 
  miniac = 1000; 
  noiac  = 0; 
  for (i=0;i<NTotal;i++) { 
    sumiac = sumiac + par[i].CountIac; 
    if (par[i].CountIac>maxiac){ 
      maxiac = par[i].CountIac; 
      maxp = i; 
    }
    if (par[i].CountIac<miniac) {
      miniac = par[i].CountIac; 
      minp = i; 
    }
    if (par[i].CountIac==0) noiac  = noiac + 1;
  }
  
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    if (LiuParam.int_stat) {
      printf(" >> Statistics: interactions per particle:\n"); 
      printf("**** Particle: %d  maximal interactions: %d\n",maxp,maxiac) ;
      printf("**** Particle: %d minimal interactions: %d\n",minp,miniac) ;
      printf("**** Average :%f\n",double(sumiac)/double(NTotal+NVirt)) ;
      printf("**** Total pairs : %d\n",nbiac) ;
      printf("**** Particles with no interactions:%d\n",noiac) ;
    }
  }
//MkDebug("DB4 "); 
//MkDebug("    ~direct_bnd_find()\n");
}


//---------------------------------------------------------------------- 
// Subroutine to determine the information of virtual particles 
// Here only the Monaghan type virtual particles for the 2D shear 
// cavity driven problem are generated. 
//   CurrentTimeStep : Current time step                                 [in] 
//   NTotal : Number of particles                                  [in] 
//   nvirt  : Number of virtual particles                         [out] 
//   Hsml   : Smoothing Length                                 [in|out] 
//   Mass   : Particle masses                                  [in|out] 
//   X      : Coordinates of all particles                     [in|out] 
//   VX     : Velocities of all particles                      [in|out] 
//   Rho    : Density                                          [in|out] 
//   U      : internal energy                                  [in|out] 
//   ParticleType   : type of particles                               [in|out] 
 
void MkLiuSPH::Virt_Part()  // should delete it later 
{
  int i, j, d, im, mp,np,nvirt; 
  double xl, dx, v_inf, rho_ref=1000; 
  FILE *fp1,*fp2, *fp3;
  MkLiuParticles &par = LiuParticles;

   i= j= d= im= mp=0; 
   xl= dx= v_inf=0; 
 
//MkDebug("    virt_part()\n");

  if (LiuParam.vp_input){
    //    LoadVPFromFile();
  }
  else  {
    nvirt = 0; 
    mp = 20;
    np = 20;
    xl = 1.0; 
    dx = xl / mp; 
    v_inf = 0;//1.e-3; 

//   Monaghan type virtual particle on the Upper side 
  //MkDebug("    virt_part()VP1 nvirt is %d ",nvirt); 
    for (i = 0;i< 2*mp+1;i++) {
      par[NTotal + nvirt].X = i*dx/2;
      par[NTotal + nvirt].Y = xl+2*dx;   
      par[NTotal + nvirt].XVel = v_inf; 
      par[NTotal + nvirt].YVel = 0.; 
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_part()VP2 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Lower side 
 
    for (i = 0;i< 2*mp+1;i++) {
      par[NTotal + nvirt].X = i*dx/2;  
      par[NTotal + nvirt].Y = 0.-2*dx;   
      par[NTotal + nvirt].XVel = 0.; 
      par[NTotal + nvirt].YVel = 0.; 
      nvirt = nvirt + 1; 
    }
  //MkDebug("    virt_part()VP3 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Left side 
 
    for (i = 0;i< 2*np-1;i++) {
      par[NTotal + nvirt].X = 0.-2*dx;
      par[NTotal + nvirt].Y = (i+1)*dx/2;
      par[NTotal + nvirt].XVel = 0.; 
      par[NTotal + nvirt].YVel = 0.;
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_part()VP4 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Right side 
 
    for (i = 0;i< 2*np-1;i++) {
      par[NTotal + nvirt].X = xl+2*dx;
      par[NTotal + nvirt].Y = (i+1)*dx/2;
      par[NTotal + nvirt].XVel = 0.; 
      par[NTotal + nvirt].YVel = 0.; 
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_part()VP5 nvirt is %d ",nvirt); 
    for (i = 0;i< nvirt;i++) {
      par[NTotal + i].Rho = rho_ref;
      par[NTotal + i].Mass = 0.01*(par[NTotal + i].Rho  * dx * dx); 
      par[NTotal + i].Eta = 100000000000;
      par[NTotal + i].Press = 0.; 
      par[NTotal + i].Energy = 357.1; 
      par[NTotal + i].ParticleType = -2; 
      par[NTotal + i].SmoothLen = dx; 
      par[NTotal + i].Radius = dx/2.0;
    }
  } 
  if ((CurrentTimeStep%LiuParam.save_step)==0) {

  }
 //MkDebug("    virt_part()VP6 nvirt is %d ",nvirt); 
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    if (LiuParam.int_stat) {
      printf(" >> Statistics: Virtual boundary particles:\n"); 
      printf("          Number of virtual particles:%d\n",NVirt);
    }
  }

//MkDebug("    ~virt_part()\n");
}
//---------------------------------------------------------------------- 
// Subroutine to define the fluid particle viscosity 
  
//   NTotal  : Number of particles                                 [in] 
//   IType    : Type of particle                                   [in] 
//   X       : Coordinates of all particles                        [in] 
//   Rho     : Density                                             [in] 
//   Eta     : Dynamic viscosity                                  [out] 

void MkLiuSPH::Viscosity() 
{       
  int i; 
  static double vis = 1.0e-3;
  MkLiuParticles &par = LiuParticles;

  for (i=0;i<NTotal+NVirt;i++) {
    if (par[i].GetParticleType()==1) par[i].SetEta(0.0);
    else if (par[i].GetParticleType()==2) par[i].SetEta(vis);
  }
}

//---------------------------------------------------------------------- 
//  Gamma law EOS: subroutine to calculate the pressure and sound   
  
//  rho    : Density                                              [in] 
//  u      : Internal energy                                      [in] 
//  p      : Pressure                                            [out] 
//  c      : sound velocity                                      [out] 

double MkLiuSPH::P_ideal_gas(double rho, double u) 
{       
  double gamma=1.4;
           
//   For air (idea gas) 
//   See Equ.(3.82) 
 
  return (gamma-1) * rho * u;
}

double MkLiuSPH::C_ideal_gas(double rho, double u) 
{       
  double gamma=1.4;
           
//   For air (idea gas) 
//   See Equ.(3.82) 
 
  return sqrt((gamma-1) * u);
}

//---------------------------------------------------------------------- 
//   Artificial equation of state for the artificial compressibility  
 
//  rho    : Density                                              [in] 
//  u      : Internal energy                                      [in] 
//  p      : Pressure                                            [out] 
//  c      : sound velocity                                      [out] 
//  Equation of state for artificial compressibility    
    
double MkLiuSPH::P_art_water(double rho) 
{       
//  Artificial EOS, Form 1 (Monaghan, 1994)  
//  See Equ.(4.88) 
  double gamma=7.; 
  double rho0=Rho_Ref;
  double b = 20;//2000; 
  double p =b*(pow(rho/rho0,gamma)-1);
  //MkDebug("P_art_water() %12.6f %12.6f %12.6f \n",rho,rho0,p);
 
//  Artificial EOS, Form 2 (Morris, 1997) 
//  See Equ.(4.89) 
//  double p,c = 0.01; 
//  p = c*c * rho;
  return p;
}

double MkLiuSPH::C_art_water(double rho) 
{       
//  Artificial EOS, Form 1 (Monaghan, 1994)  
//  See Equ.(4.88) 
  double c = 1480.0;
 
//  Artificial EOS, Form 2 (Morris, 1997) 
//  See Equ.(4.89) 
//  double c = 0.01; 
  return c;
}

 
