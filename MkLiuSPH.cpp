#include "MkLiuSPH.hpp"

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
  Dim = 0;
  NVirt=0;
  NTotal=0;
      
  MaxInteraction = 0;
  NIac=0;        

  MaxTimeStep=0; 
  CurrentTimeStep=0;   
  Dt=0;       

  LiuPairs.Clear();    

  LiuParticles.Clear();
  LiuKernel.Clear();
}

void MkLiuSPH::Draw()
{
  LiuParticles.Draw();
} 
 
void MkLiuSPH::Run()
{
  static int i, j, k, d;
  static double time, temp_rho, temp_u;
  static bool is_first = true;

  int maxn=NTotal+NVirt;
//MkDebug("MkLiuSPH::Run() Dim %d, maxn %d\n",Dim, maxn);

  MkDouble v_min(Dim, maxn), u_min(maxn), rho_min(maxn); 

//MkDebug("MkLiuSPH::Run()\n");
  if(is_first) {
    i = j = k = d = 0;
    time = temp_rho = temp_u = 0;

    for (i = 0;i< NTotal;i++) {
      LiuParticles[i].XVel[0] = 0;
      LiuParticles[i].YVel[0] = 0;
      LiuParticles[i].ZVel[0] = 0;
    }
  }
//MkDebug("MkLiuSPH::Run() #1 \n");
  if (CurrentTimeStep> MaxTimeStep) {is_first = false;return;}
  //  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    printf("______________________________________________\n");
    printf("  current number of time step = %d     current time=%f\n",CurrentTimeStep, double(time+Dt));
    printf("______________________________________________\n");
    //  }
//MkDebug("MkLiuSPH::Run() #2 \n");
//   If not first time step, then update thermal energy, density and  
//   velocity half a time step   

  if (!is_first) {
    for (i = 0;i< NTotal;i++) {
      u_min(i) = LiuParticles[i].Energy;
      temp_u=0.;
      if (Dim==1) temp_u=-LiuParam.nsym*LiuParticles[i].Press*LiuParticles[i].XVel[0]/LiuParticles[i].X/LiuParticles[i].Rho;
      LiuParticles[i].Energy = LiuParticles[i].Energy + (Dt/2.)* (LiuParticles[i].DUDt+temp_u);
      if(LiuParticles[i].Energy < 0) LiuParticles[i].Energy = 0.;
             
      if (!LiuParam.summation_density) {
	rho_min(i) = LiuParticles[i].Rho;
	temp_rho=0.;
	if (Dim==1) temp_rho=-LiuParam.nsym*LiuParticles[i].Rho*LiuParticles[i].XVel[0]/LiuParticles[i].X;
	LiuParticles[i].Rho = LiuParticles[i].Rho +(Dt/2.)*( LiuParticles[i].DRhoDt+ temp_rho);
      }

      v_min(0, i) = LiuParticles[i].XVel[0];
      LiuParticles[i].XVel[0] = LiuParticles[i].XVel[0] + (Dt/2.)*LiuParticles[i].DVXDt;
      v_min(1, i) = LiuParticles[i].YVel[0];
      LiuParticles[i].YVel[0] = LiuParticles[i].YVel[0] + (Dt/2.)*LiuParticles[i].DVYDt;
      v_min(2, i) = LiuParticles[i].ZVel[0];
      LiuParticles[i].ZVel[0] = LiuParticles[i].ZVel[0] + (Dt/2.)*LiuParticles[i].DVZDt;
    }
  }
//---  Definition of variables out of the function vector:     

//MkDebug("MkLiuSPH::Run() #3 \n");
  Single_Step();
//MkDebug("MkLiuSPH::Run() #4 \n");  
  if (is_first) {
   
    for (i=0;i<NTotal;i++) {
      temp_u=0.;
      if (Dim==1) temp_u=-LiuParam.nsym*LiuParticles[i].Press*LiuParticles[i].XVel[0]/LiuParticles[i].X/LiuParticles[i].Rho;
      LiuParticles[i].Energy = LiuParticles[i].Energy + (Dt/2.)*(LiuParticles[i].DUDt + temp_u);
      if(LiuParticles[i].Energy<0)  LiuParticles[i].Energy = 0.;
          
      if (!LiuParam.summation_density ){
	temp_rho=0.;
	if (Dim==1) temp_rho=-LiuParam.nsym*LiuParticles[i].Rho*LiuParticles[i].XVel[0]/LiuParticles[i].X;
	LiuParticles[i].Rho = LiuParticles[i].Rho + (Dt/2.)* (LiuParticles[i].DRhoDt+temp_rho);
      }
          
      LiuParticles[i].XVel[0] = LiuParticles[i].XVel[0] + (Dt/2.) * LiuParticles[i].DVXDt + LiuParticles[i].XAVel[0];
      LiuParticles[i].X = LiuParticles[i].X + Dt * LiuParticles[i].XVel[0];
      LiuParticles[i].YVel[0] = LiuParticles[i].YVel[0] + (Dt/2.) * LiuParticles[i].DVYDt + LiuParticles[i].YAVel[0];
      LiuParticles[i].Y = LiuParticles[i].Y + Dt * LiuParticles[i].YVel[0];
      LiuParticles[i].ZVel[0] = LiuParticles[i].ZVel[0] + (Dt/2.) * LiuParticles[i].DVZDt + LiuParticles[i].ZAVel[0];
      LiuParticles[i].Z = LiuParticles[i].Z + Dt * LiuParticles[i].ZVel[0];

    }
  }               
  else {
    for (i=0;i<NTotal;i++) {
      temp_u=0.;
      if (Dim==1) temp_u=-LiuParam.nsym*LiuParticles[i].Press*LiuParticles[i].XVel[0]/LiuParticles[i].X/LiuParticles[i].Rho;
      LiuParticles[i].Energy = u_min(i) + Dt*(LiuParticles[i].DUDt+temp_u);
      if(LiuParticles[i].Energy<0)  LiuParticles[i].Energy = 0.;
             
      if (!LiuParam.summation_density ) {
	temp_rho=0.;
	if (Dim==1) temp_rho=-LiuParam.nsym*LiuParticles[i].Rho*LiuParticles[i].XVel[0]/LiuParticles[i].X;
	LiuParticles[i].Rho = rho_min(i) + Dt*(LiuParticles[i].DRhoDt+temp_rho);
      }
                 
      LiuParticles[i].XVel[0] = v_min(0, i) + Dt * LiuParticles[i].DVXDt + LiuParticles[i].XAVel[0];
      LiuParticles[i].X = LiuParticles[i].X + Dt * LiuParticles[i].XVel[0];
      LiuParticles[i].YVel[0] = v_min(d, i) + Dt * LiuParticles[i].DVYDt + LiuParticles[i].YAVel[0];
      LiuParticles[i].Y = LiuParticles[i].Y + Dt * LiuParticles[i].YVel[0];
      LiuParticles[i].ZVel[0] = v_min(d, i) + Dt * LiuParticles[i].DVZDt + LiuParticles[i].ZAVel[0];
      LiuParticles[i].Z = LiuParticles[i].Z + Dt * LiuParticles[i].ZVel[0];
    }
         
  }
//MkDebug("MkLiuSPH::Run() #5 \n"); 
  if ((CurrentTimeStep%LiuParam.save_step)==0)  Output();
  //  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    printf("\n");
    //      123456789ABC 123456789ABC 123456789ABC 123456789ABC  
    printf("      x           y         velocity       dvx\n");
    printf("%12.6f %12.6f %12.6f %12.6f \n",LiuParticles[LiuParam.moni_particle].X, LiuParticles[LiuParam.moni_particle].Y,LiuParticles[LiuParam.moni_particle].XVel[0],LiuParticles[LiuParam.moni_particle].DVXDt);
    //  }

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
 
//MkDebug("\n  single_step()\n");
//MkDebug("  S1 ");
  for ( i=0;i<NTotal;i++) {
    LiuParticles[i].AVDUDt = 0.; 
    LiuParticles[i].AHDUDt = 0.; 
    LiuParticles[i].INDVXDt = 0.;
    LiuParticles[i].ARDVXDt = 0.; 
    LiuParticles[i].EXDVXDt = 0.;
    LiuParticles[i].INDVYDt = 0.;
    LiuParticles[i].ARDVYDt = 0.; 
    LiuParticles[i].EXDVYDt = 0.;
    LiuParticles[i].INDVZDt = 0.;
    LiuParticles[i].ARDVZDt = 0.; 
    LiuParticles[i].EXDVZDt = 0.;
  }
//MkDebug("  S2 ");  
//---  Positions of virtual (boundary) particles:  
 
  if (LiuParam.virtual_part) {
    Virt_Part(); 
  }
//MkDebug("  S3 ");      
//---  Interaction parameters, calculating neighboring particles 
//   and optimzing smoothing length 
 
  if (LiuParam.nnps==1) Direct_Find();

  //MkDebug("  S4 ");
//---  Density approximation or change rate 

  if (LiuParam.summation_density) Sum_Density();
  else Con_Density();
 
//---  Dynamic viscosity: 
 
//MkDebug("  S5 ");
  if (LiuParam.visc) Viscosity(); 
        
//---  Internal forces: 
  
  Int_Force(); 
               
//MkDebug("  S6 ");    

//---  Artificial viscosity: 
 
  if (LiuParam.visc_artificial) Art_Visc();
//MkDebug("  S7 ");
       
//---  External forces: 
 
//MkDebug("  S8 ");

  if (LiuParam.ex_force) Ext_Force();
 
//MkDebug("  S9 ");

//   Calculating the neighboring particles and undating DATA.HSML 
       
  if (LiuParam.sle!=0) H_Upgrade(); 
 
  if (LiuParam.heat_artificial) Art_Heat();
      
//   Calculating average velocity of each partile for avoiding penetration 
 
  if (LiuParam.average_velocity) Av_Vel();
 
//---  Convert velocity, force, and energy to f and dfdt   
 
  for (i=0;i<NTotal;i++) {
    LiuParticles[i].DVXDt = LiuParticles[i].INDVXDt + 
      LiuParticles[i].EXDVXDt + LiuParticles[i].ARDVXDt;
    LiuParticles[i].DVYDt = LiuParticles[i].INDVYDt + 
      LiuParticles[i].EXDVYDt + LiuParticles[i].ARDVYDt;
    LiuParticles[i].DVZDt = LiuParticles[i].INDVZDt + 
      LiuParticles[i].EXDVZDt + LiuParticles[i].ARDVZDt;
    LiuParticles[i].DUDt  = LiuParticles[i].DUDt + 
      LiuParticles[i].INDUDt + LiuParticles[i].AVDUDt + LiuParticles[i].AHDUDt;
  }
  // //MkDebug("  S10 ");
  if ((CurrentTimeStep%LiuParam.print_step)==0) {
    //      123456789abc 123456789abc 123456789abc
    printf("\n") ;
    printf("**** Information for particle **** %d\n",LiuParam.moni_particle);
    printf("internal a   artifical a  external a    total a \n");
    printf("%12.6f %12.6f %12.6f %12.6f\n",LiuParticles[LiuParam.moni_particle].INDVXDt,LiuParticles[LiuParam.moni_particle].ARDVXDt,LiuParticles[LiuParam.moni_particle].EXDVXDt,LiuParticles[LiuParam.moni_particle].DVXDt);
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
  MkDouble vcc(NTotal),dvx(Dim);

  i=j=k=d =0;
  dx= vr= rr= h= mc= mrho= mhsml= hvcc= mui= muj= muij= rdwdx= g1=g2=0;
       
  //---  Parameter for the artificial heat conduction:
      
  g1=0.1; 
  g2=1.0;
  for (i=0;i<NTotal+NVirt;i++) { 
    vcc(i) = 0.e0;
    LiuParticles[i].AHDUDt = 0.e0; 
  }
     
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 

    dvx(0) = LiuParticles[j].XVel[0] - LiuParticles[i].XVel[0];
    dvx(1) = LiuParticles[j].YVel[0] - LiuParticles[i].YVel[0];
    dvx(2) = LiuParticles[j].ZVel[0] - LiuParticles[i].ZVel[0];

    hvcc = dvx(0)*LiuPairs[k].DWDX(1-1);

    for (d=1;d<Dim;d++) { 
      hvcc = hvcc + dvx(d)*LiuPairs[k].DWDX(d) ;
    }    
    vcc(i) = vcc(i) + LiuParticles[j].Mass*hvcc/LiuParticles[j].Rho; 
    vcc(j) = vcc(j) + LiuParticles[i].Mass*hvcc/LiuParticles[i].Rho; 
  }
    
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    mhsml= (LiuParticles[i].SmoothLen+LiuParticles[j].SmoothLen)/2.; 
    LiuKernel.SetSmoothLen(mhsml);

    mrho = 0.5e0*(LiuParticles[i].Rho + LiuParticles[j].Rho);
    rr = 0.e0; 
    rdwdx = 0.e0; 

    dx = LiuParticles[i].X -  LiuParticles[j].X;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*LiuPairs[k].DWDX(d-1);
    dx = LiuParticles[i].Y -  LiuParticles[j].Y;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*LiuPairs[k].DWDX(d-1);
    dx = LiuParticles[i].Z -  LiuParticles[j].Z;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*LiuPairs[k].DWDX(d-1);


    mui=g1*LiuParticles[i].SmoothLen*LiuParticles[i].SoundSpeed + g2*LiuParticles[i].SmoothLen*LiuParticles[i].SmoothLen*(fabs(vcc(i))-vcc(i));
    muj=g1*LiuParticles[j].SmoothLen*LiuParticles[j].SoundSpeed + g2*LiuParticles[j].SmoothLen*LiuParticles[j].SmoothLen*(fabs(vcc(j))-vcc(j));
    muij= 0.5*(mui+muj);
    h = muij/(mrho*(rr+0.01*mhsml*mhsml))*rdwdx ;
    LiuParticles[i].AHDUDt = LiuParticles[i].AHDUDt + LiuParticles[j].Mass*h*(LiuParticles[i].Energy-LiuParticles[j].Energy) ;
    LiuParticles[j].AHDUDt = LiuParticles[j].AHDUDt + LiuParticles[i].Mass*h*(LiuParticles[j].Energy-LiuParticles[i].Energy); 
  }
 
  for (i=0;i<NTotal+NVirt;i++) { 
    LiuParticles[i].AHDUDt = 2.0e0*LiuParticles[i].AHDUDt;           
    MkDebug("%d, %12.6f\n",i,LiuParticles[i].AHDUDt);
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
  MkDouble dvx(Dim);

   i=j=k=d =0;
   dx= alpha= beta= etq= piv= muv= vr= rr= h= mc= mrho= mhsml=0;

//  Parameter for the artificial viscosity: 
//  Shear viscosity 
  alpha = 1.e0; 
      
//  Bulk viscosity 
  beta  = 1.e0;  
       
//  Parameter to avoid singularities 
  etq   = 0.1e0;
            
  for (i=0;i<NTotal+NVirt;i++){ 
    LiuParticles[i].ARDVXDt = 0.e0;
    LiuParticles[i].ARDVYDt = 0.e0;
    LiuParticles[i].ARDVZDt = 0.e0;
    LiuParticles[i].AVDUDt = 0.e0;
  }
      
//  Calculate SPH sum for artificial viscosity 
       
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    mhsml= (LiuParticles[i].SmoothLen+LiuParticles[j].SmoothLen)/2.; 
    LiuKernel.SetSmoothLen(mhsml);
    vr = 0.e0; 
    rr = 0.e0; 

    dvx(0) = LiuParticles[i].XVel[0] - LiuParticles[j].XVel[0]; 
    dx     = LiuParticles[i].X - LiuParticles[j].X; 
    vr     = vr + dvx(0)*dx; 
    rr     = rr + dx*dx; 
    dvx(1) = LiuParticles[i].YVel[0] - LiuParticles[j].YVel[0]; 
    dx     = LiuParticles[i].Y - LiuParticles[j].Y; 
    vr     = vr + dvx(1)*dx; 
    rr     = rr + dx*dx; 
    dvx(2) = LiuParticles[i].ZVel[0] - LiuParticles[j].ZVel[0]; 
    dx     = LiuParticles[i].Z - LiuParticles[j].Z; 
    vr     = vr + dvx(2)*dx; 
    rr     = rr + dx*dx; 

//  Artificial viscous force only if v_ij * r_ij  0 
 
    if (vr<0.e0) { 
 
//  Calculate muv_ij = data.Hsml v_ij * r_ij / ( r_ij^2 + data.Hsml^2 etq^2 ) 
             
      muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq);
           
//  Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij 
 
      mc   = 0.5e0*(LiuParticles[i].SoundSpeed + LiuParticles[j].SoundSpeed);
      mrho = 0.5e0*(LiuParticles[i].Rho + LiuParticles[j].Rho); 
      piv  = (beta*muv - alpha*mc)*muv/mrho;
 
//  Calculate SPH sum for artificial viscous force 
 
      h = -piv*LiuPairs[k].DWDX(1-1);
      LiuParticles[i].ARDVXDt = LiuParticles[i].ARDVXDt + LiuParticles[j].Mass*h;
      LiuParticles[j].ARDVXDt = LiuParticles[j].ARDVXDt - LiuParticles[i].Mass*h;
      LiuParticles[i].AVDUDt = LiuParticles[i].AVDUDt - LiuParticles[j].Mass*dvx(0)*h;
      LiuParticles[j].AVDUDt = LiuParticles[j].AVDUDt - LiuParticles[i].Mass*dvx(0)*h;
      h = -piv*LiuPairs[k].DWDX(2-1);
      LiuParticles[i].ARDVYDt = LiuParticles[i].ARDVYDt + LiuParticles[j].Mass*h;
      LiuParticles[j].ARDVYDt = LiuParticles[j].ARDVYDt - LiuParticles[i].Mass*h;
      LiuParticles[i].AVDUDt = LiuParticles[i].AVDUDt - LiuParticles[j].Mass*dvx(1)*h;
      LiuParticles[j].AVDUDt = LiuParticles[j].AVDUDt - LiuParticles[i].Mass*dvx(1)*h;
      h = -piv*LiuPairs[k].DWDX(3-1);
      LiuParticles[i].ARDVZDt = LiuParticles[i].ARDVZDt + LiuParticles[j].Mass*h;
      LiuParticles[j].ARDVZDt = LiuParticles[j].ARDVZDt - LiuParticles[i].Mass*h;
      LiuParticles[i].AVDUDt = LiuParticles[i].AVDUDt - LiuParticles[j].Mass*dvx(2)*h;
      LiuParticles[j].AVDUDt = LiuParticles[j].AVDUDt - LiuParticles[i].Mass*dvx(2)*h;

    }
  }
 
//  Change of specific internal energy: 
 
  for (i=0;i<NTotal+NVirt;i++) { 
    LiuParticles[i].AVDUDt = 0.5e0*LiuParticles[i].AVDUDt;
    MkDebug("i:%d, ARDVXDt:%12.6f",i,LiuParticles[i].ARDVXDt );
  }
  getch();
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
    LiuParticles[i].TDSDt = 0.e0; 
    LiuParticles[i].INDUDt = 0.e0; 
    LiuParticles[i].INDVXDt = 0.e0;
    LiuParticles[i].INDVYDt = 0.e0;
    LiuParticles[i].INDVZDt = 0.e0;
  }
 
//   Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c 
 
  if (LiuParam.visc) {
    for (k=0;k<NIac;k++) { 
      i = LiuPairs[k].I; 
      j = LiuPairs[k].J; 
      dvx(0) = LiuParticles[j].XVel[0] - LiuParticles[i].XVel[0];
      dvx(1) = LiuParticles[j].YVel[0] - LiuParticles[i].YVel[0];
      dvx(2) = LiuParticles[j].ZVel[0] - LiuParticles[i].ZVel[0];

      if (Dim==1) {
	hxx = 2.e0*dvx(1)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==2) {
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) -  dvx(1)*LiuPairs[k].DWDX(2-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - dvx(0)*LiuPairs[k].DWDX(1-1);
      }
      else if (Dim==3){
	hxx = 2.e0*dvx(0)*LiuPairs[k].DWDX(1-1) - dvx(1)*LiuPairs[k].DWDX(2-1)- dvx(2)*LiuPairs[k].DWDX(3-1);
	hxy = dvx(0)*LiuPairs[k].DWDX(2-1) + dvx(1)*LiuPairs[k].DWDX(1-1);
	hxz = dvx(0)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(1-1);
	hyy = 2.e0*dvx(1)*LiuPairs[k].DWDX(2-1) - dvx(0)*LiuPairs[k].DWDX(1-1)- dvx(2)*LiuPairs[k].DWDX(3-1);
	hyz = dvx(1)*LiuPairs[k].DWDX(3-1) + dvx(2)*LiuPairs[k].DWDX(2-1);
	hzz = 2.e0*dvx(2)*LiuPairs[k].DWDX(3-1) - dvx(0)*LiuPairs[k].DWDX(1-1) - dvx(1)*LiuPairs[k].DWDX(2-1);
      }
      hxx = 2.e0/3.e0*hxx;
      hyy = 2.e0/3.e0*hyy;
      hzz = 2.e0/3.e0*hzz;
      if (Dim==1) {
	txx(i) = txx(i) + LiuParticles[j].Mass*hxx/LiuParticles[j].Rho;
	txx(j) = txx(j) + LiuParticles[i].Mass*hxx/LiuParticles[i].Rho;
      }
      else if (Dim==2) {
	txx(i) = txx(i) + LiuParticles[j].Mass*hxx/LiuParticles[j].Rho;
	txx(j) = txx(j) + LiuParticles[i].Mass*hxx/LiuParticles[i].Rho;
	txy(i) = txy(i) + LiuParticles[j].Mass*hxy/LiuParticles[j].Rho;
	txy(j) = txy(j) + LiuParticles[i].Mass*hxy/LiuParticles[i].Rho;             
	tyy(i) = tyy(i) + LiuParticles[j].Mass*hyy/LiuParticles[j].Rho; 
	tyy(j) = tyy(j) + LiuParticles[i].Mass*hyy/LiuParticles[i].Rho;
      }
      else if (Dim==3) {
	txx(i) = txx(i) + LiuParticles[j].Mass*hxx/LiuParticles[j].Rho;
	txx(j) = txx(j) + LiuParticles[i].Mass*hxx/LiuParticles[i].Rho;    
	txy(i) = txy(i) + LiuParticles[j].Mass*hxy/LiuParticles[j].Rho;
	txy(j) = txy(j) + LiuParticles[i].Mass*hxy/LiuParticles[i].Rho;
	txz(i) = txz(i) + LiuParticles[j].Mass*hxz/LiuParticles[j].Rho;
	txz(j) = txz(j) + LiuParticles[i].Mass*hxz/LiuParticles[i].Rho;
	tyy(i) = tyy(i) + LiuParticles[j].Mass*hyy/LiuParticles[j].Rho;
	tyy(j) = tyy(j) + LiuParticles[i].Mass*hyy/LiuParticles[i].Rho;
	tyz(i) = tyz(i) + LiuParticles[j].Mass*hyz/LiuParticles[j].Rho;
	tyz(j) = tyz(j) + LiuParticles[i].Mass*hyz/LiuParticles[i].Rho;
	tzz(i) = tzz(i) + LiuParticles[j].Mass*hzz/LiuParticles[j].Rho;
	tzz(j) = tzz(j) + LiuParticles[i].Mass*hzz/LiuParticles[i].Rho;
      }
 
//   Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz: 
 
      hvcc = 0.;
      for (d=0;d<Dim;d++) {
	hvcc = hvcc + dvx(d)*LiuPairs[k].DWDX(d);
      }
      vcc(i) = vcc(i) + LiuParticles[j].Mass*hvcc/LiuParticles[j].Rho;
      vcc(j) = vcc(j) + LiuParticles[i].Mass*hvcc/LiuParticles[i].Rho;
    }
  }
 
  for (i=0;i<NTotal+NVirt;i++) {
//   Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab 
    if (LiuParam.visc) {
      if (Dim==1) {
	LiuParticles[i].TDSDt = txx(i)*txx(i);
      }
      else if (Dim==2) {
	LiuParticles[i].TDSDt = txx(i)*txx(i) + 2.e0*txy(i)*txy(i)+ tyy(i)*tyy(i);
      }
      else if (Dim==3) {
	LiuParticles[i].TDSDt = txx(i)*txx(i)+ tyy(i)*tyy(i)+ tzz(i)*tzz(i)
                      + 2.e0*txy(i)*txy(i)+ 2.e0*txz(i)*txz(i)+ 2.e0*tyz(i)*tyz(i)  ;
      }
      LiuParticles[i].TDSDt = 0.5e0*LiuParticles[i].Eta/LiuParticles[i].Rho*LiuParticles[i].TDSDt;
    }
 
//   Pressure from equation of state 
 
    if (abs(LiuParticles[i].ParticleType)==1) {
      LiuParticles[i].Press = P_ideal_gas(LiuParticles[i].Rho, LiuParticles[i].Energy);
      LiuParticles[i].SoundSpeed = C_ideal_gas(LiuParticles[i].Rho, LiuParticles[i].Energy);
    }
    else if (abs(LiuParticles[i].ParticleType)==2) {
      LiuParticles[i].Press = P_art_water(LiuParticles[i].Rho);
      LiuParticles[i].SoundSpeed = C_art_water(LiuParticles[i].Rho);
    }
  } 
//    Calculate SPH sum for pressure force -p,a/rho 
//    and viscous force (eta Tab),b/rho 
//    and the internal energy change de/dt due to -p/rho vc,c 
 
  for (k=0;k<NIac;k++) {
    i = LiuPairs[k].I;
    j = LiuPairs[k].J;
    he = 0.e0;
       
//   For SPH algorithm 1 
 
    rhoij = 1.e0/(LiuParticles[i].Rho*LiuParticles[j].Rho);
    if(LiuParam.pa_sph==1){
      
      for (d=0;d<Dim ;d++) {
         
//   Pressure part 
                     
	h = -(LiuParticles[i].Press + LiuParticles[j].Press)*LiuPairs[k].DWDX(d);
	if (d==0) he = he + (LiuParticles[j].XVel[0] - LiuParticles[i].XVel[0])*h;
	if (d==1) he = he + (LiuParticles[j].YVel[0] - LiuParticles[i].YVel[0])*h;
	if (d==2) he = he + (LiuParticles[j].ZVel[0] - LiuParticles[i].ZVel[0])*h;
 
//   Viscous force 
 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txx(i) + LiuParticles[j].Eta*txx(j))*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2) {
	      h = h + (LiuParticles[i].Eta*txy(i) + LiuParticles[j].Eta*txy(j))*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h = h + (LiuParticles[i].Eta*txz(i) + LiuParticles[j].Eta*txz(j))*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txy(i) + LiuParticles[j].Eta*txy(j))*LiuPairs[k].DWDX(1-1) 
	          + (LiuParticles[i].Eta*tyy(i) + LiuParticles[j].Eta*tyy(j))*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (LiuParticles[i].Eta*tyz(i) + LiuParticles[j].Eta*tyz(j))*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txz(i) + LiuParticles[j].Eta*txz(j))*LiuPairs[k].DWDX(1-1) 
	          + (LiuParticles[i].Eta*tyz(i) + LiuParticles[j].Eta*tyz(j))*LiuPairs[k].DWDX(2-1)
	          + (LiuParticles[i].Eta*tzz(i) + LiuParticles[j].Eta*tzz(j))*LiuPairs[k].DWDX(3-1);
	  }
	}
	h = h*rhoij;
	//        MkDebug("Int_Force i:%d, j:%d, h:%f, rhoij:%f\n",i,j,h,rhoij);

	if (d==0) { 
	  LiuParticles[i].INDVXDt = LiuParticles[i].INDVXDt + LiuParticles[j].Mass*h;
	  LiuParticles[j].INDVXDt = LiuParticles[j].INDVXDt - LiuParticles[i].Mass*h;
	}
	if (d==1) { 
	  LiuParticles[i].INDVYDt = LiuParticles[i].INDVYDt + LiuParticles[j].Mass*h;
	  LiuParticles[j].INDVYDt = LiuParticles[j].INDVYDt - LiuParticles[i].Mass*h;
	}
	if (d==2) { 
	  LiuParticles[i].INDVZDt = LiuParticles[i].INDVZDt + LiuParticles[j].Mass*h;
	  LiuParticles[j].INDVZDt = LiuParticles[j].INDVZDt - LiuParticles[i].Mass*h;
	}
      }

      he = he*rhoij;
      LiuParticles[i].INDUDt = LiuParticles[i].INDUDt + LiuParticles[j].Mass*he;
      LiuParticles[j].INDUDt = LiuParticles[j].INDUDt + LiuParticles[i].Mass*he;
    }
//   For SPH algorithm 2 
           
    else if (LiuParam.pa_sph==2){
      for (d=0;d<Dim;d++) {
	h = -(LiuParticles[i].Press/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Press/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(d);

	if (d==1) he = he + (LiuParticles[j].XVel[0] - LiuParticles[i].XVel[0])*h;
	if (d==2) he = he + (LiuParticles[j].YVel[0] - LiuParticles[i].YVel[0])*h;
	if (d==3) he = he + (LiuParticles[j].ZVel[0] - LiuParticles[i].ZVel[0])*h;

 //   Viscous force 
	if (LiuParam.visc) {
	  if (d==0) {
//   x-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txx(i)/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Eta*txx(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(1-1);
	    if (Dim>=2){
	      h = h + (LiuParticles[i].Eta*txy(i)/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Eta*txy(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(2-1);
	      if (Dim==3) {
		h = h + (LiuParticles[i].Eta*txz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Eta*txz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(3-1);
	      }
	    }
	  }             
	  else if (d==1) {
//   y-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txy(i)/LiuParticles[i].Rho/LiuParticles[i].Rho
                  + LiuParticles[j].Eta*txy(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(1-1)
	          + (LiuParticles[i].Eta*tyy(i)/LiuParticles[i].Rho/LiuParticles[i].Rho
                  + LiuParticles[j].Eta*tyy(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(2-1);
	    if (Dim==3) {
	      h = h + (LiuParticles[i].Eta*tyz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho+LiuParticles[j].Eta*tyz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(3-1);
	    }
	  }
	  else if (d==2) {
//   z-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho +LiuParticles[j].Eta*txz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(1-1) 
	          + (LiuParticles[i].Eta*tyz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho +LiuParticles[j].Eta*tyz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(2-1)
	          + (LiuParticles[i].Eta*tzz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho +LiuParticles[j].Eta*tzz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*LiuPairs[k].DWDX(3-1);
	  }
	}
	//        MkDebug("Int_Force i:%d, j:%d, h:%f, Mass[i]:%f,Mass[i]*h:%f\n",i,j,h,LiuParticles[i].Mass,LiuParticles[i].Mass*h);
	if(d==0) {
	  LiuParticles[i].INDVXDt = LiuParticles[i].INDVXDt + LiuParticles[j].Mass*h;
	  LiuParticles[j].INDVXDt = LiuParticles[j].INDVXDt - LiuParticles[i].Mass*h;
	}
	if(d==1) {
	  LiuParticles[i].INDVYDt = LiuParticles[i].INDVYDt + LiuParticles[j].Mass*h;
	  LiuParticles[j].INDVYDt = LiuParticles[j].INDVYDt - LiuParticles[i].Mass*h;
	}
	if(d==2) {
	  LiuParticles[i].INDVZDt = LiuParticles[i].INDVZDt + LiuParticles[j].Mass*h;
	  LiuParticles[j].INDVZDt = LiuParticles[j].INDVZDt - LiuParticles[i].Mass*h;
	}
      }
      LiuParticles[i].INDUDt = LiuParticles[i].INDUDt + LiuParticles[j].Mass*he;
      LiuParticles[j].INDUDt = LiuParticles[j].INDUDt + LiuParticles[i].Mass*he;
    }
  }
 
//   Change of specific internal energy de/dt = T ds/dt - p/rho vc,c: 
 
  for (i=0;i<NTotal+NVirt;i++) {
    LiuParticles[i].INDUDt = LiuParticles[i].TDSDt + 0.5e0*LiuParticles[i].INDUDt;
    //    if(fabs(LiuParticles[i].INDVXDt)>0.000001) MkDebug("i:%d, INDVXDt:%12.6f\n",i,LiuParticles[i].INDVXDt);
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

   i= j= k= d =0;
   rr= f= rr0= dd= p1= p2=0;
            
  for (i = 0;i< NTotal+NVirt;i++) { 
    LiuParticles[i].EXDVXDt = 0.;
    LiuParticles[i].EXDVYDt = 0.;
    LiuParticles[i].EXDVZDt = 0.;
  }
         
//  Consider self-gravity or not ? 
 
  if (LiuParam.self_gravity) { 
    for (i = 0; i<NTotal+NVirt;i++) { 
      if (Dim == 1) LiuParticles[i].EXDVXDt = -9.8;
      if (Dim == 2) LiuParticles[i].EXDVYDt = -9.8;
      if (Dim == 3) LiuParticles[i].EXDVZDt = -9.8;
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
    if(LiuParticles[i].ParticleType>0&&LiuParticles[j].ParticleType<0) {   
      rr = 0.;       

      dx(0) =  LiuParticles[i].X - LiuParticles[j].X; 
      rr = rr + dx(0)*dx(0); 
      dx(1) =  LiuParticles[i].Y - LiuParticles[j].Y; 
      rr = rr + dx(1)*dx(1); 
      dx(2) =  LiuParticles[i].Z - LiuParticles[j].Z; 
      rr = rr + dx(2)*dx(2); 

      rr = sqrt(rr); 
      if(rr<rr0) { 
	f = (pow(rr0/rr,p1)-pow(rr0/rr,p2))/rr*rr ;

	LiuParticles[i].EXDVXDt = LiuParticles[i].EXDVXDt + dd*dx(0)*f;
	LiuParticles[i].EXDVYDt = LiuParticles[i].EXDVYDt + dd*dx(1)*f;
	LiuParticles[i].EXDVZDt = LiuParticles[i].EXDVZDt + dd*dx(2)*f;
      } 
    }         
  }   
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

   i=j=k=d=0;        
   vcc= epsilon=0; 
       
//  epsilon --- a small constants chosen by experience, may lead to instability. 
//  for example, for the 1 dimensional shock tube problem, the E = 0.3 
 
  epsilon = 0.3;
       
  for (i = 0;i< NTotal;i++ ) { 
    LiuParticles[i].XAVel[0] = 0.0;
    LiuParticles[i].YAVel[0] = 0.0;
    LiuParticles[i].ZAVel[0] = 0.0;
  }
      
  for (k=0;k<NIac;k++) {        
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J;

    dvx(0) = LiuParticles[i].XVel[0] - LiuParticles[i].XVel[0];
    LiuParticles[i].XAVel[0] = LiuParticles[i].XAVel[0] - 2*LiuParticles[j].Mass*dvx(0)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*LiuPairs[k].W;
    LiuParticles[j].XAVel[0] = LiuParticles[j].XAVel[0] + 2*LiuParticles[i].Mass*dvx(0)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*LiuPairs[k].W;
    dvx(1) = LiuParticles[i].YVel[0] - LiuParticles[i].YVel[0];
    LiuParticles[i].YAVel[0] = LiuParticles[i].YAVel[0] - 2*LiuParticles[j].Mass*dvx(1)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*LiuPairs[k].W;
    LiuParticles[j].YAVel[0] = LiuParticles[j].YAVel[0] + 2*LiuParticles[i].Mass*dvx(1)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*LiuPairs[k].W;
    dvx(2) = LiuParticles[i].ZVel[0] - LiuParticles[i].ZVel[0];
    LiuParticles[i].ZAVel[0] = LiuParticles[i].ZAVel[0] - 2*LiuParticles[j].Mass*dvx(2)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*LiuPairs[k].W;
    LiuParticles[j].ZAVel[0] = LiuParticles[j].ZAVel[0] + 2*LiuParticles[i].Mass*dvx(2)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*LiuPairs[k].W;

  }         
  for (i = 0;i< NTotal;i++) {
    LiuParticles[i].XAVel[0] = epsilon * LiuParticles[i].XAVel[0];
    LiuParticles[i].YAVel[0] = epsilon * LiuParticles[i].YAVel[0];
    LiuParticles[i].ZAVel[0] = epsilon * LiuParticles[i].ZAVel[0];
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

      dvx(0) = LiuParticles[j].XVel[0] - LiuParticles[i].XVel[0];
      dvx(1) = LiuParticles[j].YVel[0] - LiuParticles[i].YVel[0];
      dvx(2) = LiuParticles[j].ZVel[0] - LiuParticles[i].ZVel[0];

      hvcc = dvx(0)*LiuPairs[k].DWDX(1-1);
      for (d=1;d<Dim;d++) {
	hvcc = hvcc + dvx(d)*LiuPairs[k].DWDX(d-1);
      }
      vcc(i) = vcc(i) + LiuParticles[j].Mass*hvcc/LiuParticles[j].Rho; 
      vcc(j) = vcc(j) + LiuParticles[i].Mass*hvcc/LiuParticles[i].Rho;
    }
         
    for (i = 0;i< NTotal;i++) {
      dhsml(i) = (LiuParticles[i].SmoothLen/Dim)*vcc(i);
      LiuParticles[i].SmoothLen = LiuParticles[i].SmoothLen + Dt*dhsml(i);
      if (LiuParticles[i].SmoothLen<=0) LiuParticles[i].SmoothLen = LiuParticles[i].SmoothLen - Dt*dhsml(i);  
    }
  }
     
  else if(LiuParam.sle==1) {
    fac = 2.0;
    for (i = 0;i< NTotal;i++) {
      LiuParticles[i].SmoothLen = fac * (LiuParticles[i].Mass/pow(LiuParticles[i].Rho,(1./Dim)));
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
  if (LiuParam.shearcavity) Shear_Cavity();
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
 
  NTotal=400; 
  space_x=0.6/80.;       

  //  MkDebug("  shock_tube() : \n");

  for (i=0;i<NTotal;i++) {
    LiuParticles[i].Mass=0.75/400.; 
    LiuParticles[i].SmoothLen=0.015; 
    LiuParticles[i].ParticleType=1;

    LiuParticles[i].X = LiuParticles[i].Y = LiuParticles[i].Z = 0;
    LiuParticles[i].XVel[0] = LiuParticles[i].YVel[0] = LiuParticles[i].ZVel[0] = 0;
  }

  for (i=0;i<320;i++) {
    LiuParticles[i].X=-0.6+space_x/4.*(i-1); 
  }
       
  for (i=320;i<NTotal;i++) {
    LiuParticles[i].X=0.+space_x*(i-320);
  }         

  for (i=0;i<NTotal;i++) {
    if (LiuParticles[i].X<=1.e-8)  {
	LiuParticles[i].Energy=2.5; 
	LiuParticles[i].Rho=1. ;
	LiuParticles[i].Press=1. ;
      }
    if (LiuParticles[i].X>1.e-8){
      LiuParticles[i].Energy=1.795; 
      LiuParticles[i].Rho=0.25 ;
      LiuParticles[i].Press=0.1795 ;
    }   
  }
  //  MkDebug("  ~shock_tube(): \n ");
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
  double xl, yl, dx, dy; 
 
//   Giving data.Mass and smoothing length as well as other data. 
 
//MkDebug("shear_cavity()\n");
  m = 41; 
  n = 41; 
  mp = m-1; 
  np = n-1; 
  NTotal = mp * np; 
  xl = 1.e-3; 
  yl = 1.e-3; 
  dx = xl/mp; 
  dy = yl/np; 
 
  for (i = 0;i < mp;i++) {
    for (j = 0;j < np;j++) {
      k = j + i*np; 
      LiuParticles[k].X = (i)*dx + dx/2.0;
      LiuParticles[k].Y = (j)*dy + dy/2.0;
    }
  }
 
  for (i = 0;i < mp*np ;i++) {
    LiuParticles[i].XVel[0] = 0.; 
    LiuParticles[i].YVel[0] = 0.;       
    LiuParticles[i].Rho  = 1000.;    
    LiuParticles[i].Mass = dx*dy*LiuParticles[i].Rho;   
    LiuParticles[i].Press= 0.;    
    LiuParticles[i].Energy=357.1; 
    LiuParticles[i].ParticleType = 2; 
    LiuParticles[i].SmoothLen = dx;  
  }
//MkDebug("~shear_cavity()\n");
}

void MkLiuSPH::Shake()
{
  for (int i = 0;i < NTotal ;i++) {
    LiuParticles[i].X += 5000*1e-8 - (rand()%10000)*1e-8; 
    LiuParticles[i].Y += 5000*1e-8 - (rand()%10000)*1e-8; 
  }
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

   i= j= k= d=0;
   selfdens= r=0;
 
  MkDebug("\n      sum_density()\n");
//    wi(data.MaxN)---integration of the kernel itself 
   MkDebug("      SD1 ");         
  for (d=0;d<Dim;d++) { 
    hv(d) = 0.e0; 
  }
  MkDebug("SD2 ");          
//     Self density of each particle: Wii (Kernel for distance 0) 
//     and take contribution of particle itself: 
 
  r=0.;
       
//     Firstly calculate the integration of the kernel over the space 
 
  for (i=0;i<NTotal+NVirt ;i++) {
    LiuKernel.SetSmoothLen(LiuParticles[i].GetSmoothLen());
    selfdens = LiuKernel.W(r);
    wi(i)=selfdens*LiuParticles[i].Mass/LiuParticles[i].Rho; 
  }
  MkDebug("SD3 ");           
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    wi(i) = wi(i) + LiuParticles[j].Mass/LiuParticles[j].Rho*LiuPairs[k].W; 
    wi(j) = wi(j) + LiuParticles[i].Mass/LiuParticles[i].Rho*LiuPairs[k].W; 
  }
 
//  Secondly calculate the rho integration over the space 
  MkDebug("SD4 ");           
  for (i=0;i<NTotal+NVirt;i++) { 
    LiuKernel.SetSmoothLen(LiuParticles[i].GetSmoothLen());
    selfdens = LiuKernel.W(r);
    LiuParticles[i].Rho = selfdens*LiuParticles[i].Mass; 
  }
  MkDebug("SD5 \n");           
//  Calculate SPH sum for rho: 
  for (k=0;k<NIac;k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    LiuParticles[i].Rho = LiuParticles[i].Rho + LiuParticles[j].Mass*LiuPairs[k].W; 
    LiuParticles[j].Rho = LiuParticles[j].Rho + LiuParticles[i].Mass*LiuPairs[k].W; 
  }
  for (k=0 ;k<(int) min(10,NIac);k++) { 
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 
    MkDebug("k-%d, Roh(%d):%f,Roh(%d):%f,W:%f \n",k,i,LiuParticles[i].Rho,j,LiuParticles[j].Rho,LiuPairs[k].W);
  }
  MkDebug("SD6 \n");           

//  Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w) 
      
  if (LiuParam.nor_density) {  

    for (i=0;i< NTotal+NVirt;i++) { 
      LiuParticles[i].Rho=LiuParticles[i].Rho/wi(i);
    }

    for (k=NIac-10;k<NIac;k++) { 
      i = LiuPairs[k].I; 
      j = LiuPairs[k].J; 
      MkDebug("k-%d, Roh(%d):%f,Roh(%d):%f,wi:%f,wj:%f \n",k,i,LiuParticles[i].Rho,j,LiuParticles[j].Rho,wi(i), wi(j));
  }


  }  
  MkDebug("\n      ~sum_density()\n");       
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

   i=j=k=d=0;
   vcc=0;
       
  for (i = 0;i< NTotal+NVirt;i++) { 
    LiuParticles[i].DRhoDt = 0.; 
  }
      
  for (k=0;k<NIac;k++) {       
    i = LiuPairs[k].I; 
    j = LiuPairs[k].J; 

    dvx(0) = LiuParticles[i].XVel[0] - LiuParticles[j].XVel[0];
    dvx(1) = LiuParticles[i].YVel[0] - LiuParticles[j].YVel[0];
    dvx(2) = LiuParticles[i].ZVel[0] - LiuParticles[j].ZVel[0];

    vcc = dvx(0)*LiuPairs[k].DWDX(1-1);
    for (d=1;d<Dim;d++) { 
      vcc = vcc + dvx(d)*LiuPairs[k].DWDX(d-1); 
    }
    LiuParticles[i].DRhoDt = LiuParticles[i].DRhoDt + LiuParticles[j].Mass*vcc; 
    LiuParticles[j].DRhoDt = LiuParticles[j].DRhoDt + LiuParticles[i].Mass*vcc;
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

  if (LiuParam.skf==1) scale_k = 2;  
  else scale_k = 3; 

  for (i=0;i<NTotal+NVirt-1;i++) {      
    for (j = i+1;j< NTotal+NVirt;j++){ 
      float Xi,Xj,Yi,Yj;
      Xi = LiuParticles[i].X;
      Xj = LiuParticles[j].X;
      Yi = LiuParticles[i].Y;
      Yj = LiuParticles[j].Y;

      dx = LiuParticles[i].X - LiuParticles[j].X; 
      dy = LiuParticles[i].Y - LiuParticles[j].Y; 
      dz = LiuParticles[i].Z - LiuParticles[j].Z; 
      dist = sqrt(dx*dx+dy*dy+dz*dz);
      mhsml = (LiuParticles[i].SmoothLen+LiuParticles[j].SmoothLen)/2.; 
      if (dist < scale_k*mhsml) {
	MkDebug("%d (%f,%f), %d(%f,%f) \n",i,Xi,Yi,j,Xj,Yj);
	niac = niac + 1;  
      }
    }
  }
  MkDebug("Number of pair count %d\n",niac);getch();
  return niac;
}

void MkLiuSPH::Direct_Find() 
{ 
  int i, j, d,  sumiac, maxiac, miniac, niac, noiac, maxp, minp, scale_k;  
  double dx,dy,dz, dist, mhsml;      

  i= j= d=  sumiac= maxiac= miniac= niac = noiac= maxp= minp= scale_k=0;  
  dist= mhsml=0;      

//    Smoothing kernel function  
//    skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//        = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//        = 3, Quintic kernel (Morris 1997) 

  MkDebug("    direct_find()\n");
  MkDebug("    D1 ");

  if (LiuParam.skf==1) scale_k = 2;
  else scale_k = 3 ; 

  NIac = Pair_Count();
  LiuPairs.Initialize(NIac+1);
      
  for (i=0;i<NTotal+NVirt;i++) LiuParticles[i].CountIac = 0; 
       
  for (i=0;i<NTotal+NVirt-1;i++) {      
    for (j = i+1;j< NTotal+NVirt;j++){ 

      dx = LiuParticles[i].X - LiuParticles[j].X; 
      dy = LiuParticles[i].Y - LiuParticles[j].Y; 
      dz = LiuParticles[i].Z - LiuParticles[j].Z; 
      dist = sqrt(dx*dx+dy*dy+dz*dz);
      mhsml = (LiuParticles[i].SmoothLen+LiuParticles[j].SmoothLen)/2.; 
      LiuKernel.SetSmoothLen(mhsml);
      if (dist < scale_k*mhsml) { 
	niac = niac + 1;

	LiuParticles[i].CountIac = LiuParticles[i].CountIac + 1; 
	LiuParticles[j].CountIac = LiuParticles[j].CountIac + 1; 

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
      }
    }
  }
  MkDebug("D3 "); 
//    Statistics for the interaction 
 
  sumiac = 0; 
  maxiac = 0; 
  miniac = 1000; 
  noiac  = 0; 
  for (i=0;i<NTotal+NVirt;i++) { 
    sumiac = sumiac + LiuParticles[i].CountIac; 
    if (LiuParticles[i].CountIac>maxiac){ 
      maxiac = LiuParticles[i].CountIac; 
      maxp = i; 
    }
    if (LiuParticles[i].CountIac<miniac) {
      miniac = LiuParticles[i].CountIac; 
      minp = i; 
    }
    if (LiuParticles[i].CountIac==0) noiac  = noiac + 1;
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
  MkDebug("D4 "); 
  MkDebug("    ~direct_find()\n");
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
 
void MkLiuSPH::Virt_Part()  
{
  int i, j, d, im, mp,np,nvirt; 
  double xl, dx, v_inf; 
  FILE *fp1,*fp2, *fp3;

   i= j= d= im= mp=0; 
   xl= dx= v_inf=0; 
 
//MkDebug("    virt_part()\n");

  if (LiuParam.vp_input){
    //    LoadVPFromFile();
  }
  else  {
    nvirt = 0; 
    mp = 40;
    np = 40;
    xl = 1.0e-3; 
    dx = xl / mp; 
    v_inf = 1.e-3; 

//   Monaghan type virtual particle on the Upper side 
  //MkDebug("    virt_part()VP1 nvirt is %d ",nvirt); 
    for (i = 0;i< 2*mp+1;i++) {
      LiuParticles[NTotal + nvirt].X = i*dx/2;
      LiuParticles[NTotal + nvirt].Y = xl;   
      LiuParticles[NTotal + nvirt].XVel[0] = v_inf; 
      LiuParticles[NTotal + nvirt].YVel[0] = 0.; 
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_part()VP2 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Lower side 
 
    for (i = 0;i< 2*mp+1;i++) {
      LiuParticles[NTotal + nvirt].X = i*dx/2;  
      LiuParticles[NTotal + nvirt].Y = 0.;   
      LiuParticles[NTotal + nvirt].XVel[0] = 0.; 
      LiuParticles[NTotal + nvirt].YVel[0] = 0.; 
      nvirt = nvirt + 1; 
    }
  //MkDebug("    virt_part()VP3 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Left side 
 
    for (i = 0;i< 2*np-1;i++) {
      LiuParticles[NTotal + nvirt].X = 0.;
      LiuParticles[NTotal + nvirt].Y = (i+1)*dx/2;
      LiuParticles[NTotal + nvirt].XVel[0] = 0.; 
      LiuParticles[NTotal + nvirt].YVel[0] = 0.;
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_part()VP4 nvirt is %d ",nvirt); 
//   Monaghan type virtual particle on the Right side 
 
    for (i = 0;i< 2*np-1;i++) {
      LiuParticles[NTotal + nvirt].X = xl;
      LiuParticles[NTotal + nvirt].Y = (i+1)*dx/2;
      LiuParticles[NTotal + nvirt].XVel[0] = 0.; 
      LiuParticles[NTotal + nvirt].YVel[0] = 0.; 
      nvirt = nvirt + 1;
    }
  //MkDebug("    virt_part()VP5 nvirt is %d ",nvirt); 
    for (i = 0;i< nvirt;i++) {
      LiuParticles[NTotal + i].Rho = 1000.;
      LiuParticles[NTotal + i].Mass = LiuParticles[NTotal + i].Rho  * dx * dx; 
      LiuParticles[NTotal + i].Press = 0.; 
      LiuParticles[NTotal + i].Energy = 357.1; 
      LiuParticles[NTotal + i].ParticleType = -2; 
      LiuParticles[NTotal + i].SmoothLen = dx; 
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
  for (i=0;i<NTotal+NVirt;i++) {
    if (LiuParticles[i].GetParticleType()==1) LiuParticles[i].SetEta(0.0);
    else if (LiuParticles[i].GetParticleType()==2) LiuParticles[i].SetEta(1.0e-3);
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
//   gamma=7. 
//   rho0=1000.        
//   b = 1.013e5 
//   p = b*((rho/rho0)**gamma-1)       
//   c = 1480. 
 
//  Artificial EOS, Form 2 (Morris, 1997) 
//  See Equ.(4.89) 
  double p,c = 0.01; 
  p = c*c * rho;
  return p;
}

double MkLiuSPH::C_art_water(double rho) 
{       
//  Artificial EOS, Form 1 (Monaghan, 1994)  
//  See Equ.(4.88) 
//   gamma=7. 
//   rho0=1000.        
//   b = 1.013e5 
//   p = b*((rho/rho0)**gamma-1)       
//   c = 1480. 
 
//  Artificial EOS, Form 2 (Morris, 1997) 
//  See Equ.(4.89) 
  double c = 0.01; 
  return c;
}

 
