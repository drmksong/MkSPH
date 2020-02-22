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
  LiuGrid.Initialize(gx,gy,gz);
  LiuKernel.SetupKernelFunc(knlCubicSpline,/*Dim*/ 3,/*SmoothLen*/ 0.3);// 
}

void MkLiuSPH::Clear()
{
  Dim = 0;
  NVirt=0;
  NTotal=0;
      
  MaxInteraction = 0;
  NIac=0;        

  MaxTimeStep=0; 
  TimeStep=0;   
  Dt=0;       

  Pair_I.Clear();    
  Pair_J.Clear();    

  LiuParticles.Clear();
  LiuGrid.Clear;
  LiuKernel.Clear();
}
/*
bool MkLiuSPH::LoadFromFile()
{
  FILE *fp1,*fp2, *fp3;
  int i, d, im;

  if (FileName.length==0) {
    MkDebug("Please set the filename  \n");
    return false;
  }

  fp1 = fopen(xvFileName,"r");
  fp2 = fopen(stateFileName,"r");
  fp3 = fopen(otherFileName,"r");
       
  printf("  **************************************************\n"); 
  printf("      Loading initial particle configuration...   \n");
  fscanf (fp1,"%d", &NTotal);
  printf("      Total number of particles  %f\n ", NTotal)    	 ;
  printf("  **************************************************\n")	 ;
  for (i = 1;i<= NTotal;i++) {
    fscanf(fp1,"%d",&im);
    fscanf(fp1,"%f %f %f",&LiuParticles[i].X,&LiuParticles[i].Y,&LiuParticles[i].Z);
    for (d=1;d<=Dim;d++) fscanf(fp1,"%f",&data.VX(d, i));
    fscanf(fp2,"%d %f %f %f %f",&im, &LiuParticles[i].Mass, &LiuParticles[i].Rho, &LiuParticles[i].Press, &LiuParticles[i].Energy);         
    fscanf(fp3,"%d %d %f",&im, &LiuParticles[i].ParticleType, &LiuParticles[i].SmoothLen);
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
}

bool MkLiuSPH::LoadVPFromFile()
{

  int i, j, d, im, mp; 
  double xl, dx, v_inf; 
  FILE *fp1,*fp2, *fp3;

  if (FileName.length==0) {
    MkDebug("Please set the filename  \n");
    return false;
  }

  i= j= d= im= mp=0; 
  xl= dx= v_inf=0; 
                      
  fp1 = fopen(vp_xvFileName,"r"); 
  fp2 = fopen(vp_stateFileName,"r"); 
  fp3 = fopen(vp_otherFileName,"r");
  fscanf(fp1,"%d",&NVirt); 
  for (j = 1;j<= NVirt;j++) {
    i = NTotal + j;
    fscanf(fp1,"%d",&im);
    fscanf(fp1,"%f %f %f",&LiuParticles[i].X,&LiuParticles[i].Y,&LiuParticles[i].Z);
    for (d=1;d<=Dim;d++) fscanf(fp1,"%f",&data.VX(d, i));
    fscanf(fp2,"%d %f %f %f %f",im, LiuParticles[i].Mass, LiuParticles[i].Rho, LiuParticles[i].Press, LiuParticles[i].Energy);
    fscanf(fp3,"%d %d %f",im, LiuParticles[i].ParticleType, LiuParticles[i].SmoothLen);
  }
  fclose(fp1); 
  fclose(fp2); 
  fclose(fp3);
}

bool MkLiuSPH::SaveToFile()
{
  int i, d;
  FILE *fp1,*fp2, *fp3;
       
  fp1 = fopen(xvFileName,"w"); 
  fp2 = fopen(stateFileName,"w");
  fp3 = fopen(otherFileName,"w");
      
  fprintf(fp1,"%d\n", NTotal); 
  for (i = 1;i<= NTotal;i++) {
    fprintf(fp1,"%6d ",i);
    fprintf(fp1,"%14.8f %14.8f %14.8f ", LiuParticles[i].X, LiuParticles[i].Y, LiuParticles[i].Z);
    for(d=1;d<=Dim;d++) fprintf(fp1,"%14.8f ",(data.VX(d, i)));
    fprintf(fp1,"\n ");
    fprintf(fp2,"%6d %14.8f %14.8f %14.8f %14.8f \n", i, LiuParticles[i].Mass, LiuParticles[i].Rho, LiuParticles[i].Press, LiuParticles[i].Energy); 
    fprintf(fp3,"%6d %4d %14.8f", i, LiuParticles[i].ParticleType, LiuParticles[i].SmoothLen);                                
  }
       
  fclose(fp1); 
  fclose(fp2);
  fclose(fp3);

}

bool MkLiuSPH::SaveVPToFile()
{
  int i, d;
  FILE *fp1,*fp2, *fp3;

  fp1 = fopen(vp_xvFileName,"w"); 
  fp2 = fopen(vp_stateFileName,"w"); 
  fp3 = fopen(vp_otherFileName,"w");

  fprintf(fp1,"\n", NVirt) ;
  for (i = NTotal + 1;i<= NTotal + NVirt;i++) {
    fprintf(fp1,"%6d ",i);
    fprintf(fp1,"%14.8f %14.8f %14.8f ", LiuParticles[i].X, LiuParticles[i].Y, LiuParticles[i].Z);
    for (d=1;d<=Dim;d++) fprintf(fp1,"%14.8f ", data.VX(d, i));
    fprintf(fp1,"\n");
    fprintf(fp2,"%6d %14.8f %14.8f %14.8f %14.8f\n", i, LiuParticles[i].Mass, LiuParticles[i].Rho, LiuParticles[i].Press, LiuParticles[i].Energy);
    fprintf(fp3,"%6d %4d %14.8f %14.8f\n", i, LiuParticles[i].ParticleType, LiuParticles[i].SmoothLen);                         
  }
  fclose(fp1); 
  fclose(fp2); 
  fclose(fp3);
}
*/


void MkLiuSPH::Draw()
{
  LiuParticles.Draw();
  LiuGrid.Draw();
} 

void MkLiuSPH::Run()
{
  
}

//---------------------------------------------------------------------- 
//    X-- coordinates of particles                       [input/output] 
//    data.VX-- velocities of particles                       [input/output] 
//    data.Mass-- mass of particles                                  [input] 
//    data.Rho-- dnesities of particles                       [input/output] 
//    data.P-- pressure  of particles                         [input/output] 
//    data.U-- internal energy of particles                   [input/output] 
//    data.C-- sound velocity of particles                          [output] 
//    data.IType-- types of particles                               [input] 
//         =1   ideal gas 
//         =2   water 
//         =3   tnt 
//    data.Hsml-- smoothing lengths of particles              [input/output] 
//    NTotal-- total particle number                            [input]   
//    MaxTimeStep-- maximum timesteps                           [input] 
//    Dt-- timestep                                             [input] 
    
void MkLiuSPH::Time_Integration(SPH_Data &data) 
{
  int i, j, k, d, current_ts, nstart;
  double  time, temp_rho, temp_u;

  MkDouble  x_min(Dim+1, data.MaxN+1), v_min(Dim+1, data.MaxN+1), u_min(data.MaxN+1);
  MkDouble rho_min(data.MaxN+1); 

  //  MkDouble  dvx(Dim+1, data.MaxN+1);
  i= j= k= d= current_ts= nstart=0;
  time= temp_rho= temp_u=0;
           
  //  MkDebug("time_integration()\n");
     
  for (i = 1;i<= NTotal;i++) {
    LiuParticles[i].XVel[0] = 0;
    LiuParticles[i].YVel[0] = 0;
    LiuParticles[i].ZVel[0] = 0;
  }

  for (CurrentTimeStep = nstart+1;CurrentTimeStep<= nstart+MaxTimeStep;CurrentTimeStep++) {
    current_ts=current_ts+1;
    if ((CurrentTimeStep%print_step)==0) {
      printf("______________________________________________\n");
      printf("  current number of time step = %d     current time=%f\n",CurrentTimeStep, float(time+Dt)) ;
      printf("______________________________________________\n");
    }
       
//   If not first time step, then update thermal energy, density and  
//   velocity half a time step   

//    MkDebug("  T1 ");
 
    if (CurrentTimeStep != 1) {
 
      for (i = 1;i<= NTotal;i++) {
	u_min(i) = LiuParticles[i].Energy;
	temp_u=0.;
	if (Dim==1) temp_u=-nsym*LiuParticles[i].Press*data.VX(1,i)/LiuParticles[i].X/LiuParticles[i].Rho;
	LiuParticles[i].Energy = LiuParticles[i].Energy + (Dt/2.)* (data.DUDt(i)+temp_u);
	if(LiuParticles[i].Energy<0)  LiuParticles[i].Energy = 0.;
             
	if (!summation_density) {
	  rho_min(i) = LiuParticles[i].Rho;
	  temp_rho=0.; 
	  if (Dim==1) temp_rho=-nsym*LiuParticles[i].Rho*data.VX(1,i)/LiuParticles[i].X;
	  LiuParticles[i].Rho = LiuParticles[i].Rho +(Dt/2.)*( data.DRhoDt(i)+ temp_rho);
	}

	for (d = 1;d<= Dim;d++) {
	  v_min(d, i) = data.VX(d, i) ;
	  data.VX(d, i) = data.VX(d, i) + (Dt/2.)*data.DVXDt(d, i) ;
	}
      }
    }
//---  Definition of variables out of the function vector:     

//    MkDebug("  T2 ");       

    single_step(data);
                 
    //    MkDebug("  T3 ");
  
    if (CurrentTimeStep == 1) { 
        
      for (i=1;i<=NTotal;i++) {
	temp_u=0.;
	if (Dim==1) temp_u=-nsym*LiuParticles[i].Press*data.VX(1,i)/LiuParticles[i].X/LiuParticles[i].Rho;
	LiuParticles[i].Energy = LiuParticles[i].Energy + (Dt/2.)*(data.DUDt(i) + temp_u); 
	if(LiuParticles[i].Energy<0)  LiuParticles[i].Energy = 0.;
          
	if (!summation_density ){
	  temp_rho=0.; 
	  if (Dim==1) temp_rho=-nsym*LiuParticles[i].Rho*data.VX(1,i)/LiuParticles[i].X ;
	  LiuParticles[i].Rho = LiuParticles[i].Rho + (Dt/2.)* (data.DRhoDt(i)+temp_rho) ;
	}
          
	//    MkDebug("  T4 ");
	for (d = 1;d<= Dim;d++) {
	  data.VX(d, i) = data.VX(d, i) + (Dt/2.) * data.DVXDt(d, i) + data.AveVel(d, i); 
	  LiuParticles[i].X = LiuParticles[i].X + Dt * data.VX(d, i);
	}
      }
    }               
    else {
      //MkDebug("  T5 ");
                     
      for (i=1;i<=NTotal;i++) {
	temp_u=0.;
	if (Dim==1) temp_u=-nsym*LiuParticles[i].Press*data.VX(1,i)/LiuParticles[i].X/LiuParticles[i].Rho;
	LiuParticles[i].Energy = u_min(i) + Dt*(data.DUDt(i)+temp_u);
	if(LiuParticles[i].Energy<0)  LiuParticles[i].Energy = 0.;
             
	if (!summation_density ) {
	  temp_rho=0.;
	  if (Dim==1) temp_rho=-nsym*LiuParticles[i].Rho*data.VX(1,i)/LiuParticles[i].X;
	  LiuParticles[i].Rho = rho_min(i) + Dt*(data.DRhoDt(i)+temp_rho);
	}
                 
	for (d = 1;d<= Dim;d++) {
	  data.VX(d, i) = v_min(d, i) + Dt * data.DVXDt(d, i) + data.AveVel(d, i);
	  LiuParticles[i].X = LiuParticles[i].X + Dt * data.VX(d, i);
	}
      }
         
    }      
    time = time + Dt;
 
    if ((CurrentTimeStep%save_step)==0) {
      output(data);
    }

    //    MkDebug("  T10 ");
  
    if ((CurrentTimeStep%print_step)==0) {
      //      123456789abc 123456789abc 123456789abc 123456789abc 
      printf("\n");
      printf("x            velocity     dvx\n");
      printf("%12.6f %12.6f %12.6f \nf",LiuParticles[moni_particle].X(), data.VX(1,moni_particle),data.DVXDt(1,moni_particle));
    }
  }
  nstart=current_ts;
  //  MkDebug("~time_integration()\n");

}
 
//---------------------------------------------------------------------- 
// Subroutine to determine the right hand side of a differential  
// equation in a single step for performing time integration  
 
// In this routine and its subroutines the SPH algorithms are performed. 
//   CurrentTimeStep: Current timestep number                            [in] 
//   Dt       : Timestep                                           [in] 
//   NTotal   :  Number of particles                               [in] 
//   data.Hsml     :  Smoothing Length                                  [in] 
//   data.Mass     :  Particle masses                                   [in] 
//   data.X        :  Particle position                                 [in] 
//   data.VX       :  Particle velocity                                 [in] 
//   data.U        :  Particle internal energy                          [in] 
//   data.Rho      :  Density                                       [in/out] 
//   p        :  Pressure                                         [out] 
//   t        :  Temperature                                   [in/out] 
//   data.TDSDt    :  Production of viscous entropy t*ds/dt            [out] 
//   dx       :  dx = data.VX = dx/dt                                  [out] 
//   dvx      :  dvx = dvx/dt, force per unit mass                [out] 
//   data.DUDt       :  du  = du/dt                                      [out] 
//   data.DSDt       :  ds  = ds/dt                                      [out]      
//   data.DRhoDt   :  drhodt =  drho/dt                                  [out] 
//   data.IType    :  Type of particle                                 [in] 
//   av       :  Monaghan average velocity                        [out] 
 
void MkLiuSPH::Single_Step(SPH_Data &data)  
{
  int i, d;
  MkInt ns(data.MaxN+1);
  MkDouble indvxdt(Dim+1,data.MaxN+1);
  MkDouble exdvxdt(Dim+1,data.MaxN+1),ardvxdt(Dim+1,data.MaxN+1),avdudt(data.MaxN+1),ahdudt(data.MaxN+1);
 
  MkDebug("\n  single_step()\n");
  MkDebug("  S1 ");
  for ( i=1;i<=NTotal;i++) {
    avdudt(i) = 0.; 
    ahdudt(i) = 0.; 
    for ( d=1;d<=Dim;d++) {
      indvxdt(d,i) = 0.; 
      ardvxdt(d,i) = 0.; 
      exdvxdt(d,i) = 0.;
    }
  }
  MkDebug("  S2 ");  
//---  Positions of virtual (boundary) particles:  
 
  NVirt = 0; 
  if (virtual_part) {
    virt_part(data); 
  }
  MkDebug("  S3 ");      
//---  Interaction parameters, calculating neighboring particles 
//   and optimzing smoothing length 
 
  if (nnps==1) { direct_find(data) ; ns.CopyFrom(data.CountIac);  }
  else if (nnps==2) { link_list(data); ns.CopyFrom(data.CountIac);  }
//  else if (nnps==3) 
//    tree_search(CurrentTimeStep, NTotal+nvirt,data.Hsml,x,niac,pair_i,pair_j,w,dwdx,ns);
MkDebug("  S4 ");
//---  Density approximation or change rate 

  if (summation_density) sum_density(data);
  else con_density(data);
 
//---  Dynamic viscosity: 
 
MkDebug("  S5 ");
  if (visc) Viscosity(); 
        
//---  Internal forces: 
  
  int_force(data) ; indvxdt.CopyFrom(data.INDVXDt); data.DUDt.CopyFrom(data.INDUDt);
               
  MkDebug("  S6 ");    

//---  Artificial viscosity: 
 
  if (visc_artificial) {art_visc(data); ardvxdt.CopyFrom(data.ARDVXDt);avdudt.CopyFrom(data.AVDUDt);}
  MkDebug("  S7 ");
       
//---  External forces: 
 
MkDebug("  S8 ");

  if (ex_force) {ext_force(data);exdvxdt.CopyFrom(data.EXDVXDt);}
 
  MkDebug("  S9 ");

//   Calculating the neighboring particles and undating DATA.HSML 
       
  if (sle!=0) h_upgrade(data); 
 
  if (heat_artificial) {art_heat(data);ahdudt.CopyFrom(data.AHDUDt);}
      
//   Calculating average velocity of each partile for avoiding penetration 
 
  if (average_velocity) av_vel(data);
 
//---  Convert velocity, force, and energy to f and dfdt   
 
  for (i=1;i<=NTotal;i++) {
    for (d=1;d<=Dim;d++) {
      data.DVXDt(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i);
    }
    data.DUDt(i) = data.DUDt(i) + avdudt(i) + ahdudt(i);
  }
  //   MkDebug("  S10 ");
  if ((CurrentTimeStep%print_step)==0) {
    //      123456789abc 123456789abc 123456789abc
    printf("\n") ;
    printf("**** Information for particle **** %d\n",moni_particle);
    printf("internal a   artifical a  external a    total a \n");
    printf("%12.6f %12.6f %12.6f %12.6f\n",indvxdt(1,moni_particle),ardvxdt(1,moni_particle),exdvxdt(1,moni_particle),data.DVXDt(1,moni_particle));
  }

  //  MkDebug("------------------ Data Dumping %d --------------------\n",CurrentTimeStep);
  //data.Dump();
  //MkDebug("------------------ ~Data Dumping --------------------\n");

  MkDebug("  ~single_step()\n");
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
  MkDouble vcc(NTotal+1),dvx(Dim+1);

  i=j=k=d =0;
  dx= vr= rr= h= mc= mrho= mhsml= hvcc= mui= muj= muij= rdwdx= g1=g2=0;
       
  //---  Parameter for the artificial heat conduction:
      
  g1=0.1; 
  g2=1.0;
  for (i=1;i<=NTotal+NVirt;i++) { 
    vcc(i) = 0.e0;
    data.AHDUDt(i) = 0.e0; 
  }
     
  for (k=1;k<=NIac;k++) { 
    i = Pair_I(k); 
    j = Pair_J(k); 
    for (d=1;d<=Dim;d++) { 
      dvx(d) = data.VX(d,j) - data.VX(d,i);  
    }        
    hvcc = dvx(1)*data.DWDX(1,k);
    for (d=2;d<=Dim;d++) { 
      hvcc = hvcc + dvx(d)*data.DWDX(d,k) ;
    }    
    vcc(i) = vcc(i) + LiuParticles[j].Mass*hvcc/LiuParticles[j].Rho; 
    vcc(j) = vcc(j) + LiuParticles[i].Mass*hvcc/LiuParticles[i].Rho; 
  }
    
  for (k=1;k<=NIac;k++) { 
    i = Pair_I(k); 
    j = Pair_J(k); 
    mhsml= (LiuParticles[i].SmoothLen+LiuParticles[j].SmoothLen)/2.; 
    mrho = 0.5e0*(LiuParticles[i].Rho + LiuParticles[j].Rho);
    rr = 0.e0; 
    rdwdx = 0.e0; 
/*    for (d=1;d<=Dim;d++) { 
      dx = data.X(d,i) -  data.X(d,j);
      rr = rr + dx*dx; 
      rdwdx  = rdwdx + dx*data.DWDX(d,k);
      }
*/  
//  not satisfactory
    dx = LiuParticles[i].X -  LiuParticles[j].X;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*data.DWDX(d,k);
    dx = LiuParticles[i].Y -  LiuParticles[j].Y;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*data.DWDX(d,k);
    dx = LiuParticles[i].Z -  LiuParticles[j].Z;
    rr = rr + dx*dx; 
    rdwdx  = rdwdx + dx*data.DWDX(d,k);


    mui=g1*LiuParticles[i].SmoothLen*data.C(i) + g2*LiuParticles[i].SmoothLen*LiuParticles[i].SmoothLen*(fabs(vcc(i))-vcc(i));
    muj=g1*LiuParticles[j].SmoothLen*data.C(j) + g2*LiuParticles[j].SmoothLen*LiuParticles[j].SmoothLen*(fabs(vcc(j))-vcc(j));
    muij= 0.5*(mui+muj);
    h = muij/(mrho*(rr+0.01*mhsml*mhsml))*rdwdx ;
    data.AHDUDt(i) = data.AHDUDt(i) + LiuParticles[j].Mass*h*(LiuParticles[i].Energy-LiuParticles[j].Energy) ;
    data.AHDUDt(j) = data.AHDUDt(j) + LiuParticles[i].Mass*h*(LiuParticles[j].Energy-LiuParticles[i].Energy); 
  }
 
  for (i=1;i<=NTotal+NVirt;i++) { 
    data.AHDUDt(i) = 2.0e0*data.AHDUDt(i);           
  }
}

 
//---------------------------------------------------------------------- 
//  Subroutine to calculate the artificial viscosity (Monaghan, 1992)  
//  See Equ.(4.66) Equ.(4.62) 
 
//  NTotal : Number of particles (including virtual particles)    [in] 
//  data.Hsml   : Smoothing Length                                     [in] 
//  data.Mass   : Particle masses                                      [in] 
//  data.X      : Coordinates of all particles                         [in] 
//  data.VX     : Velocities of all particles                          [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  data.Rho    : Density                                              [in] 
//  data.C      : Temperature                                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  data.W      : Kernel for all interaction pairs                     [in] 
//  data.DWDX   : Derivative of kernel( with respect to x, y and z      [in] 
//  data.ARDVXDt  : Acceleration with respect to x, y and z             [out]  
//  data.AVDUDt   : Change of specific internal energy                  [out] 

void MkLiuSPH::Art_Visc(SPH_Data &data)
{
  int i,j,k,d ;
  double dx, alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml;
  MkDouble dvx(3+1);

   i=j=k=d =0;
   dx= alpha= beta= etq= piv= muv= vr= rr= h= mc= mrho= mhsml=0;

//  Parameter for the artificial viscosity: 
//  Shear viscosity 
  alpha = 1.e0; 
      
//  Bulk viscosity 
  beta  = 1.e0;  
       
//  Parameter to avoid singularities 
  etq   = 0.1e0;
            
  for (i=1;i<=NTotal+NVirt;i++){ 
    for (d=1;d<=Dim;d++) { 
      data.ARDVXDt(d,i) = 0.e0;
    }
    data.AVDUDt(i) = 0.e0;
  }
      
//  Calculate SPH sum for artificial viscosity 
       
  for (k=1;k<=NIac;k++) { 
    i = Pair_I(k); 
    j = Pair_J(k); 
    mhsml= (LiuParticles[i].SmoothLen+LiuParticles[j].SmoothLen)/2.; 
    vr = 0.e0; 
    rr = 0.e0; 
    /*
    for (d=1;d<=Dim;d++) { 
      dvx(d) = data.VX(d,i) - data.VX(d,j); 
      dx     =  data.X(d,i) -  data.X(d,j); 
      vr     = vr + dvx(d)*dx; 
      rr     = rr + dx*dx; 
    }
    */
    dvx(1) = LiuParticles[i].XVel[0] - LiuParticles[j].XVel[0]; 
    dx     = LiuParticles[i].X - LiuParticles[j].X; 
    vr     = vr + dvx(1)*dx; 
    rr     = rr + dx*dx; 
    dvx(2) = LiuParticles[i].YVel[0] - LiuParticles[j].YVel[0]; 
    dx     = LiuParticles[i].Y - LiuParticles[j].Y; 
    vr     = vr + dvx(1)*dx; 
    rr     = rr + dx*dx; 
    dvx(2) = LiuParticles[i].ZVel[0] - LiuParticles[j].ZVel[0]; 
    dx     = LiuParticles[i].Z - LiuParticles[j].Z; 
    vr     = vr + dvx(1)*dx; 
    rr     = rr + dx*dx; 

//  Artificial viscous force only if v_ij * r_ij  0 
 
    if (vr<0.e0) { 
 
//  Calculate muv_ij = data.Hsml v_ij * r_ij / ( r_ij^2 + data.Hsml^2 etq^2 ) 
             
      muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq);
           
//  Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij 
 
      mc   = 0.5e0*(data.C(i) + data.C(j)); 
      mrho = 0.5e0*(LiuParticles[i].Rho + LiuParticles[j].Rho); 
      piv  = (beta*muv - alpha*mc)*muv/mrho;
 
//  Calculate SPH sum for artificial viscous force 
 
      for (d=1;d<=Dim;d++) { 
	h = -piv*data.DWDX(d,k); 
	data.ARDVXDt(d,i) = data.ARDVXDt(d,i) + LiuParticles[j].Mass*h; 
	data.ARDVXDt(d,j) = data.ARDVXDt(d,j) - LiuParticles[i].Mass*h; 
	data.AVDUDt(i) = data.AVDUDt(i) - LiuParticles[j].Mass*dvx(d)*h; 
	data.AVDUDt(j) = data.AVDUDt(j) - LiuParticles[i].Mass*dvx(d)*h; 
      }
    }
  }
 
//  Change of specific internal energy: 
 
  for (i=1;i<=NTotal+NVirt;i++) { 
    data.AVDUDt(i) = 0.5e0*data.AVDUDt(i);
  }
 
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
//    data.Hsml   : Smoothing Length                                     [in] 
//    data.Mass   : Particle masses                                      [in] 
//    data.VX     : Velocities of all particles                          [in] 
//    nia//  : Number of interaction pairs                          [in] 
//    data.Rho    : Density                                              [in] 
//    data.Eta    : Dynamic viscosity                                    [in] 
//    Pair_I : List of first partner of interaction pair            [in] 
//    Pair_J : List of second partner of interaction pair           [in] 
//    data.DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//    data.IType  : Type of particle (material types)                    [in] 
//    data.U      : Particle internal energy                             [in] 
//    data.X      : Particle coordinates                                 [in] 
//    t      : Particle temperature                             [in/out] 
//   data.C      : Particle sound speed                                [out] 
//   p      : Particle pressure                                   [out] 
//   data.INDVXDt  : Acceleration with respect to x, y and z             [out]  
//   data.TDSDt  : Production of viscous entropy                       [out] 
//   data.INDUDt   : Change of specific internal energy                  [out] 


void MkLiuSPH::Int_Force(SPH_Data &data) 
{
  int i, j, k, d;
  double hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij;
  MkDouble dvx(Dim+1), txx(data.MaxN+1), tyy(data.MaxN+1),tzz(data.MaxN+1), txy(data.MaxN+1), txz(data.MaxN+1), tyz(data.MaxN+1), vcc(data.MaxN+1);

//   Initialization of shear tensor, velocity divergence,  
//   viscous energy, internal energy, acceleration  

   i= j= k= d=0;
   hxx= hyy= hzz= hxy= hxz= hyz= h= hvcc= he= rhoij=0;

  for (i=1;i<=NTotal+NVirt;i++) {
    txx(i) = 0.e0; 
    tyy(i) = 0.e0; 
    tzz(i) = 0.e0; 
    txy(i) = 0.e0; 
    txz(i) = 0.e0; 
    tyz(i) = 0.e0; 
    vcc(i) = 0.e0; 
    data.TDSDt(i) = 0.e0; 
    data.INDUDt(i) = 0.e0; 
    for (d=1;d<=Dim;d++) { 
      data.INDVXDt(d,i) = 0.e0;
    }
  }
 
//   Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c 
 
  if (visc) {
    for (k=1;k<=NIac;k++) { 
      i = Pair_I(k); 
      j = Pair_J(k); 
      for (d=1;d<=Dim;d++) {
	dvx(d) = data.VX(d,j) - data.VX(d,i);
      }
      if (Dim==1) {
	hxx = 2.e0*dvx(1)*data.DWDX(1,k);
      }
      else if (Dim==2) {
	hxx = 2.e0*dvx(1)*data.DWDX(1,k) -  dvx(2)*data.DWDX(2,k);
	hxy = dvx(1)*data.DWDX(2,k) + dvx(2)*data.DWDX(1,k);
	hyy = 2.e0*dvx(2)*data.DWDX(2,k) - dvx(1)*data.DWDX(1,k);
      }
      else if (Dim==3){
	hxx = 2.e0*dvx(1)*data.DWDX(1,k) - dvx(2)*data.DWDX(2,k)- dvx(3)*data.DWDX(3,k);
	hxy = dvx(1)*data.DWDX(2,k) + dvx(2)*data.DWDX(1,k);
	hxz = dvx(1)*data.DWDX(3,k) + dvx(3)*data.DWDX(1,k);
	hyy = 2.e0*dvx(2)*data.DWDX(2,k) - dvx(1)*data.DWDX(1,k)- dvx(3)*data.DWDX(3,k);
	hyz = dvx(2)*data.DWDX(3,k) + dvx(3)*data.DWDX(2,k);
	hzz = 2.e0*dvx(3)*data.DWDX(3,k) - dvx(1)*data.DWDX(1,k) - dvx(2)*data.DWDX(2,k);
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
      for (d=1;d<=Dim;d++) {
	hvcc = hvcc + dvx(d)*data.DWDX(d,k);
      }
      vcc(i) = vcc(i) + LiuParticles[j].Mass*hvcc/LiuParticles[j].Rho;
      vcc(j) = vcc(j) + LiuParticles[i].Mass*hvcc/LiuParticles[i].Rho;
    }
  }
 
  for (i=1;i<=NTotal+NVirt;i++) {
//   Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab 
    if (visc) {
      if (Dim==1) {
	data.TDSDt(i) = txx(i)*txx(i);
      }
      else if (Dim==2) {
	data.TDSDt(i) = txx(i)*txx(i) + 2.e0*txy(i)*txy(i)+ tyy(i)*tyy(i);
      }
      else if (Dim==3) {
	data.TDSDt(i) = txx(i)*txx(i)+ tyy(i)*tyy(i)+ tzz(i)*tzz(i)
                      + 2.e0*txy(i)*txy(i)+ 2.e0*txz(i)*txz(i)+ 2.e0*tyz(i)*tyz(i)  ;
      }
      data.TDSDt(i) = 0.5e0*LiuParticles[i].Eta/LiuParticles[i].Rho*data.TDSDt(i);
    }
 
//   Pressure from equation of state 
 
    if (abs(LiuParticles[i].ParticleType)==1) {
      LiuParticles[i].Press = P_ideal_gas(LiuParticles[i].Rho, LiuParticles[i].Energy);
      data.C(i) = C_ideal_gas(LiuParticles[i].Rho, LiuParticles[i].Energy);
    }
    else if (abs(LiuParticles[i].ParticleType)==2) {
      LiuParticles[i].Press = P_art_water(LiuParticles[i].Rho);
      data.C(i) = C_art_water(LiuParticles[i].Rho);
    }
  } 
//    Calculate SPH sum for pressure force -p,a/rho 
//    and viscous force (eta Tab),b/rho 
//    and the internal energy change de/dt due to -p/rho vc,c 
 
  for (k=1;k<=NIac;k++) {
    i = Pair_I(k);
    j = Pair_J(k);
    he = 0.e0;
       
//   For SPH algorithm 1 
 
    rhoij = 1.e0/(LiuParticles[i].Rho*LiuParticles[j].Rho);
    if(pa_sph==1){
      for (d=1;d<=Dim ;d++) {
         
//   Pressure part 
                     
	h = -(LiuParticles[i].Press + LiuParticles[j].Press)*data.DWDX(d,k);
	he = he + (data.VX(d,j) - data.VX(d,i))*h;
 
//   Viscous force 
 
	if (visc) {
	  if (d==1) {
//   x-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txx(i) + LiuParticles[j].Eta*txx(j))*data.DWDX(1,k);
	    if (Dim>=2) {
	      h = h + (LiuParticles[i].Eta*txy(i) + LiuParticles[j].Eta*txy(j))*data.DWDX(2,k);
	      if (Dim==3) {
		h = h + (LiuParticles[i].Eta*txz(i) + LiuParticles[j].Eta*txz(j))*data.DWDX(3,k);
	      }
	    }
	  }
	  else if (d==2) {
//   y-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txy(i) + LiuParticles[j].Eta*txy(j))*data.DWDX(1,k) + (LiuParticles[i].Eta*tyy(i) + LiuParticles[j].Eta*tyy(j))*data.DWDX(2,k);
	    if (Dim==3) {
	      h = h + (LiuParticles[i].Eta*tyz(i) + LiuParticles[j].Eta*tyz(j))*data.DWDX(3,k);
	    }
	  }
	  else if (d==3) {
//   z-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txz(i) + LiuParticles[j].Eta*txz(j))*data.DWDX(1,k) + (LiuParticles[i].Eta*tyz(i) + LiuParticles[j].Eta*tyz(j))*data.DWDX(2,k)+ (LiuParticles[i].Eta*tzz(i) + LiuParticles[j].Eta*tzz(j))*data.DWDX(3,k);
	  }
	}
	h = h*rhoij;
	data.INDVXDt(d,i) = data.INDVXDt(d,i) + LiuParticles[j].Mass*h;
	data.INDVXDt(d,j) = data.INDVXDt(d,j) - LiuParticles[i].Mass*h;
      }
      he = he*rhoij;
      data.INDUDt(i) = data.INDUDt(i) + LiuParticles[j].Mass*he;
      data.INDUDt(j) = data.INDUDt(j) + LiuParticles[i].Mass*he;
    }
//   For SPH algorithm 2 
           
    else if (pa_sph==2){
      for (d=1;d<=Dim;d++) {
	h = -(LiuParticles[i].Press/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Press/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(d,k);
	he = he + (data.VX(d,j) - data.VX(d,i))*h;
 //   Viscous force 
	if (visc) {
	  if (d==1) {
//   x-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txx(i)/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Eta*txx(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(1,k);
	    if (Dim>=2){
	      h = h + (LiuParticles[i].Eta*txy(i)/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Eta*txy(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(2,k);
	      if (Dim==3) {
		h = h + (LiuParticles[i].Eta*txz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho + LiuParticles[j].Eta*txz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(3,k);
	      }
	    }
	  }             
	  else if (d==2) {
//   y-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txy(i)/LiuParticles[i].Rho/LiuParticles[i].Rho+LiuParticles[j].Eta*txy(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(1,k)
	          + (LiuParticles[i].Eta*tyy(i)/LiuParticles[i].Rho/LiuParticles[i].Rho+LiuParticles[j].Eta*tyy(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(2,k);
	    if (Dim==3) {
	      h = h + (LiuParticles[i].Eta*tyz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho+LiuParticles[j].Eta*tyz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(3,k);
	    }
	  }
	  else if (d==3) {
//   z-coordinate of acceleration 
	    h = h + (LiuParticles[i].Eta*txz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho +LiuParticles[j].Eta*txz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(1,k) 
	          + (LiuParticles[i].Eta*tyz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho +LiuParticles[j].Eta*tyz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(2,k)
	          + (LiuParticles[i].Eta*tzz(i)/LiuParticles[i].Rho/LiuParticles[i].Rho +LiuParticles[j].Eta*tzz(j)/LiuParticles[j].Rho/LiuParticles[j].Rho)*data.DWDX(3,k);
	  }
	}
	data.INDVXDt(d,i) = data.INDVXDt(d,i) + LiuParticles[j].Mass*h;
	data.INDVXDt(d,j) = data.INDVXDt(d,j) - LiuParticles[i].Mass*h;
      }
      data.INDUDt(i) = data.INDUDt(i) + LiuParticles[j].Mass*he;
      data.INDUDt(j) = data.INDUDt(j) + LiuParticles[i].Mass*he;
    }
  }
 
//   Change of specific internal energy de/dt = T ds/dt - p/rho vc,c: 
 
  for (i=1;i<=NTotal+NVirt;i++) {
    data.INDUDt(i) = data.TDSDt(i) + 0.5e0*data.INDUDt(i);
  }
}

//-------------------------------------------------------------------------- 
//  Subroutine to calculate the external forces, e.g. gravitational forces.       
//  The forces from the interactions with boundary virtual particles  
//  are also calculated here as external forces. 
 
//  here as the external force.  
//  NTotal  : Number of particles                                 [in] 
//  data.Mass    : Particle masses                                     [in] 
//  data.X       : Coordinates of all particles                        [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  data.IType   : type of particles                                   [in] 
//  data.Hsml   : Smoothing Length                                     [in] 
//  data.EXDVXDt   : Acceleration with respect to x, y and z            [out]  

void MkLiuSPH::Ext_Force(SPH_Data &data) 
{ 

  int i, j, k, d ;
  double rr, f, rr0, dd, p1, p2;
  MkDouble dx(Dim+1);

   i= j= k= d =0;
   rr= f= rr0= dd= p1= p2=0;
            
  for (i = 1;i<= NTotal+NVirt;i++) { 
    for (d = 1;d<= Dim;d++){ 
      data.EXDVXDt(d, i) = 0.; 
    }
  }
         
//  Consider self-gravity or not ? 
 
  if (self_gravity) { 
    for (i = 1; i<=NTotal+NVirt;i++) { 
      data.EXDVXDt(Dim, i) = -9.8; 
    }
  }  
 
//  Boundary particle force and penalty anti-penetration force.  
  rr0 = 1.25e-5; 
  dd = 1.e-2; 
  p1 = 12; 
  p2 = 4; 
       
  for ( k=1;k<=NIac;k++){ 
    i = Pair_I(k); 
    j = Pair_J(k);   
    if(LiuParticles[i].ParticleType>0&&LiuParticles[j].ParticleType<0) {   
      rr = 0.;       
      /*
      for (d=1;d<=Dim;d++){ 
	dx(d) =  data.X(d,i) -  data.X(d,j); 
	rr = rr + dx(d)*dx(d); 
      }
      */
      dx(1) =  LiuParticles[i].X - LiuParticles[j].X; 
      rr = rr + dx(1)*dx(1); 
      dx(2) =  LiuParticles[i].Y - LiuParticles[j].Y; 
      rr = rr + dx(2)*dx(2); 
      dx(3) =  LiuParticles[i].Z - LiuParticles[j].Z; 
      rr = rr + dx(1)*dx(1); 

      rr = sqrt(rr); 
      if(rr<rr0) { 
	f = (pow(rr0/rr,p1)-pow(rr0/rr,p2))/rr*rr ;
	for (d = 1; d<=Dim;d++) { 
	  data.EXDVXDt(d, i) = data.EXDVXDt(d, i) + dd*dx(d)*f ;
	}
      } 
    }         
  }   
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the average velocity to correct velocity 
//  for preventing.penetration (monaghan, 1992) 
 
//  NTotal : Number of particles                                  [in] 
//  data.Mass   : Particle masses                                      [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  data.W      : Kernel for all interaction pairs                     [in] 
//  data.VX     : Velocity of each particle                            [in] 
//  data.Rho    : Density of each particle                             [in] 
//  data.AveVel     : Average velocityof each particle                    [out] 

void MkLiuSPH::Av_Vel(SPH_Data &data)
{ 
  int i,j,k,d;        
  double  vcc, epsilon; 
  MkDouble dvx(Dim+1);

   i=j=k=d=0;        
   vcc= epsilon=0; 
       
//  epsilon --- a small constants chosen by experience, may lead to instability. 
//  for example, for the 1 dimensional shock tube problem, the E = 0.3 
 
  epsilon = 0.3;
       
  for (i = 1;i<= NTotal;i++ ) { 
    for (d = 1;d<= Dim;d++) { 
      data.AveVel(d,i) = 0.; 
    }
  }
      
  for (k=1;k<=NIac;k++) {        
    i = Pair_I(k); 
    j = Pair_J(k);
    for (d=1;d<=Dim;d++) { 
      dvx(d) = data.VX(d,i) - data.VX(d,j);
      data.AveVel(d, i) = data.AveVel(d,i) - 2*LiuParticles[j].Mass*dvx(d)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*data.W(k);
      data.AveVel(d, j) = data.AveVel(d,j) + 2*LiuParticles[i].Mass*dvx(d)/(LiuParticles[i].Rho+LiuParticles[j].Rho)*data.W(k);                       
    }
  }         
  for (i = 1;i<= NTotal;i++) {
    for (d = 1;d<= Dim;d++){ 
      data.AveVel(d,i) = epsilon * data.AveVel(d,i); 
    }
  } 
}

//----------------------------------------------------------------------       
//  Subroutine to established a pair linked list by sorting grid cell. 
//  It is suitable for a homogeneous particle distribution with the  
//  same smoothing length in an instant. A fixed number of particles 
//  lie in each cell.  
 
//    NTotal   : Number of particles                                [in] 
//    data.Hsml     : Smoothing Length                                   [in] 
//    grid     : array of grid cells                               [out] 
//    ngridx   : Number of sorting grid cells in x, y, z-direction [out] 
//    ghsmlx   : Smoothing length measured in cells of the grid    [out] 
//    data.MaxGridX : Maximum x-, y- and z-coordinate of grid range     [out] 
//    grid.MinGridX : Minimum x-, y- and z-coordinate of grid range     [out] 
//    dgeomx   : x-, y- and z-expansion of grid range              [out] 
 
//    Parameter used for sorting grid cells in the link list algorithm 
//    grid.MaxNGridX  : Maximum number of sorting grid cells in x-direction 
//    grid.MaxNGridY  : Maximum number of sorting grid cells in y-direction 
//    grid.MaxNGridZ  : Maximum number of sorting grid cells in z-direction 
//    Determining maximum number of sorting grid cells: 
//    (For an homogeneous particle distribution:) 
//    1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1 
//    2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1 
//    3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3) 

void MkLiuSPH::Init_Grid(SPH_Data &data,SPH_Grid &grid) 
{
  int i, j, k, d;
  MkDouble maxng(Dim+1), ngrid(3+1);
  double nppg=0; 

   i= j= k= d=0;
 
//    Averaged number of particles per grid cell 
 
  nppg = 3.e0;
 
//    Initialize parameters: Maximum number of grid cells 
 
  maxng(1) = grid.MaxNGridX; 
  if (Dim>=2) {
    maxng(2) = grid.MaxNGridY;
    if (Dim==3) { 
      maxng(3) = grid.MaxNGridZ; 
    }
  }
       
  for (d=1;d<=3;d++) {
    ngrid(d) = 1;
  }
       
//    Range of sorting grid 
 
  grid.MaxGridX(1) = x_maxgeom; 
  grid.MinGridX(1) = x_mingeom; 
  if (Dim>=2) { 
    grid.MaxGridX(2) = y_maxgeom; 
    grid.MinGridX(2) = y_mingeom; 
    if (Dim==3) { 
      grid.MaxGridX(3) = z_maxgeom;
      grid.MinGridX(3) = z_mingeom;
    }
  }
 
  for (d=1;d<=Dim;d++) {
    grid.DGeomX(d) = grid.MaxGridX(d) - grid.MinGridX(d);
  }
 
//    Number of grid cells in x-, y- and z-direction: 
 
  if (Dim==1) {
    grid.NGridX(1) = min(int(NTotal/nppg) + 1,maxng(1));
  }
  else if (Dim==2) {
    grid.NGridX(1) = min(int(sqrt(NTotal*grid.DGeomX(1)/(grid.DGeomX(2)*nppg))) + 1,maxng(1)); 
    grid.NGridX(2) = min(int(grid.NGridX(1)*grid.DGeomX(2)/grid.DGeomX(1)) + 1,maxng(2));
  } 
  else if (Dim==3) {
    grid.NGridX(1) = min(int(pow(NTotal*grid.DGeomX(1)*grid.DGeomX(1)/(grid.DGeomX(2)*grid.DGeomX(3)*nppg),(1.e0/3.e0))) + 1,maxng(1)); 
    grid.NGridX(2) = min(int(grid.NGridX(1)*grid.DGeomX(2)/grid.DGeomX(1)) + 1,maxng(2)); 
    grid.NGridX(3) = min(int(grid.NGridX(1)*grid.DGeomX(3)/grid.DGeomX(1)) + 1,maxng(3)); 
  }
 
//    Smoothing Length measured in grid cells: 
 
  for (d=1;d<=Dim;d++) {
    grid.GHsmlX(d) = int(float(grid.NGridX(d))*data.Hsml(1)/grid.DGeomX(d)) + 1 ;
  }
 
  for (d=1;d<=Dim;d++) {
    ngrid(d) = grid.NGridX(d);
  }
 
//    Initialize grid 
 
  for (i=1;i<=ngrid(1);i++) {
    for (j=1;j<=ngrid(2);j++) { 
      for (k=1;k<=ngrid(3);k++) {
	grid.Grid(i,j,k) = 0; 
      }
    }
  }
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the coordinates (xgcell) of the cell of  
//  the sorting  grid, in which the particle with coordinates (x) lies. 
 
//    x        : Coordinates of particle                            [in]     
//    ngridx   : Number of sorting grid cells in x, y, z-direction  [in] 
//    maxgridx : Maximum x-, y- and z-coordinate of grid range      [in] 
//    mingridx : Minimum x-, y- and z-coordinate of grid range      [in] 
//    dgeomx   : x-, y- and z-expansion of grid range               [in] 
//    xgcell   : x-, y- and z-coordinte of sorting grid cell       [out] 
 
void MkLiuSPH::Grid_Geom(int i,MkDouble x,MkInt ngridx,MkDouble maxgridx,MkDouble mingridx,MkDouble dgeomx,MkInt &xgcell,SPH_Data &data) 
{
  int d; 
 
  for (d=1;d<=3;d++) { 
    xgcell(d) = 1; 
  }
 
  for (d=1;d<=Dim;d++) { 
    if ((x(d)>maxgridx(d))||(x(d)<mingridx(d))) { 
      printf(" >>> ERROR << : Particle out of range\n"); 
      printf("    Particle position: x(%d,%d) =%f \n",i,d,x(d)) ;
      printf("    Range: [xmin,xmax](%f) = [%f,%f]\n",d,mingridx(d),maxgridx(d));
      exit(-1);//stop
    } 
    else {  
      xgcell(d) = int(float(ngridx(d))/dgeomx(d)* (x(d)-mingridx(d)) + 1.e0) ;
    }
  }
}

//----------------------------------------------------------------------- 
//  Subroutine to evolve smoothing length 
 
//  Dt     : time step                                            [in] 
//  NTotal : Number of particles                                  [in] 
//  data.Mass   : Particle masses                                      [in] 
//  data.VX     : Velocities of all particles                          [in] 
//  data.Rho    : Density                                              [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  data.DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//  data.Hsml   : Smoothing Length                                 [in/out] 

void MkLiuSPH::H_Upgrade(SPH_Data &data) 
{    
  int i,j,k,d; 
  double fac, hvcc;
  MkDouble dvx(Dim+1), vcc(data.MaxN+1), dhsml(data.MaxN+1);

   i=j=k=d=0; 
   fac= hvcc=0;
 
  if (sle==0 ) {      
 //---  Keep smoothing length unchanged.  
     return; 
  }      
  else if (sle==2) {
       
//---  dh/dt = (-1/Dim)*(h/rho)*(drho/dt). 
 
    for (i=1;i<=NTotal ;i++) {
      vcc(i) = 0.e0 ;
    }
       
    for (k=1;k<=NIac;k++) {
      i = Pair_I(k); 
      j = Pair_J(k) ;
      for (d=1;d<=Dim;d++) {
	dvx(d) = data.VX(d,j) - data.VX(d,i) ;
      }
      hvcc = dvx(1)*data.DWDX(1,k);
      for (d=2;d<=Dim;d++) {
	hvcc = hvcc + dvx(d)*data.DWDX(d,k);
      }
      vcc(i) = vcc(i) + LiuParticles[j].Mass*hvcc/LiuParticles[j].Rho; 
      vcc(j) = vcc(j) + LiuParticles[i].Mass*hvcc/LiuParticles[i].Rho;
    }
         
    for (i = 1;i<= NTotal;i++) {
      dhsml(i) = (LiuParticles[i].SmoothLen/Dim)*vcc(i);
      LiuParticles[i].SmoothLen = LiuParticles[i].SmoothLen + Dt*dhsml(i);
      if (LiuParticles[i].SmoothLen<=0) LiuParticles[i].SmoothLen = LiuParticles[i].SmoothLen - Dt*dhsml(i);  
    }
  }
     
  else if(sle==1) {
    fac = 2.0;
    for (i = 1;i<= NTotal;i++) {
      LiuParticles[i].SmoothLen = fac * (LiuParticles[i].Mass/pow(LiuParticles[i].Rho,(1./Dim)));
    }
        
  }
}

//---------------------------------------------------------------------- 
//   Subroutine for loading or generating initial particle information 
 
//   data.X-- coordinates of particles                                 [out] 
//   data.VX-- velocities of particles                                 [out] 
//   data.Mass-- mass of particles                                     [out] 
//   data.Rho-- dnesities of particles                                 [out] 
//   p-- pressure  of particles                                   [out] 
//   data.U-- internal energy of particles                             [out] 
//   data.IType-- types of particles                                   [out] 
//   data.Hsml-- smoothing lengths of particles                        [out] 
//   ntotal-- total particle number                               [out] 
 
void MkLiuSPH::Input(SPH_Data &data) 
{
  int i, d, im, m, n, mp, np,;
  FILE *fp1,*fp2, *fp3;
 
//   load initial particle information from external disk file 
   int yesorno=1;      

   //  while(yesorno) {
    printf("  ***************************************************\n"); 
    printf("          Please input the maximal time steps \n");
    printf("  ***************************************************\n");
    scanf("%d", &MaxTimeStep);
    MkDebug("    maxtimestep = %d \n",MaxTimeStep);
    time_integration(data);
    output(data);
    printf("  ***************************************************\n"); 
    printf("  Are you going to run more time steps ? (0=No, 1=yes)\n"); 
    printf("  ***************************************************\n");
    scanf ("%d ", &yesorno);
    //  }

  MkDebug("input():");
  if(config_input){
    LoadFromFile();
  }      
  else  {
           
    //    MkDebug("not config \n");                        
    fp1 = fopen("ini_xv.dat","w") ;
    fp2 = fopen("ini_state.dat","w") ;
    fp3 = fopen("ini_other.dat","w")  ;

    if (shocktube) Shock_Tube();
    if (shearcavity) {
      m = 41;
      n = 41;
      mp = m-1;
      np = n-1;
      NTotal = mp * np;
      NVirt = 2*mp + 2*np-2;
      LiuParticles.Initialize(NTotal+NVirt);

      Shear_Cavity();
      Virt_Part();
    }

    MkDebug("shear cavity %d\n",shearcavity);
    for (i = 1;i<= NTotal ;i++) {
      fprintf(fp1,"%6d ", i);
      fprintf(fp1,"%14.8f %14.8f %14.8f ", LiuParticles[i].X, LiuParticles[i].Y, LiuParticles[i].Z);
      for (d=1;d<=Dim;d++) fprintf(fp1,"%14.8f ", data.VX(d, i));
      fprintf(fp1,"\n");
      fprintf(fp2,"%6d %14.8f %14.8f %14.8f %14.8f \n", i, LiuParticles[i].Mass, LiuParticles[i].Rho, LiuParticles[i].Press, LiuParticles[i].Energy);
      fprintf(fp3,"%6d %3d %14.8f \n", i, LiuParticles[i].ParticleType, LiuParticles[i].SmoothLen);
    }
    printf("  **************************************************\n");
    printf("      Initial particle configuration generated   \n");       
    printf("      Total number of particles %d\n ", NTotal);	 
    printf("  **************************************************\n");
  }
 
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  //  MkDebug("~input():\n");
}        

//----------------------------------------------------------------------            
//   Subroutine for saving particle information to external disk file 
 
//   data.X-- coordinates of particles                                  [in] 
//   data.VX-- velocities of particles                                  [in] 
//   data.Mass-- mass of particles                                      [in] 
//   data.Rho-- dnesities of particles                                  [in] 
//   data.P-- pressure  of particles                                    [in] 
//   data.U-- internal energy of particles                              [in] 
//   data.C-- sound velocity of particles                               [in] 
//   data.IType-- types of particles                                    [in] 
//   data.Hsml-- smoothing lengths of particles                         [in] 
//   NTotal-- total particle number                                [in] 
 
void MkLiuSPH::Output()  
{
  SaveToFile();
}


//----------------------------------------------------------------------      
//   This subroutine is used to generate initial data for the  
//   1 d noh shock tube problem 
//   Data.X-- coordinates of particles                                 [out] 
//   data.VX-- velocities of particles                                 [out] 
//   data.Mass-- mass of particles                                     [out] 
//   data.Rho-- dnesities of particles                                 [out] 
//   data.P-- pressure  of particles                                   [out] 
//   data.U-- internal energy of particles                             [out] 
//   data.IType-- types of particles                                   [out] 
//        =1   ideal gas 
//   data.Hsml-- smoothing lengths of particles                        [out] 
//   ntotal-- total particle number                               [out] 
        
void MkLiuSPH::Shock_Tube(SPH_Data &data) 
{
  int i, d; 
  double space_x;
 
  NTotal=400; 
  space_x=0.6/80.;       

  //  MkDebug("  shock_tube() : \n");

  for (i=1;i<=NTotal;i++) {
    LiuParticles[i].Mass=0.75/400.; 
    LiuParticles[i].SmoothLen=0.015; 
    LiuParticles[i].ParticleType=1;
    /*
    for (d = 1;d<= Dim;d++) {
      data.X(d,i) = 0.;
      data.VX(d,i) = 0. ;
    }
    */
    LiuParticles[i].X = LiuParticles[i].Y = LiuParticles[i].Z = 0;
    LiuParticles[i].XVel[0] = LiuParticles[i].YVel[0] = LiuParticles[i].ZVel[0] = 0;
  }

  for (i=1;i<=320;i++) {
    LiuParticles[i].X=-0.6+space_x/4.*(i-1); 
  }
       
  for (i=320+1;i<=NTotal;i++) {
    LiuParticles[i].X=0.+space_x*(i-320);
  }         

  for (i=1;i<=NTotal;i++) {
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
//   data.X-- coordinates of particles                                 [out] 
//   data.VX-- velocities of particles                                 [out] 
//   data.Mass-- mass of particles                                     [out] 
//   data.Rho-- dnesities of particles                                 [out] 
//   data.P-- pressure  of particles                                   [out] 
//   data.U-- internal energy of particles                             [out] 
//   data.IType-- types of particles                                   [out] 
//        =2   water 
//   h-- smoothing lengths of particles                           [out] 
//   NTotal-- total particle number                               [out] 
void MkLiuSPH::Shear_Cavity(SPH_Data &data) 
{
  int i, j, d, m, n, mp, np, k; 
  double xl, yl, dx, dy; 
 
//   Giving data.Mass and smoothing length as well as other data. 
 
//  MkDebug("shear_cavity()\n");
  m = 41; 
  n = 41; 
  mp = m-1; 
  np = n-1; 
  NTotal = mp * np; 
  xl = 1.e-3; 
  yl = 1.e-3; 
  dx = xl/mp; 
  dy = yl/np; 
 
  for (i = 1;i<=mp;i++) {
    for (j = 1;j<= np;j++) {
      k = j + (i-1)*np; 
      //data.X(1, k) = (i-1)*dx + dx/2.; 
      //data.X(2, k) = (j-1)*dy + dy/2. ;
      LiuParticles[k].X = (i-1)*dx + dx/2.0;
      LiuParticles[k].Y = (j-1)*dy + dy/2.0;
    }
  }
 
  for (i = 1;i<= mp*np ;i++) {
    data.VX(1, i) = 0.; 
    data.VX(2, i) = 0.;       
    data.Rho (i) = 1000.;    
    LiuParticles[i].Mass = dx*dy*LiuParticles[i].Rho;   
    LiuParticles[i].Press= 0.;    
    LiuParticles[i].Energy=357.1; 
    LiuParticles[i].ParticleType = 2; 
    LiuParticles[i].SmoothLen = dx; 
  }
  //  MkDebug("~shear_cavity()\n");
}

//---------------------------------------------------------------------- 
// Subroutine to calculate the smoothing funciton for each particle and 
// the interaction parameters used by the SPH algorithm. Interaction  
// pairs are determined by using a sorting grid linked list   
 
//   CurrentTimeStep : Current time step                                 [in] 
//   NTotal    : Number of particles                               [in] 
//   data.Hsml      : Smoothing Length, same for all particles          [in] 
//   data.X         : Coordinates of all particles                      [in] 
//   NIac      : Number of interaction pairs                      [out] 
//   Pair_I    : List of first partner of interaction pair        [out] 
//   Pair_J    : List of second partner of interaction pair       [out] 
//   data.W         : Kernel for all interaction pairs                 [out] 
//   data.DWDX      : Derivative of kernel with respect to x, y and z  [out] 
//   data.CountIac  : Number of neighboring particles                  [out] 
//   Parameter used for sorting grid cells in the link list algorithm 
//   grid.MaxNGridX  : Maximum number of sorting grid cells in x-direction 
//   grid.MaxNGridY  : Maximum number of sorting grid cells in y-direction 
//   grid.MaxNGridZ  : Maximum number of sorting grid cells in z-direction 
//   Determining maximum number of sorting grid cells: 
//   (For an homogeneous particle distribution:) 
//   1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1 
//   2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1 
//   3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3) 
void MkLiuSPH::Link_List(SPH_Data &data) 
{
  int i, j, d, scale_k, sumiac, maxiac, noiac, miniac, maxp,minp,xcell,ycell,zcell;
  double hsml2,dr,r;

  MkInt gcell(3+1), celldata(data.MaxN+1),minxcell(3+1),maxxcell(3+1), dnxgcell(Dim+1),dpxgcell(Dim+1);
  MkDouble dx(Dim+1),tdwdx(Dim+1), dgeomx(Dim+1);
  SPH_Grid grid; 

  int maxngx = 100, maxngy = 100, maxngz = 1;

   i= j= d= scale_k= sumiac= maxiac= noiac= miniac= maxp=minp=xcell=ycell=zcell=0;
   hsml2=dr=r=0;

  grid.Initialize(Dim,maxngx,maxngy,maxngz);

  if (skf==1) scale_k = 2;
  else if (skf==2) scale_k = 3; 
  else if (skf==3) scale_k = 3;
      
  for (i=1;i<=NTotal+NVirt;i++) {
    data.CountIac(i) = 0;
  }
 
//   Initialize grid:   
 
  init_grid(data,grid); 
       
//   Position particles on grid and create linked list: 
       
  for (i=1;i<=NTotal+NVirt;i++) {
    grid_geom(i,LiuParticles[i].X,grid.NGridX,grid.MaxGridX,grid.MinGridX,grid.DGeomX,gcell,data);
    for (d=1;d<=Dim;d++) {
      data.XGCell(d,i) = gcell(d);
    }
    celldata(i) = grid.Grid(gcell(1),gcell(2),gcell(3));
    grid.Grid(gcell(1),gcell(2),gcell(3)) = i; 
  }
 
//   Determine interaction parameters: 
 
  NIac = 0;
  for (i=1;i<=NTotal+NVirt-1;i++) {
 //   Determine range of grid to go through: 
    for (d=1;d<=3;d++) {
      minxcell(d) = 1;
      maxxcell(d) = 1;
    }
    for (d=1;d<=Dim;d++) {
      dnxgcell(d) = data.XGCell(d,i) - grid.GHsmlX(d);
      dpxgcell(d) = data.XGCell(d,i) + grid.GHsmlX(d);
      minxcell(d) = max(dnxgcell(d),1); 
      maxxcell(d) = min(dpxgcell(d),grid.NGridX(d));
    }
 //   Search grid: 
    for (zcell=minxcell(3);zcell<=maxxcell(3);zcell++) {
      for (ycell=minxcell(2);ycell<=maxxcell(2);ycell++) {
	for (xcell=minxcell(1);xcell<=maxxcell(1);xcell++) {
	  j = grid.Grid(xcell,ycell,zcell) ;
	  while (j>i) {
	    /*
	    for (d=1;d<=Dim;d++) {
	      dx(d) = data.X(d,i) - data.X(d,j);
	      dr    = dr + dx(d)*dx(d);
	    }
	    */

	    dx(1) = LiuParticles[i].X - LiuParticles[j].X;
	    dr    = dx(1)*dx(1); 
	    dx(2) = LiuParticles[i].Y - LiuParticles[j].Y;
	    dr    = dr+ dx(2)*dx(2); 
	    dx(3) = LiuParticles[i].Z - LiuParticles[j].Z;
	    dr    = dr+dx(3)*dx(3); 

	    if (sqrt(dr)<scale_k*data.Hsml(1)) {
	      if (NIac<max_interaction) {
 
//   Neighboring pair list, and totalinteraction number and 
//   the interaction number for each particle  
 
		NIac = NIac + 1; 
		Pair_I(NIac) = i; 
		Pair_J(NIac) = j; 
		r = sqrt(dr); 
		data.CountIac(i) = data.CountIac(i) + 1; 
		data.CountIac(j) = data.CountIac(j) + 1; 
                            
//--- Kernel and derivations of kernel 
 
//		kernel(r,dx,data.Hsml(1),data.W(NIac),tdwdx,data); 
		data.W(NIac) = LiuKernel.W(r);
		for (d = 1;d<= Dim;d++) {
		  data.DWDX(d,NIac)= LiuKernel.dWdX(r,dx,d-1);
		}
	      }      
	      else {
		printf (" >>> Error << : too many interactions");
		exit(-1);//stop 
	      }
	    }
	    j = celldata(j);
	  }
	} 
      }
    }
  }
  
 
//   Statistics for the interaction 
 
  sumiac = 0;
  maxiac = 0; 
  miniac = 1000; 
  noiac  = 0; 
  for (i=1;i<=NTotal+NVirt;i++) {
    sumiac = sumiac + data.CountIac(i);
    if (data.CountIac(i)>maxiac) { 
      maxiac = data.CountIac(i);
      maxp = i;
    }
    if (data.CountIac(i)<miniac) {
      miniac = data.CountIac(i); 
      minp = i;
    }
    if (data.CountIac(i)==0) noiac  = noiac + 1; 
  }
  
  if ((CurrentTimeStep%print_step)==0) {
    if (int_stat) {
      printf(" >> Statistics: interactions per particle:\n");
      printf("**** Particle:%f, maximal interactions:%f\n",maxp, maxiac); 
      printf("**** Particle:%f, minimal interactions:%f\n",minp, miniac); 
      printf("**** Average :%f\n",float(sumiac)/float(NTotal+NVirt)); 
      printf("**** Total pairs : %f",NIac); 
      printf("**** Particles with no interactions: %f",noiac); 
    }
  }
}


 
 

//---------------------------------------------------------------------- 
//  Subroutine to calculate the density with SPH summation algorithm. 
//  See Equ.(4.35) 
 
//    NTotal : Number of particles                                  [in] 
//    data.Hsml   : Smoothing Length                                     [in] 
//    data.Mass   : Particle masses                                      [in] 
//    NIac   : Number of interaction pairs                          [in] 
//    Pair_I : List of first partner of interaction pair            [in] 
//    Pair_J : List of second partner of interaction pair           [in] 
//    data.W      : Kernel for all interaction pairs                     [in] 
//    data.IType   : type of particles                                   [in] 
//    data.X       : Coordinates of all particles                        [in] 
//    data.Rho    : Density                                             [out] 

void MkLiuSPH::Sum_Density(SPH_Data &data) 
{
  int i, j, k, d;
  double selfdens, r;
  MkDouble hv(Dim+1), wi(data.MaxN+1);

   i= j= k= d=0;
   selfdens= r=0;
 
   //  MkDebug("\n      sum_density()\n");
//    wi(data.MaxN)---integration of the kernel itself 
//   MkDebug("      SD1 ");         
  for (d=1;d<=Dim;d++) { 
    hv(d) = 0.e0; 
  }
  //  MkDebug("SD2 ");          
//     Self density of each particle: Wii (Kernel for distance 0) 
//     and take contribution of particle itself: 
 
  r=0.;
       
//     Firstly calculate the integration of the kernel over the space 
 
  for (i=1;i<=NTotal+NVirt ;i++) {
    //    kernel(r,hv,LiuParticles[i].SmoothLen,selfdens,hv,data) ;
    selfdens = LiuKernel.W(r);
    wi(i)=selfdens*LiuParticles[i].Mass/LiuParticles[i].Rho; 
  }
  //  MkDebug("SD3 ");           
  for (k=1;k<=NIac;k++) { 
    i = Pair_I(k); 
    j = Pair_J(k); 
    wi(i) = wi(i) + LiuParticles[j].Mass/LiuParticles[j].Rho*data.W(k); 
    wi(j) = wi(j) + LiuParticles[i].Mass/LiuParticles[i].Rho*data.W(k); 
  }
 
//  Secondly calculate the rho integration over the space 
//  MkDebug("SD4 ");           
  for (i=1;i<=NTotal+NVirt;i++) { 
    //    kernel(r,hv,LiuParticles[i].SmoothLen,selfdens,hv,data); 
    selfdens = LiuKernel.W(r);
    LiuParticles[i].Rho = selfdens*LiuParticles[i].Mass; 
  }
  //  MkDebug("SD5 ");           
//  Calculate SPH sum for rho: 
  for (k=1;k<=NIac;k++) { 
    i = Pair_I(k); 
    j = Pair_J(k); 
    LiuParticles[i].Rho = LiuParticles[i].Rho + LiuParticles[j].Mass*data.W(k); 
    LiuParticles[j].Rho = LiuParticles[j].Rho + LiuParticles[i].Mass*data.W(k); 
  }
  //  MkDebug("SD6 ");           

//  Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w) 
      
  if (nor_density) {  

    for (i=1;i<= NTotal+NVirt;i++) { 
      LiuParticles[i].Rho=LiuParticles[i].Rho/wi(i);
    }
  }  
  //  MkDebug("\n      ~sum_density()\n");       
}
       
//---------------------------------------------------------------------- 
//  Subroutine to calculate the density with SPH continuity approach. 
//  See Equ.(4.34) 
 
//  NTotal : Number of particles                                  [in] 
//  Mass   : Particle masses                                      [in] 
//  NIac   : Number of interaction pairs                          [in] 
//  Pair_I : List of first partner of interaction pair            [in] 
//  Pair_J : List of second partner of interaction pair           [in] 
//  data.DWDX   : derivation of Kernel for all interaction pairs       [in] 
//  data.VX     : Velocities of all particles                          [in] 
//  data.IType   : type of particles                                   [in] 
//  data.X      : Coordinates of all particles                         [in] 
//  Rho    : Density                                              [in] 
//  data.DRhoDt : Density change rate of each particle                [out]    

void MkLiuSPH::Con_Density(SPH_Data &data)
{ 
       
  int i,j,k,d;
  double vcc;
  MkDouble dvx(Dim+1);

   i=j=k=d=0;
   vcc=0;
       
  for (i = 1;i<= NTotal+NVirt;i++) { 
    data.DRhoDt(i) = 0.; 
  }
      
  for (k=1;k<=NIac;k++) {       
    i = Pair_I(k); 
    j = Pair_J(k); 
    for (d=1;d<=Dim;d++) { 
      dvx(d) = data.VX(d,i) - data.VX(d,j);
    }
    vcc = dvx(1)*data.DWDX(1,k);
    for (d=2;d<=Dim;d++) { 
      vcc = vcc + dvx(d)*data.DWDX(d,k); 
    }
    data.DRhoDt(i) = data.DRhoDt(i) + LiuParticles[j].Mass*vcc; 
    data.DRhoDt(j) = data.DRhoDt(j) + LiuParticles[i].Mass*vcc;
   
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
//    data.Hsml      : Smoothing Length                                  [in] 
//    data.X         : Coordinates of all particles                      [in] 
//    NIac      : Number of interaction pairs                      [out] 
//    Pair_I    : List of first partner of interaction pair        [out] 
//    Pair_J    : List of second partner of interaction pair       [out] 
//    data.W         : Kernel for all interaction pairs                 [out] 
//    data.DWDX      : Derivative of kernel with respect to x, y and z  [out] 
//    data.CountIac  : Number of neighboring particles                  [out] 

void MkLiuSPH::Direct_Find(SPH_Data &data) 
{ 
  int i, j, d,  sumiac, maxiac, miniac, noiac, maxp, minp, scale_k;  
  double  driac, r, mhsml;      
  MkDouble dxiac(Dim+1), tdwdx(Dim+1);

   i= j= d=  sumiac= maxiac= miniac= noiac= maxp= minp= scale_k=0;  
   driac= r= mhsml=0;      

//    Smoothing kernel function  
//    skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//        = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//        = 3, Quintic kernel (Morris 1997) 

MkDebug("    direct_find()\n");
MkDebug("    D1 ");
  if (skf==1) {  
    scale_k = 2;
  }  
  else if (skf==2){
    scale_k = 3 ; 
  }
  else if (skf==3) {
    scale_k = 3 ; 
  }
      
  for (i=1;i<=NTotal+NVirt;i++){
    data.CountIac(i) = 0; 
  }
  MkDebug("D2 ");
       
  NIac = 0;
 
  for (i=1;i<=NTotal+NVirt-1;i++) {      
    for (j = i+1;j<= NTotal+NVirt;j++){ 
      /*
      for (d=1;d<=Dim;d++) { 
	dxiac(d) = data.X(d,i) - data.X(d,j); 
	driac    = driac + dxiac(d)*dxiac(d); 
      }
      */

      dxiac(1) = LiuParticles[i].X - LiuParticles[j].X; 
      driac    = dxiac(1)*dxiac(1); 
      dxiac(2) = LiuParticles[i].Y - LiuParticles[j].Y; 
      driac    = driac + dxiac(2)*dxiac(2); 
      dxiac(3) = LiuParticles[i].Z - LiuParticles[j].Z; 
      driac    = driac + dxiac(3)*dxiac(3); 

      MkDebug("-1 ");


      MkDebug("-2 ");

      mhsml = (LiuParticles[i].SmoothLen+LiuParticles[j].SmoothLen)/2.; 
      if (sqrt(driac)<scale_k*mhsml) { 
	if (NIac<max_interaction) {     

      MkDebug("-3 ");
 
//    Neighboring pair list, and totalinteraction number and 
//    the interaction number for each particle  

	  NIac = NIac + 1; 
	  Pair_I(NIac) = i; 
	  Pair_J(NIac) = j; 
	  r = sqrt(driac); 
	  data.CountIac(i) = data.CountIac(i) + 1; 
	  data.CountIac(j) = data.CountIac(j) + 1; 

	  MkDebug("-4 ");

//    Kernel and derivations of kernel 
 
//	  kernel(r,dxiac,mhsml,data.W(NIac),tdwdx,data) ;

	  MkDebug("Why I cannot see this. ");

	  MkDebug("-5 ");
	  data.W(NIac) = LiuKernel.W(r);
	  for (d=1;d<=Dim;d++) { 
	    data.DWDX(d,NIac) = LiuKernel.dWdX(r,dxiac,d-1);
	  }
      MkDebug("-6 ");
	}
	else {
	  printf(" >>> ERROR << : Too many interactions\n");  
	  exit(-1);//stop 
	}
      } 
    }
  }
  MkDebug("D3 "); 
//    Statistics for the interaction 
 
  sumiac = 0; 
  maxiac = 0; 
  miniac = 1000; 
  noiac  = 0; 
  for (i=1;i<=NTotal+NVirt;i++) { 
    sumiac = sumiac + data.CountIac(i); 
    if (data.CountIac(i)>maxiac){ 
      maxiac = data.CountIac(i); 
      maxp = i; 
    }
    if (data.CountIac(i)<miniac) {
      miniac = data.CountIac(i); 
      minp = i; 
    }
    if (data.CountIac(i)==0) noiac  = noiac + 1;
  }
  
  if ((CurrentTimeStep%print_step)==0) {
    if (int_stat) {
      printf(" >> Statistics: interactions per particle:\n"); 
      printf("**** Particle: %d  maximal interactions: %d\n",maxp,maxiac) ;
      printf("**** Particle: %d minimal interactions: %d\n",minp,miniac) ;
      printf("**** Average :%f\n",float(sumiac)/float(NTotal+NVirt)) ;
      printf("**** Total pairs : %d\n",NIac) ;
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
//   data.Hsml   : Smoothing Length                                 [in|out] 
//   data.Mass   : Particle masses                                  [in|out] 
//   data.X      : Coordinates of all particles                     [in|out] 
//   data.VX     : Velocities of all particles                      [in|out] 
//   data.Rho    : Density                                          [in|out] 
//   data.U      : internal energy                                  [in|out] 
//   data.IType   : type of particles                               [in|out] 
 
void MkLiuSPH::Virt_Part()  
{
  int i, j, d, im, mp; 
  double xl, dx, v_inf; 
  FILE *fp1,*fp2, *fp3;

   i= j= d= im= mp=0; 
   xl= dx= v_inf=0; 
 
  if (vp_input){
    LoadVPFromFile();
  }
  else  {
    nvirt = 0; 
    mp = 40;
    xl = 1.0e-3; 
    dx = xl / mp; 
    v_inf = 1.e-3; 
 
//   Monaghan type virtual particle on the Upper side 
 
    for (i = 1;i<= 2*mp+1;i++) {
      nvirt = nvirt + 1;
      data.X(1, NTotal + nvirt) = (i-1)*dx/2;
      data.X(2, NTotal + nvirt) = xl;   
      data.VX(1, NTotal + nvirt) = v_inf; 
      data.VX(2, NTotal + nvirt) = 0.; 
    }
 
//   Monaghan type virtual particle on the Lower side 
 
    for (i = 1;i<= 2*mp+1;i++) {
      nvirt = nvirt + 1; 
      data.X(1, NTotal + nvirt) = (i-1)*dx/2;  
      data.X(2, NTotal + nvirt) = 0.;   
      data.VX(1, NTotal + nvirt) = 0.; 
      data.VX(2, NTotal + nvirt) = 0.; 
    }
 
//   Monaghan type virtual particle on the Left side 
 
    for (i = 1;i<= 2*mp-1;i++) {
      nvirt = nvirt + 1;
      data.X(1, NTotal + nvirt) = 0.;
      data.X(2, NTotal + nvirt) = i*dx/2;
      data.VX(1, NTotal + nvirt) = 0.; 
      data.VX(2, NTotal + nvirt) = 0.;
    }
 
//   Monaghan type virtual particle on the Right side 
 
    for (i = 1;i<= 2*mp-1;i++) {
      nvirt = nvirt + 1;
      data.X(1, NTotal + nvirt) = xl;
      data.X(2, NTotal + nvirt) = i*dx/2;
      data.VX(1, NTotal + nvirt) = 0.; 
      data.VX(2, NTotal + nvirt) = 0.; 
    }
 
    for (i = 1;i<= nvirt;i++) {
      data.Rho (NTotal + i) = 1000.;
      data.Mass(NTotal + i) = data.Rho (NTotal + i) * dx * dx; 
      data.P(NTotal + i) = 0.; 
      data.U(NTotal + i) = 357.1; 
      data.ParticleType(NTotal + i) = -2; 
      data.Hsml(NTotal + i) = dx; 
    }
  } 
  if ((CurrentTimeStep%save_step)==0) {
    SaveVPToFile();
  }
 
  if ((CurrentTimeStep%print_step)==0) {
    if (int_stat) {
      printf(" >> Statistics: Virtual boundary particles:\n"); 
      printf("          Number of virtual particles:%d\n",NVirt);
    }
  }
}
//---------------------------------------------------------------------- 
// Subroutine to define the fluid particle viscosity 
  
//   NTotal  : Number of particles                                 [in] 
//   data.IType    : Type of particle                                   [in] 
//   data.X       : Coordinates of all particles                        [in] 
//   data.Rho     : Density                                             [in] 
//   data.Eta     : Dynamic viscosity                                  [out] 

void MkLiuSPH::Viscosity() 
{       
  int i; 
  for (i=1;i<=NTotal+NVirt;i++) {
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

void MkLiuSPH::P_ideal_gas(double rho, double u) 
{       
  float gamma=1.4;
           
//   For air (idea gas) 
//   See Equ.(3.82) 
 
  p = (gamma-1) * rho * u;      
  return p;
}

float MkLiuSPH::C_ideal_gas(double rho, double u) 
{       
  float gamma=1.4;
           
//   For air (idea gas) 
//   See Equ.(3.82) 
 
  c = sqrt((gamma-1) * u);  
  return c;
}

//---------------------------------------------------------------------- 
//   Artificial equation of state for the artificial compressibility  
 
//  rho    : Density                                              [in] 
//  u      : Internal energy                                      [in] 
//  p      : Pressure                                            [out] 
//  c      : sound velocity                                      [out] 
//  Equation of state for artificial compressibility    
   
    
void MkLiuSPH::P_art_water(double rho) 
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
  float c = 0.01; 
  p = c*c * rho;
  return p;
}

void MkLiuSPH::C_art_water(double rho) 
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
  float c = 0.01; 
  return c;
}


//---------------------------------------------------------------------- 
// Subroutine to calculate the smoothing kernel wij and its  
// derivatives dwdxij. 
//   if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//          = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//          = 3, Quintic kernel (Morris 1997) 
 
//   r    : Distance between particles i and j                     [in] 
//   dx   : x-, y- and z-distance between i and j                  [in]   
//   data.Hsml : Smoothing length                                       [in] 
//   data.W    : Kernel for all interaction pairs                      [out] 
//   data.DWDX : Derivative of kernel with respect to x, y and z       [out] 
 
//    kernel(r,dx   ,data.Hsml(1),data.W(NIac),tdwdx,data); LiuKernel.W(r); LiuKernel.dWdX(r,dx,0);
//    kernel(r,hv   ,LiuParticles[i].SmoothLen,selfdens         ,hv   ,data) ;
//    kernel(r,hv   ,LiuParticles[i].SmoothLen,selfdens         ,hv   ,data);
//    kernel(r,dxiac,mhsml       ,data.W(NIac),tdwdx,data) ; 

/*
void MkLiuSPH::Kernel(double r,MkDouble dx,double hsml,double &w,MkDouble &dwdx,SPH_Data &data)    
{       
  int i, j, d;
  double q, dw, factor;

   i= j= d=0;
   q= dw= factor=0;
 
  q = r/hsml;
  w = 0.e0;
  for (d=1;d<=Dim;d++) {
    dwdx(d) = 0.e0 ;
  }

  if (skf==1) {
  MkDebug("K1 ");
       
    if (Dim==1) factor = 1.e0/hsml;
    else if (Dim==2) factor = 15.e0/(7.e0*pi*hsml*hsml);
    else if (Dim==3) factor = 3.e0/(2.e0*pi*hsml*hsml*hsml);
    else  {
      printf(" >>> Error << : Wrong dimension: Dim =%f\n",Dim);
      exit(0);//stop 
    }
       
    if (q>=0&&q<=1.e0) {
      w = factor * (2./3. - q*q + q*q*q / 2.);
      for (d = 1;d<= Dim;d++) {
	dwdx(d) = factor * (-2.+3./2.*q)/hsml/hsml * dx(d);
      }
    }
    else if (q>1.e0&&q<=2) {
      w = factor * 1.e0/6.e0 * (2.-q)*(2.-q)*(2.-q);
      for (d = 1;d<= Dim;d++) {
	dwdx(d) =-factor * 1.e0/6.e0 * 3.*(2.-q)*(2.-q)/hsml * (dx(d)/r);
      }
    }
    else {
      w=0.;
      for (d= 1;d<= Dim;d++) {
	dwdx(d) = 0.;
      }
    }
  }
  else if (skf==2) {
  MkDebug("K2 ");
       
    factor = 1.e0 / (pow(hsml,Dim) * pow(pi,(Dim/2.)));
    if(q>=0&&q<=3) {
      w = factor * exp(-q*q);
      for (d = 1;d<= Dim;d++) {
	dwdx(d) = w * ( -2.* dx(d)/hsml/hsml);
      }
    }
    else {
      w = 0.;
      for (d = 1;d<= Dim;d++) {
	dwdx(d) = 0.;
      }	    
    }
  }
	 
  else if (skf==3) {
  MkDebug("K3 ");
    if (Dim==1) factor = 1.e0 / (120.e0*hsml);
    else if (Dim==2) factor = 7.e0 / (478.e0*pi*hsml*hsml);
    else if (Dim==3) factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml) ;
    else {
      printf(" >>> Error << : Wrong dimension: Dim =%f\n",Dim);
      exit(-1);//stop 
    }
             
    if(q>=0&&q<=1) {
      w = factor * ( pow(3-q,5) - 6*pow(2-q,5) + 15*pow(1-q,5) );
      for (d= 1;d<= Dim;d++) {
	dwdx(d) = factor * ( (-120 + 120*q - 50*q*q)/ hsml/ hsml * dx(d) ); 
      }
    }
    else if(q>1&&q<=2) {
      w = factor * ( pow(3-q,5) - 6*pow(2-q,5) );
      for (d= 1;d<= Dim;d++) {
	dwdx(d) = factor * (-5*pow((3-q),4) + 30*pow(2-q,4))/ hsml * (dx(d)/r);
      }
    }
    else if(q>2&&q<=3) {
      w = factor * pow(3-q,5);
      for (d= 1;d<= Dim;d++) {
	dwdx(d) = factor * (-5*pow(3-q,4)) / hsml * (dx(d)/r);
      }
    }
    else {
      w = 0.;
      for (d = 1;d<= Dim;d++) {
	dwdx(d) = 0.;
      }
    }
  }		 
  MkDebug("K4 it is really strange ");
  return;
}
*/
 
