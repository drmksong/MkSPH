#include "source.h"

SPH_Data::SPH_Data()
{
  Clear();
}

SPH_Data::~SPH_Data()
{
  Clear();
}

void SPH_Data::Clear()
{
  NTotal=0;      
  MaxN=0;        
  MaxInteraction = 0;
  MaxTimeStep=0; 
  NIac=0;        
  ITimeStep=0;   
  NVirt=0;
  Dim = 0;
  Dt=0;       
  CountIac.Clear(); 
  X.Clear();      
  VX.Clear();     
  AveVel.Clear(); 
  IType.Clear();     
  Mass.Clear();   
  Rho.Clear();    
  Eta.Clear();    
  P.Clear();      
  T.Clear();      
  U.Clear();      
  C.Clear();      
  Hsml.Clear();   
  XGCell.Clear();

  DEDt.Clear();   
  DRhoDt.Clear(); 
  TDSDt.Clear();  
  DUDt.Clear();   

  DVXDt.Clear();  

  Pair_I.Clear();    
  Pair_J.Clear();    
  W.Clear();      
  DWDX.Clear();   
}

void SPH_Data::Initialize(int maxn, int dim)
{
  MaxN = maxn;
  Dim = dim;
  MaxInteraction = 100*MaxN;
  CountIac.Initialize(MaxN+1); 
  X.Initialize(Dim+1, MaxN+1);      
  VX.Initialize(Dim+1, MaxN+1);     
  AveVel.Initialize(Dim+1, MaxN+1); 
  IType.Initialize(MaxN+1);     
  Mass.Initialize(MaxN+1);   
  Rho.Initialize(MaxN+1);    
  Eta.Initialize(MaxN+1);    
  P.Initialize(MaxN+1);      
  T.Initialize(MaxN+1);      
  U.Initialize(MaxN+1);      
  C.Initialize(MaxN+1);      
  Hsml.Initialize(MaxN+1);
  XGCell.Initialize(3+1,MaxN+1);

  DEDt.Initialize(MaxN+1);   
  DRhoDt.Initialize(MaxN+1); 
  TDSDt.Initialize(MaxN+1);  
  DUDt.Initialize(MaxN+1);   

  DVXDt.Initialize(Dim+1, MaxN+1);  
   
  W.Initialize(MaxInteraction+1);      
  Pair_I.Initialize(MaxInteraction+1);    
  Pair_J.Initialize(MaxInteraction+1);    
  DWDX.Initialize(Dim+1, MaxInteraction+1);   
}

void SPH_Data::Dump()
{

  MkDebug("NTotal= %d\n",NTotal);      
  MkDebug("MaxTimeStep=%d\n",MaxTimeStep); 
  MkDebug("NIac=%d\n", NIac);        
  MkDebug("ITimeStep=%d\n",ITimeStep);   
  MkDebug("NVirt=%d\n",NVirt);
  MkDebug("Dt=%f\n",Dt);       

  MkDebug("MaxN = %d\n",MaxN);
  MkDebug("Dim = %d\n",Dim);
  MkDebug("MaxInteraction = %d\n",MaxInteraction);
  //  MkDebug("CountIac = %d\n",CountIac);
  MkDebug("X : xsz %d, ysz %d, zsz %d\n", X.getSzX(), X.getSzY(), X.getSzZ());
  MkDebug("X(1600) : xp %f, yp %f, zp %f\n", X(1,1600), X(2,1600), X(3,1600));
  MkDebug("VX : xsz %d, ysz %d, zsz %d\n", VX.getSzX(), VX.getSzY(), VX.getSzZ());
  MkDebug("VX(1600) : xv %f, yv %f, zv %f\n", VX(1,1600), VX(2,1600), VX(3,1600));
  MkDebug("AveVel : xsz %d, ysz %d, zsz %d\n", AveVel.getSzX(), AveVel.getSzY(), AveVel.getSzZ());
  MkDebug("AveVel(1600) : xav %f, yav %f, zav %f\n", AveVel(1,1600), AveVel(2,1600), AveVel(3,1600));
  /*
  MkDebug(VX);
  MkDebug(AveVel);
  MkDebug(IType);
  MkDebug(Mass);
  MkDebug(Rho);
  MkDebug(Eta);
  MkDebug(P);
  MkDebug(T);
  MkDebug(U);
  MkDebug(C);
  MkDebug(Hsml);
  MkDebug(XGCell);

  MkDebug(DEDt);
  MkDebug(DRhoDt);
  MkDebug(TDSDt);
  MkDebug(DUDt);

  MkDebug(DVXDt);
   
  MkDebug(,W);
  MkDebug(Pair_I);
  MkDebug(Pair_J);
  MkDebug(DWDX);
  */
}


SPH_Grid::SPH_Grid()
{
  Clear();     
}

SPH_Grid::~SPH_Grid()
{
  Clear();     
}

void SPH_Grid::Clear()
{
  Dim = 0;
  MaxNGridX=0;
  MaxNGridY=0;
  MaxNGridZ=0;
  NGridX.Clear();      
  MaxGridX.Clear();
  MinGridX.Clear(); 
  DGeomX.Clear();  
  GHsmlX.Clear();   
  Grid.Clear();     
}

void SPH_Grid::Initialize(int dim, int maxngx, int maxngy, int maxngz)
{
  Dim = dim;
  MaxNGridX=maxngx;
  MaxNGridY=maxngy;
  MaxNGridZ=maxngz;
  NGridX.Initialize(Dim+1);      
  MaxGridX.Initialize(Dim+1);
  MinGridX.Initialize(Dim+1);
  DGeomX.Initialize(Dim+1);
  GHsmlX.Initialize(Dim+1);
  Grid.Initialize(MaxNGridX+1, MaxNGridY+1, MaxNGridZ+1);
}

void SPH_Grid::Dump()
{

}


//---------------------------------------------------------------------- 
//  Subroutine to calculate the artificial heat(Fulk, 1994, p, a-17)  
//  See Equ.(4.74)   
 
//  data.NTotal : Number of particles                                  [in] 
//  data.Hsml   : Smoothing Length                                     [in] 
//  data.Mass   : Particle masses                                      [in] 
//  data.X      : Coordinates of all particles                         [in] 
//  data.VX     : Velocities of all particles                          [in] 
//  data.Rho    : Density                                              [in] 
//  data.U      : specific internal energy                             [in] 
//  data.C      : Sound veolcity                                       [in] 
//  data.NIac   : Number of interaction pairs                          [in] 
//  data.Pair_I : List of first partner of interaction pair            [in] 
//  data.Pair_J : List of second partner of interaction pair           [in] 
//  data.W      : Kernel for all interaction pairs                     [in] 
//  data.DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//  data.DEDt   : produced artificial heat, adding to energy Eq.      [out] 

void art_heat(SPH_Data &data)
{ 
  int i,j,k,d ;
  double dx, vr, rr, h, mc, mrho, mhsml, hvcc, mui, muj, muij, rdwdx, g1,g2;
  MkDouble vcc(data.NTotal+1),dvx(data.Dim+1);

  i=j=k=d =0;
  dx= vr= rr= h= mc= mrho= mhsml= hvcc= mui= muj= muij= rdwdx= g1=g2=0;
       
  //---  Parameter for the artificial heat conduction:
      
  g1=0.1; 
  g2=1.0;
  for (i=1;i<=data.NTotal+data.NVirt;i++) { 
    vcc(i) = 0.e0;
    data.DEDt(i) = 0.e0; 
  }
     
  for (k=1;k<=data.NIac;k++) { 
    i = data.Pair_I(k); 
    j = data.Pair_J(k); 
    for (d=1;d<=data.Dim;d++) { 
      dvx(d) = data.VX(d,j) - data.VX(d,i);  
    }        
    hvcc = dvx(1)*data.DWDX(1,k);
    for (d=2;d<=data.Dim;d++) { 
      hvcc = hvcc + dvx(d)*data.DWDX(d,k) ;
    }    
    vcc(i) = vcc(i) + data.Mass(j)*hvcc/data.Rho(j); 
    vcc(j) = vcc(j) + data.Mass(i)*hvcc/data.Rho(i); 
  }
    
  for (k=1;k<=data.NIac;k++) { 
    i = data.Pair_I(k); 
    j = data.Pair_J(k); 
    mhsml= (data.Hsml(i)+data.Hsml(j))/2.; 
    mrho = 0.5e0*(data.Rho(i) + data.Rho(j));
    rr = 0.e0; 
    rdwdx = 0.e0; 
    for (d=1;d<=data.Dim;d++) { 
      dx = data.X(d,i) -  data.X(d,j);
      rr = rr + dx*dx; 
      rdwdx  = rdwdx + dx*data.DWDX(d,k);
    }
    mui=g1*data.Hsml(i)*data.C(i) + g2*data.Hsml(i)*data.Hsml(i)*(fabs(vcc(i))-vcc(i));
    muj=g1*data.Hsml(j)*data.C(j) + g2*data.Hsml(j)*data.Hsml(j)*(fabs(vcc(j))-vcc(j));
    muij= 0.5*(mui+muj);
    h = muij/(mrho*(rr+0.01*mhsml*mhsml))*rdwdx ;
    data.DEDt(i) = data.DEDt(i) + data.Mass(j)*h*(data.U(i)-data.U(j)) ;
    data.DEDt(j) = data.DEDt(j) + data.Mass(i)*h*(data.U(j)-data.U(i)); 
  }
 
  for (i=1;i<=data.NTotal+data.NVirt;i++) { 
    data.DEDt(i) = 2.0e0*data.DEDt(i);           
  }
}

 
//---------------------------------------------------------------------- 
//  Subroutine to calculate the artificial viscosity (Monaghan, 1992)  
//  See Equ.(4.66) Equ.(4.62) 
 
//  data.NTotal : Number of particles (including virtual particles)    [in] 
//  data.Hsml   : Smoothing Length                                     [in] 
//  data.Mass   : Particle masses                                      [in] 
//  data.X      : Coordinates of all particles                         [in] 
//  data.VX     : Velocities of all particles                          [in] 
//  data.NIac   : Number of interaction pairs                          [in] 
//  data.Rho    : Density                                              [in] 
//  data.C      : Temperature                                          [in] 
//  data.Pair_I : List of first partner of interaction pair            [in] 
//  data.Pair_J : List of second partner of interaction pair           [in] 
//  data.W      : Kernel for all interaction pairs                     [in] 
//  data.DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//  data.DVXDt  : Acceleration with respect to x, y and z             [out]  
//  data.DEDt   : Change of specific internal energy                  [out] 

void art_visc(SPH_Data &data)
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
            
  for (i=1;i<=data.NTotal+data.NVirt;i++){ 
    for (d=1;d<=data.Dim;d++) { 
      data.DVXDt(d,i) = 0.e0;
    }
    data.DEDt(i) = 0.e0;
  }
      
//  Calculate SPH sum for artificial viscosity 
       
  for (k=1;k<=data.NIac;k++) { 
    i = data.Pair_I(k); 
    j = data.Pair_J(k); 
    mhsml= (data.Hsml(i)+data.Hsml(j))/2.; 
    vr = 0.e0; 
    rr = 0.e0; 
    for (d=1;d<=data.Dim;d++) { 
      dvx(d) = data.VX(d,i) - data.VX(d,j); 
      dx     =  data.X(d,i) -  data.X(d,j); 
      vr     = vr + dvx(d)*dx; 
      rr     = rr + dx*dx; 
    }
 
//  Artificial viscous force only if v_ij * r_ij  0 
 
    if (vr<0.e0) { 
 
//  Calculate muv_ij = data.Hsml v_ij * r_ij / ( r_ij^2 + data.Hsml^2 etq^2 ) 
             
      muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq);
           
//  Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij 
 
      mc   = 0.5e0*(data.C(i) + data.C(j)); 
      mrho = 0.5e0*(data.Rho(i) + data.Rho(j)); 
      piv  = (beta*muv - alpha*mc)*muv/mrho;
 
//  Calculate SPH sum for artificial viscous force 
 
      for (d=1;d<=data.Dim;d++) { 
	h = -piv*data.DWDX(d,k); 
	data.DVXDt(d,i) = data.DVXDt(d,i) + data.Mass(j)*h; 
	data.DVXDt(d,j) = data.DVXDt(d,j) - data.Mass(i)*h; 
	data.DEDt(i) = data.DEDt(i) - data.Mass(j)*dvx(d)*h; 
	data.DEDt(j) = data.DEDt(j) - data.Mass(i)*dvx(d)*h; 
      }
    }
  }
 
//  Change of specific internal energy: 
 
  for (i=1;i<=data.NTotal+data.NVirt;i++) { 
    data.DEDt(i) = 0.5e0*data.DEDt(i);
  }
 
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the average velocity to correct velocity 
//  for preventing.penetration (monaghan, 1992) 
 
//  data.NTotal : Number of particles                                  [in] 
//  data.Mass   : Particle masses                                      [in] 
//  data.NIac   : Number of interaction pairs                          [in] 
//  data.Pair_I : List of first partner of interaction pair            [in] 
//  data.Pair_J : List of second partner of interaction pair           [in] 
//  data.W      : Kernel for all interaction pairs                     [in] 
//  data.VX     : Velocity of each particle                            [in] 
//  data.Rho    : Density of each particle                             [in] 
//  data.AveVel     : Average velocityof each particle                    [out] 

void av_vel(SPH_Data &data)
{ 
  int i,j,k,d;        
  double  vcc, epsilon; 
  MkDouble dvx(data.Dim+1);

   i=j=k=d=0;        
   vcc= epsilon=0; 
       
//  epsilon --- a small constants chosen by experience, may lead to instability. 
//  for example, for the 1 dimensional shock tube problem, the E = 0.3 
 
  epsilon = 0.3;
       
  for (i = 1;i<= data.NTotal;i++ ) { 
    for (d = 1;d<= data.Dim;d++) { 
      data.AveVel(d,i) = 0.; 
    }
  }
      
  for (k=1;k<=data.NIac;k++) {        
    i = data.Pair_I(k); 
    j = data.Pair_J(k);
    for (d=1;d<=data.Dim;d++) { 
      dvx(d) = data.VX(d,i) - data.VX(d,j);
      data.AveVel(d, i) = data.AveVel(d,i) - 2*data.Mass(j)*dvx(d)/(data.Rho(i)+data.Rho(j))*data.W(k);
      data.AveVel(d, j) = data.AveVel(d,j) + 2*data.Mass(i)*dvx(d)/(data.Rho(i)+data.Rho(j))*data.W(k);                       
    }
  }         
  for (i = 1;i<= data.NTotal;i++) {
    for (d = 1;d<= data.Dim;d++){ 
      data.AveVel(d,i) = epsilon * data.AveVel(d,i); 
    }
  } 
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the density with SPH summation algorithm. 
//  See Equ.(4.35) 
 
//    data.NTotal : Number of particles                                  [in] 
//    data.Hsml   : Smoothing Length                                     [in] 
//    data.Mass   : Particle masses                                      [in] 
//    data.NIac   : Number of interaction pairs                          [in] 
//    data.Pair_I : List of first partner of interaction pair            [in] 
//    data.Pair_J : List of second partner of interaction pair           [in] 
//    data.W      : Kernel for all interaction pairs                     [in] 
//    data.IType   : type of particles                                   [in] 
//    data.X       : Coordinates of all particles                        [in] 
//    data.Rho    : Density                                             [out] 

void sum_density(SPH_Data &data) 
{
  int i, j, k, d;
  double selfdens, r;
  MkDouble hv(data.Dim+1), wi(data.MaxN+1);

   i= j= k= d=0;
   selfdens= r=0;
 
   //  MkDebug("\n      sum_density()\n");
//    wi(data.MaxN)---integration of the kernel itself 
//   MkDebug("      SD1 ");         
  for (d=1;d<=data.Dim;d++) { 
    hv(d) = 0.e0; 
  }
  //  MkDebug("SD2 ");          
//     Self density of each particle: Wii (Kernel for distance 0) 
//     and take contribution of particle itself: 
 
  r=0.;
       
//     Firstly calculate the integration of the kernel over the space 
 
  for (i=1;i<=data.NTotal+data.NVirt ;i++) {
    kernel(r,hv,data.Hsml(i),selfdens,hv,data) ;
    wi(i)=selfdens*data.Mass(i)/data.Rho(i); 
  }
  //  MkDebug("SD3 ");           
  for (k=1;k<=data.NIac;k++) { 
    i = data.Pair_I(k); 
    j = data.Pair_J(k); 
    wi(i) = wi(i) + data.Mass(j)/data.Rho(j)*data.W(k); 
    wi(j) = wi(j) + data.Mass(i)/data.Rho(i)*data.W(k); 
  }
 
//  Secondly calculate the rho integration over the space 
//  MkDebug("SD4 ");           
  for (i=1;i<=data.NTotal+data.NVirt;i++) { 
    kernel(r,hv,data.Hsml(i),selfdens,hv,data); 
    data.Rho(i) = selfdens*data.Mass(i); 
  }
  //  MkDebug("SD5 ");           
//  Calculate SPH sum for rho: 
  for (k=1;k<=data.NIac;k++) { 
    i = data.Pair_I(k); 
    j = data.Pair_J(k); 
    data.Rho(i) = data.Rho(i) + data.Mass(j)*data.W(k); 
    data.Rho(j) = data.Rho(j) + data.Mass(i)*data.W(k); 
  }
  //  MkDebug("SD6 ");           

//  Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w) 
      
  if (nor_density) {  

    for (i=1;i<= data.NTotal+data.NVirt;i++) { 
      data.Rho(i)=data.Rho(i)/wi(i);
    }
  }  
  //  MkDebug("\n      ~sum_density()\n");       
}
       
//---------------------------------------------------------------------- 
//  Subroutine to calculate the density with SPH continuity approach. 
//  See Equ.(4.34) 
 
//  data.NTotal : Number of particles                                  [in] 
//  data.Mass   : Particle masses                                      [in] 
//  data.NIac   : Number of interaction pairs                          [in] 
//  data.Pair_I : List of first partner of interaction pair            [in] 
//  data.Pair_J : List of second partner of interaction pair           [in] 
//  data.DWDX   : derivation of Kernel for all interaction pairs       [in] 
//  data.VX     : Velocities of all particles                          [in] 
//  data.IType   : type of particles                                   [in] 
//  data.X      : Coordinates of all particles                         [in] 
//  data.Rho    : Density                                              [in] 
//  data.DRhoDt : Density change rate of each particle                [out]    

void con_density(SPH_Data &data)
{ 
       
  int i,j,k,d;
  double vcc;
  MkDouble dvx(data.Dim+1);

   i=j=k=d=0;
   vcc=0;
       
  for (i = 1;i<= data.NTotal+data.NVirt;i++) { 
    data.DRhoDt(i) = 0.; 
  }
      
  for (k=1;k<=data.NIac;k++) {       
    i = data.Pair_I(k); 
    j = data.Pair_J(k); 
    for (d=1;d<=data.Dim;d++) { 
      dvx(d) = data.VX(d,i) - data.VX(d,j);
    }
    vcc = dvx(1)*data.DWDX(1,k);
    for (d=2;d<=data.Dim;d++) { 
      vcc = vcc + dvx(d)*data.DWDX(d,k); 
    }
    data.DRhoDt(i) = data.DRhoDt(i) + data.Mass(j)*vcc; 
    data.DRhoDt(j) = data.DRhoDt(j) + data.Mass(i)*vcc;
   
  }
	  
}

//---------------------------------------------------------------------- 
//  Subroutine to calculate the smoothing funciton for each particle and 
//  the interaction parameters used by the SPH algorithm. Interaction  
//  pairs are determined by directly comparing the particle distance  
//  with the corresponding smoothing length. 
//  See p.148 in Chapter 4 
 
//    data.ITimeStep : Current time step                                 [in] 
//    data.NTotal    : Number of particles                               [in] 
//    data.Hsml      : Smoothing Length                                  [in] 
//    data.X         : Coordinates of all particles                      [in] 
//    data.NIac      : Number of interaction pairs                      [out] 
//    data.Pair_I    : List of first partner of interaction pair        [out] 
//    data.Pair_J    : List of second partner of interaction pair       [out] 
//    data.W         : Kernel for all interaction pairs                 [out] 
//    data.DWDX      : Derivative of kernel with respect to x, y and z  [out] 
//    data.CountIac  : Number of neighboring particles                  [out] 

void direct_find(SPH_Data &data) 
{ 
  int i, j, d,  sumiac, maxiac, miniac, noiac, maxp, minp, scale_k;  
  double  driac, r, mhsml;      
  MkDouble dxiac(data.Dim+1), tdwdx(data.Dim+1);

   i= j= d=  sumiac= maxiac= miniac= noiac= maxp= minp= scale_k=0;  
   driac= r= mhsml=0;      

//    Smoothing kernel function  
//    skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985) 
//        = 2, Gauss kernel   (Gingold and Monaghan 1981)  
//        = 3, Quintic kernel (Morris 1997) 

//  MkDebug("    direct_find()\n");
//  MkDebug("    D1 ");
  if (skf==1) {  
    scale_k = 2;
  }  
  else if (skf==2){
    scale_k = 3 ; 
  }
  else if (skf==3) {
    scale_k = 3 ; 
  }
      
  for (i=1;i<=data.NTotal+data.NVirt;i++){
    data.CountIac(i) = 0; 
  }
  //  MkDebug("D2 ");
       
  data.NIac = 0;
 
  for (i=1;i<=data.NTotal+data.NVirt-1;i++) {      
    for (j = i+1;j<= data.NTotal+data.NVirt;j++){ 
      dxiac(1) = data.X(1,i) - data.X(1,j); 
      driac    = dxiac(1)*dxiac(1); 

      for (d=2;d<=data.Dim;d++) { 
	dxiac(d) = data.X(d,i) - data.X(d,j); 
	driac    = driac + dxiac(d)*dxiac(d); 
      }

      mhsml = (data.Hsml(i)+data.Hsml(j))/2.; 
      if (sqrt(driac)<scale_k*mhsml) { 
	if (data.NIac<max_interaction) {     
 
//    Neighboring pair list, and totalinteraction number and 
//    the interaction number for each particle  

	  data.NIac = data.NIac + 1; 
	  data.Pair_I(data.NIac) = i; 
	  data.Pair_J(data.NIac) = j; 
	  r = sqrt(driac); 
	  data.CountIac(i) = data.CountIac(i) + 1; 
	  data.CountIac(j) = data.CountIac(j) + 1; 

//    Kernel and derivations of kernel 
 
	  kernel(r,dxiac,mhsml,data.W(data.NIac),tdwdx,data) ;

	  for (d=1;d<=data.Dim;d++) { 
	    data.DWDX(d,data.NIac) = tdwdx(d);
	  }

	}
	else {
	  printf(" >>> ERROR << : Too many interactions\n");  
	  exit(-1);//stop 
	}
      } 
    }
  }
  //  MkDebug("D3 "); 
//    Statistics for the interaction 
 
  sumiac = 0; 
  maxiac = 0; 
  miniac = 1000; 
  noiac  = 0; 
  for (i=1;i<=data.NTotal+data.NVirt;i++) { 
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
  
  if ((data.ITimeStep%print_step)==0) {
    if (int_stat) {
      printf(" >> Statistics: interactions per particle:\n"); 
      printf("**** Particle: %d  maximal interactions: %d\n",maxp,maxiac) ;
      printf("**** Particle: %d minimal interactions: %d\n",minp,miniac) ;
      printf("**** Average :%f\n",float(sumiac)/float(data.NTotal+data.NVirt)) ;
      printf("**** Total pairs : %d\n",data.NIac) ;
      printf("**** Particles with no interactions:%d\n",noiac) ;
    }
  }
  //  MkDebug("D4 "); 
  //  MkDebug("    ~direct_find()\n");
}


//---------------------------------------------------------------------- 
//  Gamma law EOS: subroutine to calculate the pressure and sound   
  
//  rho    : Density                                              [in] 
//  u      : Internal energy                                      [in] 
//  p      : Pressure                                            [out] 
//  c      : sound velocity                                      [out] 

void p_gas(double rho, double u, double &p, double &c) 
{       
           
  double gamma  ;
           
//   For air (idea gas) 
//   See Equ.(3.82) 
 
  gamma=1.4; 
  p = (gamma-1) * rho * u;      
  c = sqrt((gamma-1) * u);  
      
}

//---------------------------------------------------------------------- 
//   Artificial equation of state for the artificial compressibility  
 
//  rho    : Density                                              [in] 
//  u      : Internal energy                                      [in] 
//  p      : Pressure                                            [out] 
//  c      : sound velocity                                      [out] 
//  Equation of state for artificial compressibility    
       
void p_art_water(double rho, double &p, double &c) 
{       
  double gamma, rho0; 
 
//  Artificial EOS, Form 1 (Monaghan, 1994)  
//  See Equ.(4.88) 
//   gamma=7. 
//   rho0=1000.        
//   b = 1.013e5 
//   p = b*((rho/rho0)**gamma-1)       
//   c = 1480. 
 
//  Artificial EOS, Form 2 (Morris, 1997) 
//  See Equ.(4.89) 
  c = 0.01; 
  p = c*c * rho;
}

//-------------------------------------------------------------------------- 
//  Subroutine to calculate the external forces, e.g. gravitational forces.       
//  The forces from the interactions with boundary virtual particles  
//  are also calculated here as external forces. 
 
//  here as the external force.  
//  data.NTotal  : Number of particles                                 [in] 
//  data.Mass    : Particle masses                                     [in] 
//  data.X       : Coordinates of all particles                        [in] 
//  data.Pair_I : List of first partner of interaction pair            [in] 
//  data.Pair_J : List of second partner of interaction pair           [in] 
//  data.IType   : type of particles                                   [in] 
//  data.Hsml   : Smoothing Length                                     [in] 
//  data.DVXDt   : Acceleration with respect to x, y and z            [out]  

void ext_force(SPH_Data &data) 
{ 

  int i, j, k, d ;
  double rr, f, rr0, dd, p1, p2;
  MkDouble dx(data.Dim+1);

   i= j= k= d =0;
   rr= f= rr0= dd= p1= p2=0;
            
  for (i = 1;i<= data.NTotal+data.NVirt;i++) { 
    for (d = 1;d<= data.Dim;d++){ 
      data.DVXDt(d, i) = 0.; 
    }
  }
         
//  Consider self-gravity or not ? 
 
  if (self_gravity) { 
    for (i = 1; i<=data.NTotal+data.NVirt;i++) { 
      data.DVXDt(data.Dim, i) = -9.8; 
    }
  }  
 
//  Boundary particle force and penalty anti-penetration force.  
  rr0 = 1.25e-5; 
  dd = 1.e-2; 
  p1 = 12; 
  p2 = 4; 
       
  for ( k=1;k<=data.NIac;k++){ 
    i = data.Pair_I(k); 
    j = data.Pair_J(k);   
    if(data.IType(i)>0&&data.IType(j)<0) {   
      rr = 0.;       
      for (d=1;d<=data.Dim;d++){ 
	dx(d) =  data.X(d,i) -  data.X(d,j); 
	rr = rr + dx(d)*dx(d); 
      }
      rr = sqrt(rr); 
      if(rr<rr0) { 
	f = (pow(rr0/rr,p1)-pow(rr0/rr,p2))/rr*rr ;
	for (d = 1; d<=data.Dim;d++) { 
	  data.DVXDt(d, i) = data.DVXDt(d, i) + dd*dx(d)*f ;
	}
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
 
void grid_geom(int i,MkDouble x,MkInt ngridx,MkDouble maxgridx,MkDouble mingridx,MkDouble dgeomx,MkInt &xgcell,SPH_Data &data) 
{
  int d; 
 
  for (d=1;d<=3;d++) { 
    xgcell(d) = 1; 
  }
 
  for (d=1;d<=data.Dim;d++) { 
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
 
//  data.Dt     : time step                                            [in] 
//  data.NTotal : Number of particles                                  [in] 
//  data.Mass   : Particle masses                                      [in] 
//  data.VX     : Velocities of all particles                          [in] 
//  data.Rho    : Density                                              [in] 
//  data.NIac   : Number of interaction pairs                          [in] 
//  data.Pair_I : List of first partner of interaction pair            [in] 
//  data.Pair_J : List of second partner of interaction pair           [in] 
//  data.DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//  data.Hsml   : Smoothing Length                                 [in/out] 

void h_upgrade(SPH_Data &data) 
{    
  int i,j,k,d; 
  double fac, hvcc;
  MkDouble dvx(data.Dim+1), vcc(data.MaxN+1), dhsml(data.MaxN+1);

   i=j=k=d=0; 
   fac= hvcc=0;
 
  if (sle==0 ) {      
 //---  Keep smoothing length unchanged.  
     return; 
  }      
  else if (sle==2) {
       
//---  dh/dt = (-1/data.Dim)*(h/rho)*(drho/dt). 
 
    for (i=1;i<=data.NTotal ;i++) {
      vcc(i) = 0.e0 ;
    }
       
    for (k=1;k<=data.NIac;k++) {
      i = data.Pair_I(k); 
      j = data.Pair_J(k) ;
      for (d=1;d<=data.Dim;d++) {
	dvx(d) = data.VX(d,j) - data.VX(d,i) ;
      }
      hvcc = dvx(1)*data.DWDX(1,k);
      for (d=2;d<=data.Dim;d++) {
	hvcc = hvcc + dvx(d)*data.DWDX(d,k);
      }
      vcc(i) = vcc(i) + data.Mass(j)*hvcc/data.Rho(j); 
      vcc(j) = vcc(j) + data.Mass(i)*hvcc/data.Rho(i);
    }
         
    for (i = 1;i<= data.NTotal;i++) {
      dhsml(i) = (data.Hsml(i)/data.Dim)*vcc(i);
      data.Hsml(i) = data.Hsml(i) + data.Dt*dhsml(i);
      if (data.Hsml(i)<=0) data.Hsml(i) = data.Hsml(i) - data.Dt*dhsml(i);  
    }
  }
     
  else if(sle==1) {
    fac = 2.0;
    for (i = 1;i<= data.NTotal;i++) {
      data.Hsml(i) = fac * (data.Mass(i)/pow(data.Rho(i),(1./data.Dim)));
    }
        
  }
}

//----------------------------------------------------------------------       
//  Subroutine to established a pair linked list by sorting grid cell. 
//  It is suitable for a homogeneous particle distribution with the  
//  same smoothing length in an instant. A fixed number of particles 
//  lie in each cell.  
 
//    data.NTotal   : Number of particles                                [in] 
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

void init_grid(SPH_Data &data,SPH_Grid &grid) 
{
  int i, j, k, d;
  MkDouble maxng(data.Dim+1), ngrid(3+1);
  double nppg=0; 

   i= j= k= d=0;
 
//    Averaged number of particles per grid cell 
 
  nppg = 3.e0;
 
//    Initialize parameters: Maximum number of grid cells 
 
  maxng(1) = grid.MaxNGridX; 
  if (data.Dim>=2) {
    maxng(2) = grid.MaxNGridY;
    if (data.Dim==3) { 
      maxng(3) = grid.MaxNGridZ; 
    }
  }
       
  for (d=1;d<=3;d++) {
    ngrid(d) = 1;
  }
       
//    Range of sorting grid 
 
  grid.MaxGridX(1) = x_maxgeom; 
  grid.MinGridX(1) = x_mingeom; 
  if (data.Dim>=2) { 
    grid.MaxGridX(2) = y_maxgeom; 
    grid.MinGridX(2) = y_mingeom; 
    if (data.Dim==3) { 
      grid.MaxGridX(3) = z_maxgeom;
      grid.MinGridX(3) = z_mingeom;
    }
  }
 
  for (d=1;d<=data.Dim;d++) {
    grid.DGeomX(d) = grid.MaxGridX(d) - grid.MinGridX(d);
  }
 
//    Number of grid cells in x-, y- and z-direction: 
 
  if (data.Dim==1) {
    grid.NGridX(1) = min(int(data.NTotal/nppg) + 1,maxng(1));
  }
  else if (data.Dim==2) {
    grid.NGridX(1) = min(int(sqrt(data.NTotal*grid.DGeomX(1)/(grid.DGeomX(2)*nppg))) + 1,maxng(1)); 
    grid.NGridX(2) = min(int(grid.NGridX(1)*grid.DGeomX(2)/grid.DGeomX(1)) + 1,maxng(2));
  } 
  else if (data.Dim==3) {
    grid.NGridX(1) = min(int(pow(data.NTotal*grid.DGeomX(1)*grid.DGeomX(1)/(grid.DGeomX(2)*grid.DGeomX(3)*nppg),(1.e0/3.e0))) + 1,maxng(1)); 
    grid.NGridX(2) = min(int(grid.NGridX(1)*grid.DGeomX(2)/grid.DGeomX(1)) + 1,maxng(2)); 
    grid.NGridX(3) = min(int(grid.NGridX(1)*grid.DGeomX(3)/grid.DGeomX(1)) + 1,maxng(3)); 
  }
 
//    Smoothing Length measured in grid cells: 
 
  for (d=1;d<=data.Dim;d++) {
    grid.GHsmlX(d) = int(float(grid.NGridX(d))*data.Hsml(1)/grid.DGeomX(d)) + 1 ;
  }
 
  for (d=1;d<=data.Dim;d++) {
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
 
void input(SPH_Data &data) 
{
  int i, d, im;
  FILE *fp1,*fp2, *fp3;
 
//   load initial particle information from external disk file 
 
//  MkDebug("input():");
  if(config_input){
 
    //    MkDebug("config \n");                        
    fp1 = fopen("f_xv.dat","r") ;
    fp2 = fopen("f_state.dat","r");
    fp3 = fopen("f_other.dat","r");
       
    printf("  **************************************************\n"); 
    printf("      Loading initial particle configuration...   \n");
    fscanf (fp1,"%d", &data.NTotal);
    printf("      Total number of particles  %f\n ", data.NTotal)    	 ;
    printf("  **************************************************\n")	 ;
    for (i = 1;i<= data.NTotal;i++) {
      fscanf(fp1,"%d",&im);
      for (d=1;d<=data.Dim;d++) fscanf(fp1,"%f",&data.X(d, i));
      for (d=1;d<=data.Dim;d++) fscanf(fp1,"%f",&data.VX(d, i));
      fscanf(fp2,"%d %f %f %f %f",&im, &data.Mass(i), &data.Rho(i), &data.P(i), &data.U(i));         
      fscanf(fp3,"%d %d %f",&im, &data.IType(i), &data.Hsml(i));
    }
  }      
  else  {
           
    //    MkDebug("not config \n");                        
    fp1 = fopen("ini_xv.dat","w") ;
    fp2 = fopen("ini_state.dat","w") ;
    fp3 = fopen("ini_other.dat","w")  ;

    if (shocktube) shock_tube(data);
    if (shearcavity) shear_cavity(data);

    for (i = 1;i<= data.NTotal ;i++) {
      fprintf(fp1,"%6d ", i);
      for (d=1;d<=data.Dim;d++) fprintf(fp1,"%14.8f ", data.X(d, i));
      for (d=1;d<=data.Dim;d++) fprintf(fp1,"%14.8f ", data.VX(d, i));
      fprintf(fp1,"\n");
      fprintf(fp2,"%6d %14.8f %14.8f %14.8f %14.8f \n", i, data.Mass(i), data.Rho(i), data.P(i), data.U(i));
      fprintf(fp3,"%6d %3d %14.8f \n", i, data.IType(i), data.Hsml(i));
    }
    printf("  **************************************************\n");
    printf("      Initial particle configuration generated   \n");       
    printf("      Total number of particles %d\n ", data.NTotal);	 
    printf("  **************************************************\n");
  }
 
  fclose(fp1) ;
  fclose(fp2) ;
  fclose(fp3)  ;
  //  MkDebug("~input():\n");
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
        
void shock_tube(SPH_Data &data) 
{
  int i, d; 
  double space_x;
 
  data.NTotal=400; 
  space_x=0.6/80.;       

  //  MkDebug("  shock_tube() : \n");

  for (i=1;i<=data.NTotal;i++) {
    data.Mass(i)=0.75/400.; 
    data.Hsml(i)=0.015; 
    data.IType(i)=1;
    for (d = 1;d<= data.Dim;d++) {
      data.X(d,i) = 0.;
      data.VX(d,i) = 0. ;
    }
  }

  for (i=1;i<=320;i++) {
    data.X(1,i)=-0.6+space_x/4.*(i-1); 
  }
       
  for (i=320+1;i<=data.NTotal;i++) {
    data.X(1,i)=0.+space_x*(i-320);
  }         

  for (i=1;i<=data.NTotal;i++) {
    if (data.X(1,i)<=1.e-8)  {
	data.U(i)=2.5; 
	data.Rho(i)=1. ;
	data.P(i)=1. ;
      }
    if (data.X(1,i)>1.e-8){
      data.U(i)=1.795; 
      data.Rho(i)=0.25 ;
      data.P(i)=0.1795 ;
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
//   data.NTotal-- total particle number                               [out] 
void shear_cavity(SPH_Data &data) 
{
  int i, j, d, m, n, mp, np, k; 
  double xl, yl, dx, dy; 
 
//   Giving data.Mass and smoothing length as well as other data. 
 
//  MkDebug("shear_cavity()\n");
  m = 41; 
  n = 41; 
  mp = m-1; 
  np = n-1; 
  data.NTotal = mp * np; 
  xl = 1.e-3; 
  yl = 1.e-3; 
  dx = xl/mp; 
  dy = yl/np; 
 
  for (i = 1;i<=mp;i++) {
    for (j = 1;j<= np;j++) {
      k = j + (i-1)*np; 
      data.X(1, k) = (i-1)*dx + dx/2.; 
      data.X(2, k) = (j-1)*dy + dy/2. ;
    }
  }
 
  for (i = 1;i<= mp*np ;i++) {
    data.VX(1, i) = 0.; 
    data.VX(2, i) = 0.;       
    data.Rho (i) = 1000.;    
    data.Mass(i) = dx*dy*data.Rho(i);   
    data.P(i)= 0.;    
    data.U(i)=357.1; 
    data.IType(i) = 2; 
    data.Hsml(i) = dx; 
  }
  //  MkDebug("~shear_cavity()\n");
}





//---------------------------------------------------------------------- 
//  Subroutine to calculate the internal forces on the right hand side  
//  of the Navier-Stokes equations, i.e. the pressure gradient and the 
//  gradient of the viscous stress tensor, used by the time integration.  
//  Moreover the entropy production due to viscous dissipation, tds/dt,  
//  and the change of internal energy per mass, de/dt, are calculated.  

//    data.ITimeStep: Current timestep number                            [in] 
//    data.Dt     :   Time step                                          [in] 
//    data.NTotal : Number of particles                                  [in] 
//    data.Hsml   : Smoothing Length                                     [in] 
//    data.Mass   : Particle masses                                      [in] 
//    data.VX     : Velocities of all particles                          [in] 
//    nia//  : Number of interaction pairs                          [in] 
//    data.Rho    : Density                                              [in] 
//    data.Eta    : Dynamic viscosity                                    [in] 
//    data.Pair_I : List of first partner of interaction pair            [in] 
//    data.Pair_J : List of second partner of interaction pair           [in] 
//    data.DWDX   : Derivative of kernel with respect to x, y and z      [in] 
//    data.IType  : Type of particle (material types)                    [in] 
//    data.U      : Particle internal energy                             [in] 
//    data.X      : Particle coordinates                                 [in] 
//    t      : Particle temperature                             [in/out] 
//   data.C      : Particle sound speed                                [out] 
//   p      : Particle pressure                                   [out] 
//   data.DVXDt  : Acceleration with respect to x, y and z             [out]  
//   data.TDSDt  : Production of viscous entropy                       [out] 
//   data.DEDt   : Change of specific internal energy                  [out] 


void int_force(SPH_Data &data) 
{
  int i, j, k, d;
  double hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij;
  MkDouble dvx(data.Dim+1), txx(data.MaxN+1), tyy(data.MaxN+1),tzz(data.MaxN+1), txy(data.MaxN+1), txz(data.MaxN+1), tyz(data.MaxN+1), vcc(data.MaxN+1);

//   Initialization of shear tensor, velocity divergence,  
//   viscous energy, internal energy, acceleration  

   i= j= k= d=0;
   hxx= hyy= hzz= hxy= hxz= hyz= h= hvcc= he= rhoij=0;

  for (i=1;i<=data.NTotal+data.NVirt;i++) {
    txx(i) = 0.e0; 
    tyy(i) = 0.e0; 
    tzz(i) = 0.e0; 
    txy(i) = 0.e0; 
    txz(i) = 0.e0; 
    tyz(i) = 0.e0; 
    vcc(i) = 0.e0; 
    data.TDSDt(i) = 0.e0; 
    data.DEDt(i) = 0.e0; 
    for (d=1;d<=data.Dim;d++) { 
      data.DVXDt(d,i) = 0.e0;
    }
  }
 
//   Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c 
 
  if (visc) {
    for (k=1;k<=data.NIac;k++) { 
      i = data.Pair_I(k); 
      j = data.Pair_J(k); 
      for (d=1;d<=data.Dim;d++) {
	dvx(d) = data.VX(d,j) - data.VX(d,i);
      }
      if (data.Dim==1) {
	hxx = 2.e0*dvx(1)*data.DWDX(1,k);
      }
      else if (data.Dim==2) {
	hxx = 2.e0*dvx(1)*data.DWDX(1,k) -  dvx(2)*data.DWDX(2,k);
	hxy = dvx(1)*data.DWDX(2,k) + dvx(2)*data.DWDX(1,k);
	hyy = 2.e0*dvx(2)*data.DWDX(2,k) - dvx(1)*data.DWDX(1,k);
      }
      else if (data.Dim==3){
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
      if (data.Dim==1) {
	txx(i) = txx(i) + data.Mass(j)*hxx/data.Rho(j);
	txx(j) = txx(j) + data.Mass(i)*hxx/data.Rho(i);
      }
      else if (data.Dim==2) {
	txx(i) = txx(i) + data.Mass(j)*hxx/data.Rho(j);
	txx(j) = txx(j) + data.Mass(i)*hxx/data.Rho(i);
	txy(i) = txy(i) + data.Mass(j)*hxy/data.Rho(j);
	txy(j) = txy(j) + data.Mass(i)*hxy/data.Rho(i);             
	tyy(i) = tyy(i) + data.Mass(j)*hyy/data.Rho(j); 
	tyy(j) = tyy(j) + data.Mass(i)*hyy/data.Rho(i);
      }
      else if (data.Dim==3) {
	txx(i) = txx(i) + data.Mass(j)*hxx/data.Rho(j);
	txx(j) = txx(j) + data.Mass(i)*hxx/data.Rho(i);    
	txy(i) = txy(i) + data.Mass(j)*hxy/data.Rho(j);
	txy(j) = txy(j) + data.Mass(i)*hxy/data.Rho(i);
	txz(i) = txz(i) + data.Mass(j)*hxz/data.Rho(j);
	txz(j) = txz(j) + data.Mass(i)*hxz/data.Rho(i);
	tyy(i) = tyy(i) + data.Mass(j)*hyy/data.Rho(j);
	tyy(j) = tyy(j) + data.Mass(i)*hyy/data.Rho(i);
	tyz(i) = tyz(i) + data.Mass(j)*hyz/data.Rho(j);
	tyz(j) = tyz(j) + data.Mass(i)*hyz/data.Rho(i);
	tzz(i) = tzz(i) + data.Mass(j)*hzz/data.Rho(j);
	tzz(j) = tzz(j) + data.Mass(i)*hzz/data.Rho(i);
      }
 
//   Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz: 
 
      hvcc = 0.;
      for (d=1;d<=data.Dim;d++) {
	hvcc = hvcc + dvx(d)*data.DWDX(d,k);
      }
      vcc(i) = vcc(i) + data.Mass(j)*hvcc/data.Rho(j);
      vcc(j) = vcc(j) + data.Mass(i)*hvcc/data.Rho(i);
    }
  }
 
  for (i=1;i<=data.NTotal+data.NVirt;i++) {
//   Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab 
    if (visc) {
      if (data.Dim==1) {
	data.TDSDt(i) = txx(i)*txx(i);
      }
      else if (data.Dim==2) {
	data.TDSDt(i) = txx(i)*txx(i) + 2.e0*txy(i)*txy(i)+ tyy(i)*tyy(i);
      }
      else if (data.Dim==3) {
	data.TDSDt(i) = txx(i)*txx(i)+ tyy(i)*tyy(i)+ tzz(i)*tzz(i)
                      + 2.e0*txy(i)*txy(i)+ 2.e0*txz(i)*txz(i)+ 2.e0*tyz(i)*tyz(i)  ;
      }
      data.TDSDt(i) = 0.5e0*data.Eta(i)/data.Rho(i)*data.TDSDt(i);
    }
 
//   Pressure from equation of state 
 
    if (abs(data.IType(i))==1) {
      p_gas(data.Rho(i), data.U(i), data.P(i),data.C(i));
    }
    else if (abs(data.IType(i))==2) {
      p_art_water(data.Rho(i), data.P(i), data.C(i));
    }
  } 
//    Calculate SPH sum for pressure force -p,a/rho 
//    and viscous force (eta Tab),b/rho 
//    and the internal energy change de/dt due to -p/rho vc,c 
 
  for (k=1;k<=data.NIac;k++) {
    i = data.Pair_I(k);
    j = data.Pair_J(k);
    he = 0.e0;
       
//   For SPH algorithm 1 
 
    rhoij = 1.e0/(data.Rho(i)*data.Rho(j));
    if(pa_sph==1){
      for (d=1;d<=data.Dim ;d++) {
         
//   Pressure part 
                     
	h = -(data.P(i) + data.P(j))*data.DWDX(d,k);
	he = he + (data.VX(d,j) - data.VX(d,i))*h;
 
//   Viscous force 
 
	if (visc) {
	  if (d==1) {
//   x-coordinate of acceleration 
	    h = h + (data.Eta(i)*txx(i) + data.Eta(j)*txx(j))*data.DWDX(1,k);
	    if (data.Dim>=2) {
	      h = h + (data.Eta(i)*txy(i) + data.Eta(j)*txy(j))*data.DWDX(2,k);
	      if (data.Dim==3) {
		h = h + (data.Eta(i)*txz(i) + data.Eta(j)*txz(j))*data.DWDX(3,k);
	      }
	    }
	  }
	  else if (d==2) {
//   y-coordinate of acceleration 
	    h = h + (data.Eta(i)*txy(i) + data.Eta(j)*txy(j))*data.DWDX(1,k) + (data.Eta(i)*tyy(i) + data.Eta(j)*tyy(j))*data.DWDX(2,k);
	    if (data.Dim==3) {
	      h = h + (data.Eta(i)*tyz(i) + data.Eta(j)*tyz(j))*data.DWDX(3,k);
	    }
	  }
	  else if (d==3) {
//   z-coordinate of acceleration 
	    h = h + (data.Eta(i)*txz(i) + data.Eta(j)*txz(j))*data.DWDX(1,k) + (data.Eta(i)*tyz(i) + data.Eta(j)*tyz(j))*data.DWDX(2,k)+ (data.Eta(i)*tzz(i) + data.Eta(j)*tzz(j))*data.DWDX(3,k);
	  }
	}
	h = h*rhoij;
	data.DVXDt(d,i) = data.DVXDt(d,i) + data.Mass(j)*h;
	data.DVXDt(d,j) = data.DVXDt(d,j) - data.Mass(i)*h;
      }
      he = he*rhoij;
      data.DEDt(i) = data.DEDt(i) + data.Mass(j)*he;
      data.DEDt(j) = data.DEDt(j) + data.Mass(i)*he;
    }
//   For SPH algorithm 2 
           
    else if (pa_sph==2){
      for (d=1;d<=data.Dim;d++) {
	h = -(data.P(i)/data.Rho(i)/data.Rho(i) + data.P(j)/data.Rho(j)/data.Rho(j))*data.DWDX(d,k);
	he = he + (data.VX(d,j) - data.VX(d,i))*h;
 //   Viscous force 
	if (visc) {
	  if (d==1) {
//   x-coordinate of acceleration 
	    h = h + (data.Eta(i)*txx(i)/data.Rho(i)/data.Rho(i) + data.Eta(j)*txx(j)/data.Rho(j)/data.Rho(j))*data.DWDX(1,k);
	    if (data.Dim>=2){
	      h = h + (data.Eta(i)*txy(i)/data.Rho(i)/data.Rho(i) + data.Eta(j)*txy(j)/data.Rho(j)/data.Rho(j))*data.DWDX(2,k);
	      if (data.Dim==3) {
		h = h + (data.Eta(i)*txz(i)/data.Rho(i)/data.Rho(i) + data.Eta(j)*txz(j)/data.Rho(j)/data.Rho(j))*data.DWDX(3,k);
	      }
	    }
	  }             
	  else if (d==2) {
//   y-coordinate of acceleration 
	    h = h + (data.Eta(i)*txy(i)/data.Rho(i)/data.Rho(i)+data.Eta(j)*txy(j)/data.Rho(j)/data.Rho(j))*data.DWDX(1,k)
	          + (data.Eta(i)*tyy(i)/data.Rho(i)/data.Rho(i)+data.Eta(j)*tyy(j)/data.Rho(j)/data.Rho(j))*data.DWDX(2,k);
	    if (data.Dim==3) {
	      h = h + (data.Eta(i)*tyz(i)/data.Rho(i)/data.Rho(i)+data.Eta(j)*tyz(j)/data.Rho(j)/data.Rho(j))*data.DWDX(3,k);
	    }
	  }
	  else if (d==3) {
//   z-coordinate of acceleration 
	    h = h + (data.Eta(i)*txz(i)/data.Rho(i)/data.Rho(i) +data.Eta(j)*txz(j)/data.Rho(j)/data.Rho(j))*data.DWDX(1,k) 
	          + (data.Eta(i)*tyz(i)/data.Rho(i)/data.Rho(i) +data.Eta(j)*tyz(j)/data.Rho(j)/data.Rho(j))*data.DWDX(2,k)
	          + (data.Eta(i)*tzz(i)/data.Rho(i)/data.Rho(i) +data.Eta(j)*tzz(j)/data.Rho(j)/data.Rho(j))*data.DWDX(3,k);
	  }
	}
	data.DVXDt(d,i) = data.DVXDt(d,i) + data.Mass(j)*h;
	data.DVXDt(d,j) = data.DVXDt(d,j) - data.Mass(i)*h;
      }
      data.DEDt(i) = data.DEDt(i) + data.Mass(j)*he;
      data.DEDt(j) = data.DEDt(j) + data.Mass(i)*he;
    }
  }
 
//   Change of specific internal energy de/dt = T ds/dt - p/rho vc,c: 
 
  for (i=1;i<=data.NTotal+data.NVirt;i++) {
    data.DEDt(i) = data.TDSDt(i) + 0.5e0*data.DEDt(i);
  }
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
 


void kernel(double r,MkDouble dx,double hsml,double &w,MkDouble &dwdx,SPH_Data &data)    
{       
  int i, j, d;
  double q, dw, factor;

   i= j= d=0;
   q= dw= factor=0;
 
  q = r/hsml;
  w = 0.e0;
  for (d=1;d<=data.Dim;d++) {
    dwdx(d) = 0.e0 ;
  }
 
  if (skf==1) {
       
    if (data.Dim==1) factor = 1.e0/hsml;
    else if (data.Dim==2) factor = 15.e0/(7.e0*pi*hsml*hsml);
    else if (data.Dim==3) factor = 3.e0/(2.e0*pi*hsml*hsml*hsml);
    else  {
      printf(" >>> Error << : Wrong dimension: Dim =%f\n",data.Dim);
      exit(0);//stop 
    }
       
    if (q>=0&&q<=1.e0) {
      w = factor * (2./3. - q*q + q*q*q / 2.);
      for (d = 1;d<= data.Dim;d++) {
	dwdx(d) = factor * (-2.+3./2.*q)/hsml/hsml * dx(d);
      }
    }
    else if (q>1.e0&&q<=2) {
      w = factor * 1.e0/6.e0 * (2.-q)*(2.-q)*(2.-q);
      for (d = 1;d<= data.Dim;d++) {
	dwdx(d) =-factor * 1.e0/6.e0 * 3.*(2.-q)*(2.-q)/hsml * (dx(d)/r);
      }
    }
    else {
      w=0.;
      for (d= 1;d<= data.Dim;d++) {
	dwdx(d) = 0.;
      }
    }
  }

  else if (skf==2) {
       
    factor = 1.e0 / (pow(hsml,data.Dim) * pow(pi,(data.Dim/2.)));
    if(q>=0&&q<=3) {
      w = factor * exp(-q*q);
      for (d = 1;d<= data.Dim;d++) {
	dwdx(d) = w * ( -2.* dx(d)/hsml/hsml);
      }
    }
    else {
      w = 0.;
      for (d = 1;d<= data.Dim;d++) {
	dwdx(d) = 0.;
      }	    
    }
  }
	 
  else if (skf==3) {
    if (data.Dim==1) factor = 1.e0 / (120.e0*hsml);
    else if (data.Dim==2) factor = 7.e0 / (478.e0*pi*hsml*hsml);
    else if (data.Dim==3) factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml) ;
    else {
      printf(" >>> Error << : Wrong dimension: Dim =%f\n",data.Dim);
      exit(-1);//stop 
    }
             
    if(q>=0&&q<=1) {
      w = factor * ( pow(3-q,5) - 6*pow(2-q,5) + 15*pow(1-q,5) );
      for (d= 1;d<= data.Dim;d++) {
	dwdx(d) = factor * ( (-120 + 120*q - 50*q*q)/ hsml/ hsml * dx(d) ); 
      }
    }
    else if(q>1&&q<=2) {
      w = factor * ( pow(3-q,5) - 6*pow(2-q,5) );
      for (d= 1;d<= data.Dim;d++) {
	dwdx(d) = factor * (-5*pow((3-q),4) + 30*pow(2-q,4))/ hsml * (dx(d)/r);
      }
    }
    else if(q>2&&q<=3) {
      w = factor * pow(3-q,5);
      for (d= 1;d<= data.Dim;d++) {
	dwdx(d) = factor * (-5*pow(3-q,4)) / hsml * (dx(d)/r);
      }
    }
    else {
      w = 0.;
      for (d = 1;d<= data.Dim;d++) {
	dwdx(d) = 0.;
      }
    }
  }		 
}

//---------------------------------------------------------------------- 
// Subroutine to calculate the smoothing funciton for each particle and 
// the interaction parameters used by the SPH algorithm. Interaction  
// pairs are determined by using a sorting grid linked list   
 
//   data.ITimeStep : Current time step                                 [in] 
//   data.NTotal    : Number of particles                               [in] 
//   data.Hsml      : Smoothing Length, same for all particles          [in] 
//   data.X         : Coordinates of all particles                      [in] 
//   data.NIac      : Number of interaction pairs                      [out] 
//   data.Pair_I    : List of first partner of interaction pair        [out] 
//   data.Pair_J    : List of second partner of interaction pair       [out] 
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
void link_list(SPH_Data &data) 
{
  int i, j, d, scale_k, sumiac, maxiac, noiac, miniac, maxp,minp,xcell,ycell,zcell;
  double hsml2,dr,r;

  MkInt gcell(3+1), celldata(data.MaxN+1),minxcell(3+1),maxxcell(3+1), dnxgcell(data.Dim+1),dpxgcell(data.Dim+1);
  MkDouble dx(data.Dim+1),tdwdx(data.Dim+1), dgeomx(data.Dim+1);
  SPH_Grid grid; 

  int maxngx = 100, maxngy = 100, maxngz = 1;

   i= j= d= scale_k= sumiac= maxiac= noiac= miniac= maxp=minp=xcell=ycell=zcell=0;
   hsml2=dr=r=0;

  grid.Initialize(data.Dim,maxngx,maxngy,maxngz);

  if (skf==1) scale_k = 2;
  else if (skf==2) scale_k = 3; 
  else if (skf==3) scale_k = 3;
      
  for (i=1;i<=data.NTotal+data.NVirt;i++) {
    data.CountIac(i) = 0;
  }
 
//   Initialize grid:   
 
  init_grid(data,grid); 
       
//   Position particles on grid and create linked list: 
       
  for (i=1;i<=data.NTotal+data.NVirt;i++) {
    grid_geom(i,data.X(1,i),grid.NGridX,grid.MaxGridX,grid.MinGridX,grid.DGeomX,gcell,data);
    for (d=1;d<=data.Dim;d++) {
      data.XGCell(d,i) = gcell(d);
    }
    celldata(i) = grid.Grid(gcell(1),gcell(2),gcell(3));
    grid.Grid(gcell(1),gcell(2),gcell(3)) = i; 
  }
 
//   Determine interaction parameters: 
 
  data.NIac = 0;
  for (i=1;i<=data.NTotal+data.NVirt-1;i++) {
 //   Determine range of grid to go through: 
    for (d=1;d<=3;d++) {
      minxcell(d) = 1;
      maxxcell(d) = 1;
    }
    for (d=1;d<=data.Dim;d++) {
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
	    dx(1) = data.X(1,i) - data.X(1,j);
	    dr    = dx(1)*dx(1); 
	    for (d=2;d<=data.Dim;d++) {
	      dx(d) = data.X(d,i) - data.X(d,j);
	      dr    = dr + dx(d)*dx(d);
	    }
	    if (sqrt(dr)<scale_k*data.Hsml(1)) {
	      if (data.NIac<max_interaction) {
 
//   Neighboring pair list, and totalinteraction number and 
//   the interaction number for each particle  
 
		data.NIac = data.NIac + 1; 
		data.Pair_I(data.NIac) = i; 
		data.Pair_J(data.NIac) = j; 
		r = sqrt(dr); 
		data.CountIac(i) = data.CountIac(i) + 1; 
		data.CountIac(j) = data.CountIac(j) + 1; 
                            
//--- Kernel and derivations of kernel 
 
		kernel(r,dx,data.Hsml(1),data.W(data.NIac),tdwdx,data); 
		for (d = 1;d<= data.Dim;d++) {
		  data.DWDX(d,data.NIac)=tdwdx(d);
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
  for (i=1;i<=data.NTotal+data.NVirt;i++) {
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
  
  if ((data.ITimeStep%print_step)==0) {
    if (int_stat) {
      printf(" >> Statistics: interactions per particle:\n");
      printf("**** Particle:%f, maximal interactions:%f\n",maxp, maxiac); 
      printf("**** Particle:%f, minimal interactions:%f\n",minp, miniac); 
      printf("**** Average :%f\n",float(sumiac)/float(data.NTotal+data.NVirt)); 
      printf("**** Total pairs : %f",data.NIac); 
      printf("**** Particles with no interactions: %f",noiac); 
    }
  }
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
//   data.NTotal-- total particle number                                [in] 
 
void output(SPH_Data &data)  
{
  int i, d;
  FILE *fp1,*fp2, *fp3;
       
  fp1 = fopen("f_xv.dat","w"); 
  fp2 = fopen("f_state.dat","w");
  fp3 = fopen("f_other.dat","w");
      
  fprintf(fp1,"%d\n", data.NTotal); 
  for (i = 1;i<= data.NTotal;i++) {
    fprintf(fp1,"%6d ",i);
    for(d=1;d<=data.Dim;d++) fprintf(fp1,"%14.8f ", (data.X(d, i)));
    for(d=1;d<=data.Dim;d++) fprintf(fp1,"%14.8f ",(data.VX(d, i)));
    fprintf(fp1,"\n ");
    fprintf(fp2,"%6d %14.8f %14.8f %14.8f %14.8f \n", i, data.Mass(i), data.Rho(i), data.P(i), data.U(i)); 
    fprintf(fp3,"%6d %4d %14.8f", i, data.IType(i), data.Hsml(i));                                
  }
       
  fclose(fp1); 
  fclose(fp2);
  fclose(fp3);
}

//---------------------------------------------------------------------- 
// Subroutine to determine the right hand side of a differential  
// equation in a single step for performing time integration  
 
// In this routine and its subroutines the SPH algorithms are performed. 
//   data.ITimeStep: Current timestep number                            [in] 
//   data.Dt       : Timestep                                           [in] 
//   data.NTotal   :  Number of particles                               [in] 
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
 
void single_step(SPH_Data &data)  
{
  int i, d;
  MkInt ns(data.MaxN+1);
  MkDouble indvxdt(data.Dim+1,data.MaxN+1);
  MkDouble exdvxdt(data.Dim+1,data.MaxN+1),ardvxdt(data.Dim+1,data.MaxN+1),avdudt(data.MaxN+1),ahdudt(data.MaxN+1);

  //  MkDebug("\n  single_step()\n");
  //  MkDebug("  S1 ");
  for ( i=1;i<=data.NTotal;i++) {
    avdudt(i) = 0.; 
    ahdudt(i) = 0.; 
    for ( d=1;d<=data.Dim;d++) {
      indvxdt(d,i) = 0.; 
      ardvxdt(d,i) = 0.; 
      exdvxdt(d,i) = 0.;
    }
  }
  //  MkDebug("  S2 ");  
//---  Positions of virtual (boundary) particles:  
 
  data.NVirt = 0; 
  if (virtual_part) {
    virt_part(data); 
  }
  //  MkDebug("  S3 ");      
//---  Interaction parameters, calculating neighboring particles 
//   and optimzing smoothing length 
 
  if (nnps==1) { direct_find(data) ; ns.CopyFrom(data.CountIac);  }
  else if (nnps==2) { link_list(data); ns.CopyFrom(data.CountIac);  }
//  else if (nnps==3) 
//    tree_search(data.ITimeStep, data.NTotal+nvirt,data.Hsml,x,niac,pair_i,pair_j,w,dwdx,ns);
//  MkDebug("  S4 ");
//---  Density approximation or change rate 

  if (summation_density) sum_density(data);
  else con_density(data);
 
//---  Dynamic viscosity: 
 
//  MkDebug("  S5 ");
  if (visc) viscosity(data); 
        
//---  Internal forces: 
  
  int_force(data) ; indvxdt.CopyFrom(data.DVXDt); data.DUDt.CopyFrom(data.DEDt);
               
  //  MkDebug("  S6 ");    

//---  Artificial viscosity: 
 
  if (visc_artificial) {art_visc(data); ardvxdt.CopyFrom(data.DVXDt);avdudt.CopyFrom(data.DEDt);}
  //  MkDebug("  S7 ");
       
//---  External forces: 
 
//  MkDebug("  S8 ");

  if (ex_force) {ext_force(data);exdvxdt.CopyFrom(data.DVXDt);}
 
  //  MkDebug("  S9 ");

//   Calculating the neighboring particles and undating DATA.HSML 
       
  if (sle!=0) h_upgrade(data); 
 
  if (heat_artificial) {art_heat(data);ahdudt.CopyFrom(data.DEDt);}
      
//   Calculating average velocity of each partile for avoiding penetration 
 
  if (average_velocity) av_vel(data);
 
//---  Convert velocity, force, and energy to f and dfdt   
 
  for (i=1;i<=data.NTotal;i++) {
    for (d=1;d<=data.Dim;d++) {
      data.DVXDt(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i);
    }
    data.DUDt(i) = data.DUDt(i) + avdudt(i) + ahdudt(i);
  }
  //   MkDebug("  S10 ");
  if ((data.ITimeStep%print_step)==0) {
    //      123456789abc 123456789abc 123456789abc
    printf("\n") ;
    printf("**** Information for particle **** %d\n",moni_particle);
    printf("internal a   artifical a  external a    total a \n");
    printf("%12.6f %12.6f %12.6f %12.6f\n",indvxdt(1,moni_particle),ardvxdt(1,moni_particle),exdvxdt(1,moni_particle),data.DVXDt(1,moni_particle));
  }

  //  MkDebug("------------------ Data Dumping %d --------------------\n",data.ITimeStep);
  //data.Dump();
  //MkDebug("------------------ ~Data Dumping --------------------\n");

  //  MkDebug("  ~single_step()\n");
}
 
//---------------------------------------------------------------------- 
// Subroutine to define the fluid particle viscosity 
  
//   data.NTotal  : Number of particles                                 [in] 
//   data.IType    : Type of particle                                   [in] 
//   data.X       : Coordinates of all particles                        [in] 
//   data.Rho     : Density                                             [in] 
//   data.Eta     : Dynamic viscosity                                  [out] 
 

void viscosity(SPH_Data &data) 
{       
  int i; 
  for (i=1;i<=data.NTotal+data.NVirt;i++) {
    if (abs(data.IType(i))==1) data.Eta(i)=0.;
    else if (abs(data.IType(i))==2) data.Eta(i)=1.0e-3;
  }
}

//---------------------------------------------------------------------- 
//    data.X-- coordinates of particles                       [input/output] 
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
//    data.NTotal-- total particle number                            [input]   
//    data.MaxTimeStep-- maximum timesteps                           [input] 
//    data.Dt-- timestep                                             [input] 
    
void time_integration(SPH_Data &data) 
{
  int i, j, k, d, current_ts, nstart;
  double  time, temp_rho, temp_u;

  MkDouble  x_min(data.Dim+1, data.MaxN+1), v_min(data.Dim+1, data.MaxN+1), u_min(data.MaxN+1);
  MkDouble rho_min(data.MaxN+1); 

  //  MkDouble  dvx(data.Dim+1, data.MaxN+1);
  i= j= k= d= current_ts= nstart=0;
  time= temp_rho= temp_u=0;
           
  //  MkDebug("time_integration()\n");
     
  for (i = 1;i<= data.NTotal;i++) {
    for (d = 1;d<= data.Dim;d++) {
      data.AveVel(d, i) = 0.;
    }
  }

  for (data.ITimeStep = nstart+1;data.ITimeStep<= nstart+data.MaxTimeStep;data.ITimeStep++) {
    current_ts=current_ts+1;
    if ((data.ITimeStep%print_step)==0) {
      printf("______________________________________________\n");
      printf("  current number of time step = %d     current time=%f\n",data.ITimeStep, float(time+data.Dt)) ;
      printf("______________________________________________\n");
    }
       
//   If not first time step, then update thermal energy, density and  
//   velocity half a time step   

//    MkDebug("  T1 ");
 
    if (data.ITimeStep != 1) {
 
      for (i = 1;i<= data.NTotal;i++) {
	u_min(i) = data.U(i);
	temp_u=0.;
	if (data.Dim==1) temp_u=-nsym*data.P(i)*data.VX(1,i)/data.X(1,i)/data.Rho(i);
	data.U(i) = data.U(i) + (data.Dt/2.)* (data.DUDt(i)+temp_u);
	if(data.U(i)<0)  data.U(i) = 0.;
             
	if (!summation_density) {
	  rho_min(i) = data.Rho(i);
	  temp_rho=0.; 
	  if (data.Dim==1) temp_rho=-nsym*data.Rho(i)*data.VX(1,i)/data.X(1,i);
	  data.Rho(i) = data.Rho(i) +(data.Dt/2.)*( data.DRhoDt(i)+ temp_rho);
	}

	for (d = 1;d<= data.Dim;d++) {
	  v_min(d, i) = data.VX(d, i) ;
	  data.VX(d, i) = data.VX(d, i) + (data.Dt/2.)*data.DVXDt(d, i) ;
	}
      }
    }
//---  Definition of variables out of the function vector:     

//    MkDebug("  T2 ");       

    single_step(data);
                 
    //    MkDebug("  T3 ");
  
    if (data.ITimeStep == 1) { 
        
      for (i=1;i<=data.NTotal;i++) {
	temp_u=0.;
	if (data.Dim==1) temp_u=-nsym*data.P(i)*data.VX(1,i)/data.X(1,i)/data.Rho(i);
	data.U(i) = data.U(i) + (data.Dt/2.)*(data.DUDt(i) + temp_u); 
	if(data.U(i)<0)  data.U(i) = 0.;
          
	if (!summation_density ){
	  temp_rho=0.; 
	  if (data.Dim==1) temp_rho=-nsym*data.Rho(i)*data.VX(1,i)/data.X(1,i) ;
	  data.Rho(i) = data.Rho(i) + (data.Dt/2.)* (data.DRhoDt(i)+temp_rho) ;
	}
          
	//    MkDebug("  T4 ");
	for (d = 1;d<= data.Dim;d++) {
	  data.VX(d, i) = data.VX(d, i) + (data.Dt/2.) * data.DVXDt(d, i) + data.AveVel(d, i); 
	  data.X(d, i) = data.X(d, i) + data.Dt * data.VX(d, i);
	}
      }
    }               
    else {
      //MkDebug("  T5 ");
                     
      for (i=1;i<=data.NTotal;i++) {
	temp_u=0.;
	if (data.Dim==1) temp_u=-nsym*data.P(i)*data.VX(1,i)/data.X(1,i)/data.Rho(i);
	data.U(i) = u_min(i) + data.Dt*(data.DUDt(i)+temp_u);
	if(data.U(i)<0)  data.U(i) = 0.;
             
	if (!summation_density ) {
	  temp_rho=0.;
	  if (data.Dim==1) temp_rho=-nsym*data.Rho(i)*data.VX(1,i)/data.X(1,i);
	  data.Rho(i) = rho_min(i) + data.Dt*(data.DRhoDt(i)+temp_rho);
	}
                 
	for (d = 1;d<= data.Dim;d++) {
	  data.VX(d, i) = v_min(d, i) + data.Dt * data.DVXDt(d, i) + data.AveVel(d, i);
	  data.X(d, i) = data.X(d, i) + data.Dt * data.VX(d, i);
	}
      }
         
    }      
    time = time + data.Dt;
 
    if ((data.ITimeStep%save_step)==0) {
      output(data);
    }

    //    MkDebug("  T10 ");
  
    if ((data.ITimeStep%print_step)==0) {
      //      123456789abc 123456789abc 123456789abc 123456789abc 
      printf("\n");
      printf("x            velocity     dvx\n");
      printf("%12.6f %12.6f %12.6f \nf",data.X(1,moni_particle), data.VX(1,moni_particle),data.DVXDt(1,moni_particle));
    }
  }
  nstart=current_ts;
  //  MkDebug("~time_integration()\n");

}
 
//---------------------------------------------------------------------- 
// Subroutine to determine the information of virtual particles 
// Here only the Monaghan type virtual particles for the 2D shear 
// cavity driven problem are generated. 
//   data.ITimeStep : Current time step                                 [in] 
//   data.NTotal : Number of particles                                  [in] 
//   nvirt  : Number of virtual particles                         [out] 
//   data.Hsml   : Smoothing Length                                 [in|out] 
//   data.Mass   : Particle masses                                  [in|out] 
//   data.X      : Coordinates of all particles                     [in|out] 
//   data.VX     : Velocities of all particles                      [in|out] 
//   data.Rho    : Density                                          [in|out] 
//   data.U      : internal energy                                  [in|out] 
//   data.IType   : type of particles                               [in|out] 
 
void virt_part(SPH_Data &data)  
{
  int i, j, d, im, mp; 
  double xl, dx, v_inf; 
  FILE *fp1,*fp2, *fp3;

   i= j= d= im= mp=0; 
   xl= dx= v_inf=0; 
 
  if (vp_input){
                         
    fp1 = fopen("xv_vp.dat","r"); 
    fp2 = fopen("state_vp.dat","r"); 
    fp3 = fopen("other_vp.dat","r");
    fscanf(fp1,"%d",&data.NVirt); 
    for (j = 1;j<= data.NVirt;j++) {
      i = data.NTotal + j;
      fscanf(fp1,"%d",&im);
      for (d=1;d<=data.Dim;d++) fscanf(fp1,"%f",&data.X(d, i));
      for (d=1;d<=data.Dim;d++) fscanf(fp1,"%f",&data.VX(d, i));
      fscanf(fp2,"%d %f %f %f %f",im, data.Mass(i), data.Rho(i), data.P(i), data.U(i));
      fscanf(fp3,"%d %d %f",im, data.IType(i), data.Hsml(i));
    }
    fclose(fp1); 
    fclose(fp2); 
    fclose(fp3);

  }
  else  {
    data.NVirt = 0; 
    mp = 40;
    xl = 1.0e-3; 
    dx = xl / mp; 
    v_inf = 1.e-3; 
 
//   Monaghan type virtual particle on the Upper side 
 
    for (i = 1;i<= 2*mp+1;i++) {
      data.NVirt = data.NVirt + 1;
      data.X(1, data.NTotal + data.NVirt) = (i-1)*dx/2;
      data.X(2, data.NTotal + data.NVirt) = xl;   
      data.VX(1, data.NTotal + data.NVirt) = v_inf; 
      data.VX(2, data.NTotal + data.NVirt) = 0.; 
    }
 
//   Monaghan type virtual particle on the Lower side 
 
    for (i = 1;i<= 2*mp+1;i++) {
      data.NVirt = data.NVirt + 1; 
      data.X(1, data.NTotal + data.NVirt) = (i-1)*dx/2;  
      data.X(2, data.NTotal + data.NVirt) = 0.;   
      data.VX(1, data.NTotal + data.NVirt) = 0.; 
      data.VX(2, data.NTotal + data.NVirt) = 0.; 
    }
 
//   Monaghan type virtual particle on the Left side 
 
    for (i = 1;i<= 2*mp-1;i++) {
      data.NVirt = data.NVirt + 1;
      data.X(1, data.NTotal + data.NVirt) = 0.;
      data.X(2, data.NTotal + data.NVirt) = i*dx/2;
      data.VX(1, data.NTotal + data.NVirt) = 0.; 
      data.VX(2, data.NTotal + data.NVirt) = 0.;
    }
 
//   Monaghan type virtual particle on the Right side 
 
    for (i = 1;i<= 2*mp-1;i++) {
      data.NVirt = data.NVirt + 1;
      data.X(1, data.NTotal + data.NVirt) = xl;
      data.X(2, data.NTotal + data.NVirt) = i*dx/2;
      data.VX(1, data.NTotal + data.NVirt) = 0.; 
      data.VX(2, data.NTotal + data.NVirt) = 0.; 
    }
 
    for (i = 1;i<= data.NVirt;i++) {
      data.Rho (data.NTotal + i) = 1000.;
      data.Mass(data.NTotal + i) = data.Rho (data.NTotal + i) * dx * dx; 
      data.P(data.NTotal + i) = 0.; 
      data.U(data.NTotal + i) = 357.1; 
      data.IType(data.NTotal + i) = -2; 
      data.Hsml(data.NTotal + i) = dx; 
    }
  } 
  if ((data.ITimeStep%save_step)==0) {
    fp1 = fopen("xv_vp.dat","w"); 
    fp2 = fopen("state_vp.dat","w"); 
    fp3 = fopen("other_vp.dat","w");
    fprintf(fp1,"\n", data.NVirt) ;
    for (i = data.NTotal + 1;i<= data.NTotal + data.NVirt;i++) {
      fprintf(fp1,"%6d ",i);
      for (d=1;d<=data.Dim;d++) fprintf(fp1,"%14.8f ", data.X(d, i));
      for (d=1;d<=data.Dim;d++) fprintf(fp1,"%14.8f ", data.VX(d, i));
      fprintf(fp1,"\n");
      fprintf(fp2,"%6d %14.8f %14.8f %14.8f %14.8f\n", i, data.Mass(i), data.Rho(i), data.P(i), data.U(i));
      fprintf(fp3,"%6d %4d %14.8f %14.8f\n", i, data.IType(i), data.Hsml(i));                         
    }
    fclose(fp1); 
    fclose(fp2); 
    fclose(fp3);
  }
 
  if ((data.ITimeStep%print_step)==0) {
    if (int_stat) {
      printf(" >> Statistics: Virtual boundary particles:\n"); 
      printf("          Number of virtual particles:%d\n",data.NVirt);
    }
  }
}
 
//---------------------------------------------------------------------- 
//   This is a three dimensional SPH code. the followings are the  
//   basic parameters needed in this code or calculated by this code 
 
//   data.Mass-- mass of particles                                      [in] 
//   data.NTotal-- total particle number ues                            [in] 
//   data.Dt--- Time step used in the time integration                  [in] 
//   data.IType-- types of particles                                    [in] 
//   data.X-- coordinates of particles                              [in/out] 
//   data.VX-- velocities of particles                              [in/out] 
//   data.Rho-- densities of particles                              [in/out] 
//   data.P-- pressure  of particles                                [in/out] 
//   data.U-- internal energy of particles                          [in/out] 
//   data.Hsml-- smoothing lengths of particles                     [in/out] 
//   data.C-- sound velocity of particles                              [out] 

void initLiuSPH()
{

  SPH_Data data;
  int yesorno=1;       
  
  if (shocktube)   data.Dt = 0.005;
  if (shearcavity) data.Dt = 5.e-5;
 
  MkDebug("program start!\n");

  data.Initialize(12000,3);

  MkDebug("data initialized!\n");
  input(data);

  while(yesorno) {
    printf("  ***************************************************\n"); 
    printf("          Please input the maximal time steps \n");
    printf("  ***************************************************\n");
    scanf("%d", &data.MaxTimeStep);
    MkDebug("    maxtimestep = %d \n",data.MaxTimeStep);
    time_integration(data);
    output(data);
    printf("  ***************************************************\n"); 
    printf("  Are you going to run more time steps ? (0=No, 1=yes)\n"); 
    printf("  ***************************************************\n");
    scanf ("%d ", &yesorno);
  }
  MkDebug("normal termination!\n");
  exit(-1);

}

int main(int argc, char **argv)
{ 


  //  particles();
  InitLiuSPH();
   
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

  return -1;


}
