
#include "MkSparseMatrix.h"
//---------------------------------------------------------------------------
MkSMij::MkSMij()
{
   Value = 0;
   FI = 0;
   FJ = 0;
   NextSM = NULL;
   PrevSM = NULL;
}

MkSMij::MkSMij(int i,int j)
{
   Value = 0;
   FI = i;
   FJ = j;
   NextSM = NULL;
   PrevSM = NULL;
}

MkSMij::MkSMij(int i,int j,double v)
{
   Value = v;
   FI = i;
   FJ = j;
   NextSM = NULL;
   PrevSM = NULL;
}

MkSMij::MkSMij(MkSMij *next,MkSMij *prev)
{
   Value = 0;
   FI = 0;
   FJ = 0;
   NextSM = next;
   PrevSM = prev;
}

MkSMij::MkSMij(int i, int j,MkSMij *next, MkSMij *prev)
{
   Value = 0;
   FI = i;
   FJ = j;
   NextSM = next;
   PrevSM = prev;
}

void change(MkSMij *sma, MkSMij *smb)
{
  MkSMij *tmp;
  if(!sma || !smb) return;
  tmp = sma;
  sma->NextSM = smb->NextSM;
  sma->PrevSM = smb->PrevSM;
  sma->NextSM->PrevSM = smb->NextSM->PrevSM;
  sma->PrevSM->NextSM = smb->PrevSM->NextSM;
  smb->NextSM = tmp->NextSM;
  smb->PrevSM = tmp->PrevSM;
  smb->NextSM->PrevSM = tmp->NextSM->PrevSM;
  smb->PrevSM->NextSM = tmp->PrevSM->NextSM;
}

bool MkSMij::operator==(MkSMij & sm)
{
   return FI==sm.FI && FJ==sm.FJ;
}

MkSMij & MkSMij::operator=(MkSMij &sm)
{
   Value = sm.Value;
   FI = sm.FI;
   FJ = sm.FJ;
   NextSM = sm.NextSM;
   PrevSM = sm.PrevSM;
   return *this;
}

// it is different between
// sm1 == sm2   (pointer comparison)
// *sm1 == *sm2 (call operator)
//---------------------------------------------------------------------------
MkSparseMatrix::MkSparseMatrix()
{
   M = 0;
   N = 0;
   Top=NULL;
   Cur=NULL;
   Last=NULL;
   Up = NULL;
   Low = NULL;
   MaxIteration = 1000;
   Tol = FTOL;
   FMatType = mtNormal;
}

MkSparseMatrix::~MkSparseMatrix()
{
   int cnt=0;
   MkDebug("Deallocating...\n");
   while(Top) {
     cnt++;
     MkDebug(" Top(%d) : (%d,%d),%f is being deleted\n",cnt,Top->FI,Top->FJ,Top->Value);
     Del(Top);
     for(long i=0;i<pow(3,10);i++);
   }
   MkDebug("Deallocation finished...\n");
   M = 0;
   N = 0;
   Top=NULL;
   Cur=NULL;
   Last=NULL;
   Up = NULL;
   Low = NULL;
   FMatType = mtNormal;
}

MkSparseMatrix::MkSparseMatrix(int m,int n)
{
   M = m;
   N = n;
   Top = NULL;
   Cur = NULL;
   Last = NULL;
   Up = NULL;
   Low = NULL;
   MaxIteration = 1000;
   Tol = FTOL;
   FMatType = mtNormal;
   if(M==N) 
     for (int i=0;i<M;i++) {
   	 Add(new MkSMij(i,i));
     }
   FIndex.Initialize(N);
}

MkSparseMatrix::MkSparseMatrix(MkMatrix &mat)
{
  SetMatrix(mat);
}

MkSparseMatrix::MkSparseMatrix(MkSparseMatrix &mat)
{
  SetMatrix(mat);
}

void MkSparseMatrix::SetMatrix(MkMatrix &mat)
{
  int i,j;
  M = mat.GetFI();
  N = mat.GetFJ();
  MaxIteration = 1000;
  Tol = FTOL;
  FMatType = mat.GetMatType();

  FIndex.CopyFrom(mat.GetIndex());
  for(i=0;i<(mat.GetFI()<mat.GetFJ()?mat.GetFI():mat.GetFJ());i++) {
    if(mat(i,i)>TINY) Add(i,i,mat(i,i));
  }
  for(i=1;i<mat.GetFI();i++) {
    for(j=0;i>j;j++) {
      if(mat(i,j)>FTOL) Add(i,j,mat(i,j));
    }
  }
  for(j=1;j<mat.GetFJ();j++) {
    for(i=0;i<j;i++) {
      if(mat(i,j)>FTOL) Add(i,j,mat(i,j));
    }
  }
}

void MkSparseMatrix::SetMatrix(MkSparseMatrix &mat)
{
  MkSMij *sm;
  M = mat.M;
  N = mat.N;
  MaxIteration = 1000;
  Tol = FTOL;
  FMatType = mat.FMatType;

  FIndex.CopyFrom(mat.GetIndex());
  sm = mat.Top;
  while(sm) {
    Add(sm->FI,sm->FJ,sm->Value);
    sm=sm->NextSM;
  }
}

void MkSparseMatrix::Print(char *fname)
{
  FILE *fp;
  fp = fopen(fname,"w");
  if(!fp) return;
  MkSMij *sm;
  sm = Top;

  if(Top) {
    fprintf(fp,"Top index:(%d,%d), value:%f\n",Top->FI,Top->FJ,Top->Value);
    MkDebug("Top index:(%d,%d), value:%f\n",Top->FI,Top->FJ,Top->Value);
  }
  if(Last) {
    fprintf(fp,"Last index:(%d,%d), value:%f\n",Last->FI,Last->FJ,Last->Value);
    MkDebug("Last index:(%d,%d), value:%f\n",Last->FI,Last->FJ,Last->Value);
  }
  if(Low) {
    fprintf(fp,"Low index:(%d,%d), value:%f\n",Low->FI,Low->FJ,Low->Value);
    MkDebug("Low index:(%d,%d), value:%f\n",Low->FI,Low->FJ,Low->Value);
  }
  if(Up) {
    fprintf(fp,"Up index:(%d,%d), value:%f\n",Up->FI,Up->FJ,Up->Value);
    MkDebug("Up index:(%d,%d), value:%f\n",Up->FI,Up->FJ,Up->Value);
  }
  while(sm) {
    fprintf(fp,"index:(%d,%d), value:%f\n",sm->FI,sm->FJ,sm->Value);
    MkDebug("index:(%d,%d), value:%f\n",sm->FI,sm->FJ,sm->Value);
    sm=sm->Next();
  }
  fclose(fp);
}

void MkSparseMatrix::Print()
{
  MkSMij *sm;
  sm = Top;

  if(Top) {
    MkDebug("Top index:(%d,%d), value:%f\n",Top->FI,Top->FJ,Top->Value);
  }
  if(Last) {
    MkDebug("Last index:(%d,%d), value:%f\n",Last->FI,Last->FJ,Last->Value);
  }
  if(Low) {
    MkDebug("Low index:(%d,%d), value:%f\n",Low->FI,Low->FJ,Low->Value);
  }
  if(Up) {
    MkDebug("Up index:(%d,%d), value:%f\n",Up->FI,Up->FJ,Up->Value);
  }
  while(sm) {
    MkDebug("index:(%d,%d), value:%f\n",sm->FI,sm->FJ,sm->Value);
    sm=sm->Next();
  }
}
              
void MkSparseMatrix::Identity(bool clear)
{
  MkSMij *sm;
  if(clear) {
    Clear();
    for(int i=0;i<(M<N?M:N);i++)
      Add(i,i,1.0);
  }
  else {
    Cur=Top;
    while(Cur&&Cur->FI==Cur->FJ) {
      Cur->Value = 1.0;
      Cur=Cur->Next();
    }
    Cur = Low?Low:Up;
    while(Cur) {
      Cur->Value = 0.0;
      Cur=Cur->Next();
    }
  }
}

void MkSparseMatrix::Clear()
{
  MkSMij *sm;
  Cur = Top;
  while(Cur) {
    sm = Cur;
    Cur=Cur->NextSM;
    delete sm;
    sm = NULL;
  }
  Top = Last = Low = Up = Cur = NULL;
}

bool MkSparseMatrix::Transpose()
{
  MkSMij temp,*sm;
  int tmp;
  temp = *Low;
  Low->PrevSM = Up->PrevSM;
  Up->PrevSM = temp.PrevSM;
  sm = Low;
  while(sm) {
    tmp = sm->FI;
    sm->FI = sm->FJ;
    sm->FJ = tmp;
    sm=sm->NextSM;
  }

  tmp = M;
  M = N;
  N = tmp;
  
  FMatType = FMatType == mtNormal ? mtTransposed : mtNormal;
  return true;
}

bool MkSparseMatrix::Invert()
{
  return false;
}

bool MkSparseMatrix::LUDecompose()
{
  int i,imax,j,k,sign=1;
  double big,dum,sum,temp;
  MkDouble vv,diag;
  MkSMij *sm,*low,*up;
  bool is_exist=false;

  if(FMatType == mtLUDecomposed) return true;
  if (M != N) return false;
  int n = M;

  vv.Initialize(n);
  diag.Initialize(n);
  FD=1;

  sm = Top;
  while(sm) {
    vv(sm->FI)=fabs(sm->Value)>fabs(vv(sm->FI))?vv(sm->FI)=sm->Value:vv(sm->FI);
    sm=sm->NextSM;
  }
  
  for (i=0;i<n;i++) {
    if(fabs(vv(i)) < EPS)  MkDebug("Singular matrix in routine MkMatrix::LUDecompose");
    vv(i)=1.0/vv(i);
  }
  
  for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
          sum=(*this)(i,j);
          if(fabs(sum)<TINY) false;
          low = FindFirstInRow(i);
          up = FindFirstInCol(j);
          if(!low||!up) continue; // if there is no element for low or up, no need to calculate
          while(low->FI==i && up->FJ==j) {
            if(low->FJ==up->FI) {
               sum -= low->Value*up->Value;
            }
            sign>0?low=low->NextSM:up=up->NextSM;
            sign*=-1;
          }
          if(fabs(sum)>TINY&& is_exist) (*this)(i,j) = sum;
          else if(fabs(sum)>TINY && !is_exist) Add(i,j,sum);
      }
      big=0.0;
      for (i=j;i<n;i++) {
          sum=(*this)(i,j);
          if(fabs(sum)<TINY) false;
          low = FindFirstInRow(i);
          up = FindFirstInCol(j);
          if(!low||!up) continue; // if there is no element for low or up, no need to calculate
          while(low->FI==i && up->FJ==j) {
            if(low->FJ==up->FI) {
               sum -= low->Value*up->Value;
            }
            sign>0?low=low->NextSM:up=up->NextSM;
            sign*=-1;
          }
          if(fabs(sum)>TINY&& is_exist) (*this)(i,j) = sum;
          else if(fabs(sum)>TINY && !is_exist) Add(i,j,sum);

          if ( (dum=vv(i)*fabs(sum)) >= big) {
             big=dum;
             imax=i;
          }
      }
      if (j != imax) {
         MkSMij *smi,*smj;
         smi = Top;
         smj = Top;
         while(smi) {
           if(smi->FI==imax)
             IsExist(j,smi->FJ)?swap(smi->Value,Cur->Value):(void)Move(smi,j,smi->FJ);
           smi->NextSM;
         }

         while(smj) {
           if(smj->FI==j){
             if(!IsExist(imax,smj->FJ)) Move(smj,imax,smj->FJ);
           }
         }

         FD = -(FD);
         vv(imax)=vv(j);
      }
      FIndex(j)=imax;
  }

  sm=Top;
  while(sm->FI==sm->FJ) {
    diag(sm->FI) = sm->Value;
    sm=sm->NextSM;
  }
  
  sm=Low;
  while(sm->FI>sm->FJ) {
    sm->Value/=fabs(diag(sm->FJ))>TINY?diag(sm->FJ):TINY;
    sm=sm->NextSM;
  }

  FMatType = mtLUDecomposed;
  return true;
}

bool MkSparseMatrix::LUBackSubstitute(MkVector &b)
{
    int i,ii=-1,ip,j;
    double sum;
    MkDouble diag;
    MkSMij *sm;
    int n;

    if (M != N) return false;

    n = M;
    if (FMatType != mtLUDecomposed) return false;

    diag.Initialize(n);

    sm=Top;
    while(sm->FI==sm->FJ) {
      diag(sm->FI) = sm->Value;
      sm=sm->NextSM;
    }
    
    arrangeByNormal();
    FMatType = mtBackSubstitute;

    for(sm=Low,i=0;i<n;i++) {
      ip=FIndex(i);
      sum=b(ip);
      b(ip)=b(i);
      for(sm = FindFirstInRow(i);sm->FI==i&&0<=sm->FJ&&sm->FJ<i;sm=sm->NextSM)
        sum -= sm->Value*b(sm->FJ);
      b(i)=sum;
    }

    for (sm=Low,i=n-1;i>=0;i--) {
      sum=b(i);
      for(sm = FindFirstInRow(i);sm->FI==i&&i+1<=sm->FJ&&sm->FJ<n;sm=sm->NextSM)
        sum -= sm->Value*b(sm->FJ);
      b(i)=sum/diag(i);
    }
    return true;
}

bool MkSparseMatrix::GaussSeidel(MkVector &X0,MkVector &B)
{
    int i,m=0,I,J,iter=0;
    MkDouble X(N);
    MkDouble U(N);
    MkSMij *sm;
    bool isConverged=false;

    if(FMatType == mtLUDecomposed) return false;
    if (M != N) return false;

    arrangeByGauss();

    for(i=0;i<N;i++){
      X(i) = (*this)(i,i);
      if(fabs(X(i))>EPS) B(i) = B(i)/X(i);
      else return false;
    }

    sm = Low;
    while(sm) {
      if(fabs(X(sm->FI))>EPS) sm->Value = sm->Value/X(sm->FI);
      else return false;
      sm->NextSM;
    }

    X.CopyFrom(X0.GetDouble());
    while(!isConverged&&iter<MaxIteration) {
      sm = Low;
      U.CopyFrom(X);
      while(sm) {
        for(I=sm->FI,X(I)=0;sm?I==sm->FI:false;sm=sm->NextSM) {
          X(I) -= sm->Value*X(sm->FJ);
        }
      }
      for(i=0;i<N;i++)
        X(i) += B(i);

      isConverged = true;
      for(i=0;i<N;i++) {
        if(fabs(X(i)-U(i))>Tol) {
          isConverged = false;
          break;
        }
      }
    }
    return true;
}

bool MkSparseMatrix::isSingular()
{
    int i,dum=1;

    if (M != N) return true;

    MkInt A(M);
    for (i=0;i<M;i++) A(i) = 0;

    MkSMij *cur;
    cur = Top;
    while(cur) {
      if(fabs(cur->Value)>EPS) A(cur->FI) = 1;
      cur=cur->Next();
    }

    for(i=0;i<M;i++)dum*=A(i);
    return dum==0?true:false;
}

bool MkSparseMatrix::Solve(MkVector &B)
{
    if (isSingular()) return false;
    MkVector X;
    MkSparseMatrix A(*this);
    X = B;
    A.LUDecompose();
    A.LUBackSubstitute(X);
    MkSparseMatrix A1(*this);
    A1.GaussSeidel(X,B);
    B = X;
    return true;
}

bool MkSparseMatrix::Solve(MkVector &B,SolveType solve_type)
{
    if (isSingular()) return false;

    if (solve_type == stLUD) {
       MkVector X;
       MkSparseMatrix A(*this);
       X = B;
       A.LUDecompose();
       A.LUBackSubstitute(X);
       B = X;
       return true;
    }
    else if (solve_type == stGauss) {
       MkVector X;
       MkSparseMatrix A(*this);
       X = B;
       A.GaussSeidel(X,B);
       B = X;
       return true;
    }
    else if (solve_type == stHybrid)
       return Solve(B);
    else return false;
}

bool MkSparseMatrix::Add(int ai,int aj)
{
  if(ai>M||aj>N) return false;
  return Add(new MkSMij(ai,aj,1.0));
}

bool MkSparseMatrix::Add(int ai,int aj,double v)
{
   if(ai>M||aj>N) return false;
   return Add(new MkSMij(ai,aj,v));
}

// Top -> Low -> Up -> Last

bool MkSparseMatrix::Add(MkSMij *sm)
{
   if(FMatType==mtNormal) return addByNormal(sm);
   else if(FMatType==mtGauss) return addByGauss(sm);
   else if(FMatType==mtCol) return addByCol(sm);
   else if(FMatType==mtRow) return addByRow(sm);
   return false;
}

bool MkSparseMatrix::Del(int ai,int aj)
{
  if(!IsExist(ai,aj)) return false;
  MkSMij *sm;
  sm = Cur;
  if(sm->FI==ai && sm->FJ==aj) {
    return Del(sm);
  }
  return false;
}

bool MkSparseMatrix::Del(MkSMij *sm)
{
   if(FMatType==mtNormal) return delByNormal(sm);
   else if(FMatType==mtGauss) return delByGauss(sm);
   else if(FMatType==mtCol) return delByCol(sm);
   else if(FMatType==mtRow) return delByRow(sm);
   return false;
}

bool MkSparseMatrix::Move(int bi, int bj, int ai, int aj)
{
  MkSMij *sm,*next,*prev;
  if(!IsExist(bi,bj)) return false;
  sm = Cur;
  return Move(sm,ai,aj);
}

bool MkSparseMatrix::Move(MkSMij *sm, int ai, int aj)
{
   if(FMatType==mtNormal) return moveByNormal(sm,ai,aj);
   else if(FMatType==mtGauss) return moveByGauss(sm,ai,aj);
   else if(FMatType==mtCol) return moveByCol(sm,ai,aj);
   else if(FMatType==mtRow) return moveByRow(sm,ai,aj);
   return false;
}

bool MkSparseMatrix::IsExist(MkSMij *sm)
{
   int ai,aj;
   ai = sm->FI;
   aj = sm->FJ;
   return IsExist(ai,aj);
}

bool MkSparseMatrix::IsExist(int ai,int aj)
{
  MkSMij *next,*prev,sm(ai,aj);
  if(!Top) return false;
  if(!Cur) Cur = Top;
  while(Cur) {
    next = Cur->Next();
    prev = Cur->Prev();
    if(isEQ(&sm,Cur)) return true;
    if(isGT(&sm,Cur)) {
      if(!next) return false;
      else if(isGT(&sm,next)) {Cur=next;continue;}
      else if(isLT(&sm,next)) return false;
    }
    else if(isLT(&sm,Cur)) {
      if(!prev) return false;
      else if(isLT(&sm,prev)) {Cur=prev;continue;}
      else if(isGT(&sm,prev)) {Cur=prev;return false;}
    }
  }
  if(!Cur) return false;//if it happens, it is really strange!!!
  else return true;
}

bool MkSparseMatrix::FindLoc(int ai,int aj)
{
  return IsExist(ai,aj);
}

MkSMij *MkSparseMatrix::FindFirstInRow(int ai)
{
   if(FMatType==mtNormal) return findFirstInRowByNormal(ai);
   else if(FMatType==mtGauss) return findFirstInRowByGauss(ai);
   else if(FMatType==mtCol) return findFirstInRowByCol(ai);
   else if(FMatType==mtRow) return findFirstInRowByRow(ai);
   return SMNull;
}

MkSMij *MkSparseMatrix::FindFirstInCol(int aj)
{
   if(FMatType==mtNormal) return findFirstInColByNormal(aj);
   else if(FMatType==mtGauss) return findFirstInColByGauss(aj);
   else if(FMatType==mtCol) return findFirstInColByCol(aj);
   else if(FMatType==mtCol) return findFirstInColByCol(aj);
   return SMNull;
}

MkSMij *MkSparseMatrix::FindLastInRow(int ai)
{
   if(FMatType==mtNormal) return findLastInRowByNormal(ai);
   else if(FMatType==mtGauss) return findLastInRowByGauss(ai);
   else if(FMatType==mtCol) return findLastInRowByCol(ai);
   else if(FMatType==mtRow) return findLastInRowByRow(ai);
   return SMNull;
}

MkSMij *MkSparseMatrix::FindLastInCol(int aj)
{
   if(FMatType==mtNormal) return findLastInColByNormal(aj);
   else if(FMatType==mtGauss) return findLastInColByGauss(aj);
   else if(FMatType==mtCol) return findLastInColByCol(aj);
   else if(FMatType==mtCol) return findLastInColByCol(aj);
   return SMNull;
}

int MkSparseMatrix::CountFillIns()
{
  int count;
  double Total;
  Total = M*N;
  count=0;
  MkSMij *tmp;
  tmp = Top;
  while(tmp) {
    count++;
	tmp = tmp->Next();  
  }
  FillIns = double(count)/Total;
  return count;
}

bool MkSparseMatrix::Verify()
{
  MkSMij *cur;
  if(!Top) {
    if(Last || Cur || Up || Low) return false;
    else return true;
  }
  else if(!Last) {
    if(Top || Cur || Up || Low) return false;
    else return true;
  }
  else if(Top==Low) {
    cur = Top;
    while(cur&&cur!=Up) {
      if(cur->FI <= cur->FJ) break;
      cur=cur->NextSM;
    }
    if(!cur) return true;
    if(cur->FI <= cur->FJ) return false;

    if(Up) cur = Up;
    while(cur) {
      if(cur->FI >= cur->FJ) break;
      cur=cur->NextSM;
    }
    if(!cur) return true;
    if(cur->FI >= cur->FJ) return false;
  }
  else if(Top==Up) {
    cur = Top;
    while(cur) {
      if(cur->FI >= cur->FJ) break;
      cur=cur->NextSM;
    }
    if(!cur) return true;
    if(cur->FI >= cur->FJ) return false;
  }
  else {
    cur = Top;
    while(cur&&cur!=Low&&cur!=Up) {
      if(cur->FI != cur->FJ) break;
      cur=cur->NextSM;
    }
    if(!cur) return true;
    if(cur->FI != cur->FJ) return false;

    if(Low) cur = Low;
    while(cur&&cur!=Up) {
      if(cur->FI <= cur->FJ) break;
      cur=cur->NextSM;
    }
    if(!cur) return true;
    if(cur->FI <= cur->FJ) return false;

    if(Up) cur = Up;
    while(cur) {
      if(cur->FI >= cur->FJ) break;
      cur=cur->NextSM;
    }
    if(!cur) return true;
    if(cur->FI >= cur->FJ) return false;
  }
  return true;
}

double & MkSparseMatrix::operator()(int ai, int aj)
{
   static double Zero;Zero=0;
   if(!IsExist(ai,aj)) return Zero;
   if (Cur->FI == ai && Cur->FJ == aj) return Cur->Value;
   else return Zero;
}

//private member --------------------------------------------------------------
int MkSparseMatrix::orderNormal(MkSMij *sm)
{
  int n,I,J;
  n = M>N?N:M;
  I = sm->FI;
  J = sm->FJ;
  
  if(I>=n || J>=n) return -1;
  if(I==J) return I;
  else if(I>J) return I*(I-1)/2+J+1+n-1; 
  else if(I<J) return J*(J-1)/2+I+1+(n-1)*(n-2)/2+2*(n-1);
  return -1;
}

int MkSparseMatrix::orderGauss(MkSMij *sm)
{
  int n,I,J;
  n = M>N?N:M;
  I = sm->FI;
  J = sm->FJ;
  
  if(I>=n || J>=n) return -1;
  if(I==J) return I;
  else if(I!=J) return (I+1)*(n-1)+((I>J)?J+1:J); 
  return -1;
}

int MkSparseMatrix::orderCol(MkSMij *sm)
{
  int n,I,J;
  n = M>N?N:M;
  I = sm->FI;
  J = sm->FJ;
  
  if(I>=n || J>=n) return -1;
  return J*n+I; 
}

int MkSparseMatrix::orderRow(MkSMij *sm)
{
  int n,I,J;
  n = M>N?N:M;
  I = sm->FI;
  J = sm->FJ;
  
  if(I>=n || J>=n) return -1;
  return I*n+J; 
}

bool MkSparseMatrix::isGT(MkSMij *sm1,MkSMij *sm2)
{
  if(FMatType==mtNormal) return orderNormal(sm1)>orderNormal(sm2);
  else if(FMatType==mtGauss) return orderGauss(sm1)>orderGauss(sm2);
  else if(FMatType==mtCol) return orderCol(sm1)>orderCol(sm2);
  else if(FMatType==mtRow) return orderRow(sm1)>orderRow(sm2);
  else return false;
}

bool MkSparseMatrix::isLT(MkSMij *sm1,MkSMij *sm2)
{
  if(FMatType==mtNormal) return orderNormal(sm1)<orderNormal(sm2);
  else if(FMatType==mtGauss) return orderGauss(sm1)<orderGauss(sm2);
  else if(FMatType==mtCol) return orderCol(sm1)<orderCol(sm2);
  else if(FMatType==mtRow) return orderRow(sm1)<orderRow(sm2);
  else return false;
}

bool MkSparseMatrix::isEQ(MkSMij *sm1,MkSMij *sm2)
{
  if(FMatType==mtNormal) return orderNormal(sm1)==orderNormal(sm2);
  else if(FMatType==mtGauss) return orderGauss(sm1)==orderGauss(sm2);
  else if(FMatType==mtCol) return orderCol(sm1)==orderCol(sm2);
  else if(FMatType==mtRow) return orderRow(sm1)==orderRow(sm2);
  else return false;
}

bool MkSparseMatrix::isGE(MkSMij *sm1,MkSMij *sm2)
{
  if(FMatType==mtNormal) return orderNormal(sm1)>=orderNormal(sm2);
  else if(FMatType==mtGauss) return orderGauss(sm1)>=orderGauss(sm2);
  else if(FMatType==mtCol) return orderCol(sm1)>=orderCol(sm2);
  else if(FMatType==mtRow) return orderRow(sm1)>=orderRow(sm2);
  else return false;
}

bool MkSparseMatrix::isLE(MkSMij *sm1,MkSMij *sm2)
{
  if(FMatType==mtNormal) return orderNormal(sm1)<=orderNormal(sm2);
  else if(FMatType==mtGauss) return orderGauss(sm1)<=orderGauss(sm2);
  else if(FMatType==mtCol) return orderCol(sm1)<=orderCol(sm2);
  else if(FMatType==mtRow) return orderRow(sm1)<=orderRow(sm2);
  else return false;
}

bool MkSparseMatrix::addByNormal(MkSMij *sm)
{
  int I,J;
  MkSMij *next,*prev;
  
  I = sm->FI;
  J = sm->FJ;

  if(I<0||J<0) {
    delete sm;
    return false;
  }

  MkDebug(">> Finding location for (%d,%d)\n",sm->FI,sm->FJ);

  if (IsExist(sm)) {  // Find appropriate location to add
    delete sm;
    return false;
  }

  if(Top==NULL) {
    Top = sm;
    Top->NextSM = NULL;
    Top->PrevSM = NULL;
    Cur = sm;
    Last = sm;
    if(Top->FI > Top->FJ) Low = Top;
    else if(Top->FI < Top->FJ) Up = Top;
    return true;
  }

  if(!Cur) {
    delete sm;
    return false;
  }
    
  MkDebug("Cur (%d,%d),%f\n",Cur->FI,Cur->FJ,Cur->Value);
  MkDebug("sm is (%d,%d)\n",sm->FI,sm->FJ);
  
  next = Cur->Next();
  prev = Cur->Prev();

  if(I>J) {
    if(!Low) Low = sm;
    else if(isLT(sm,Low)) Low = sm;
  } 

  if(I<J) {
    if(!Up) Up = sm;
    else if(isLT(sm,Up)) Up = sm;
  }

  if(isGT(sm,Cur)) {
    if(!next){
      Cur->NextSM = sm;
      sm->PrevSM = Cur;
      Last = sm;
      return true;
    }
    else if(isGT(sm,next)) {
      return false; // strange!!!
    }
    else if(isLT(sm,next)) {
      Cur->NextSM = sm;
      sm->NextSM = next;
      next->PrevSM = sm;
      sm->PrevSM = Cur;
      return true;
    }
    return false;
  }
  else if(isLT(sm,Cur)) {
    if(!prev) {
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      Top = sm;
      return true;
    }
    else if(isLT(sm,prev)){
      return false;// strange!!!
    }
    else if(isGT(sm,prev)){
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      prev->NextSM = sm;
      sm->PrevSM = prev;
      return true;
    }
    return false;
  }
  else return false;
}

bool MkSparseMatrix::addByGauss(MkSMij *sm)
{
  int I,J;
  MkSMij *next,*prev;
  
  I = sm->FI;
  J = sm->FJ;

  if(I<0||J<0) {
    delete sm;
    return false;
  }

  MkDebug(">> Finding location for (%d,%d)\n",sm->FI,sm->FJ);

  if (IsExist(sm)) {  // Find appropriate location to add
    delete sm;
    return false;
  }

  if(Top==NULL) {
    Top = sm;
    Top->NextSM = NULL;
    Top->PrevSM = NULL;
    Cur = sm;
    Last = sm;
    if(Top->FI > Top->FJ) Low = Top;
    else if(Top->FI < Top->FJ) Up = Top;
    return true;
  }

  if(!Cur) {
    delete sm;
    return false;
  }
    
  MkDebug("Cur (%d,%d),%f\n",Cur->FI,Cur->FJ,Cur->Value);
  MkDebug("sm is (%d,%d)\n",sm->FI,sm->FJ);
  
  next = Cur->Next();
  prev = Cur->Prev();

  if(I!=J) {
    if(!Low) Low = sm;
    else if(isLT(sm,Low)) Low = sm;
  } 

  if(isGT(sm,Cur)) {
    if(!next){
      Cur->NextSM = sm;
      sm->PrevSM = Cur;
      Last = sm;
      return true;
    }
    else if(isGT(sm,next)) {
      return false; // strange!!!
    }
    else if(isLT(sm,next)) {
      Cur->NextSM = sm;
      sm->NextSM = next;
      next->PrevSM = sm;
      sm->PrevSM = Cur;
      return true;
    }
    return false;
  }
  else if(isLT(sm,Cur)) {
    if(!prev) {
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      Top = sm;
      return true;
    }
    else if(isLT(sm,prev)){
      return false;// strange!!!
    }
    else if(isGT(sm,prev)){
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      prev->NextSM = sm;
      sm->PrevSM = prev;
      return true;
    }
    return false;
  }
  else return false;

}

bool MkSparseMatrix::addByCol(MkSMij *sm)
{
  int I,J;
  MkSMij *next,*prev;
  
  I = sm->FI;
  J = sm->FJ;

  if(I<0||J<0) {
    delete sm;
    return false;
  }

  MkDebug(">> Finding location for (%d,%d)\n",sm->FI,sm->FJ);

  if (IsExist(sm)) {  // Find appropriate location to add
    delete sm;
    return false;
  }

  if(Top==NULL) {
    Top = sm;
    Top->NextSM = NULL;
    Top->PrevSM = NULL;
    Cur = sm;
    Last = sm;
    if(Top->FI > Top->FJ) Low = Top;
    else if(Top->FI < Top->FJ) Up = Top;
    return true;
  }

  if(!Cur) {
    delete sm;
    return false;
  }
    
  MkDebug("Cur (%d,%d),%f\n",Cur->FI,Cur->FJ,Cur->Value);
  MkDebug("sm is (%d,%d)\n",sm->FI,sm->FJ);
  
  next = Cur->Next();
  prev = Cur->Prev();

  if(isGT(sm,Cur)) {
    if(!next){
      Cur->NextSM = sm;
      sm->PrevSM = Cur;
      Last = sm;
      return true;
    }
    else if(isGT(sm,next)) {
      return false; // strange!!!
    }
    else if(isLT(sm,next)) {
      Cur->NextSM = sm;
      sm->NextSM = next;
      next->PrevSM = sm;
      sm->PrevSM = Cur;
      return true;
    }
    return false;
  }
  else if(isLT(sm,Cur)) {
    if(!prev) {
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      Top = sm;
      return true;
    }
    else if(isLT(sm,prev)){
      return false;// strange!!!
    }
    else if(isGT(sm,prev)){
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      prev->NextSM = sm;
      sm->PrevSM = prev;
      return true;
    }
    return false;
  }
  else return false;
}

bool MkSparseMatrix::addByRow(MkSMij *sm)
{
  int I,J;
  MkSMij *next,*prev;
  
  I = sm->FI;
  J = sm->FJ;

  if(I<0||J<0) {
    delete sm;
    return false;
  }

  MkDebug(">> Finding location for (%d,%d)\n",sm->FI,sm->FJ);

  if (IsExist(sm)) {  // Find appropriate location to add
    delete sm;
    return false;
  }

  if(Top==NULL) {
    Top = sm;
    Top->NextSM = NULL;
    Top->PrevSM = NULL;
    Cur = sm;
    Last = sm;
    if(Top->FI > Top->FJ) Low = Top;
    else if(Top->FI < Top->FJ) Up = Top;
    return true;
  }

  if(!Cur) {
    delete sm;
    return false;
  }
    
  MkDebug("Cur (%d,%d),%f\n",Cur->FI,Cur->FJ,Cur->Value);
  MkDebug("sm is (%d,%d)\n",sm->FI,sm->FJ);
  
  next = Cur->Next();
  prev = Cur->Prev();

  if(isGT(sm,Cur)) {
    if(!next){
      Cur->NextSM = sm;
      sm->PrevSM = Cur;
      Last = sm;
      return true;
    }
    else if(isGT(sm,next)) {
      return false; // strange!!!
    }
    else if(isLT(sm,next)) {
      Cur->NextSM = sm;
      sm->NextSM = next;
      next->PrevSM = sm;
      sm->PrevSM = Cur;
      return true;
    }
    return false;
  }
  else if(isLT(sm,Cur)) {
    if(!prev) {
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      Top = sm;
      return true;
    }
    else if(isLT(sm,prev)){
      return false;// strange!!!
    }
    else if(isGT(sm,prev)){
      Cur->PrevSM = sm;
      sm->NextSM = Cur;
      prev->NextSM = sm;
      sm->PrevSM = prev;
      return true;
    }
    return false;
  }
  else return false;
}

bool MkSparseMatrix::delByNormal(MkSMij *sm)
{
  if(sm==Top) {
    Top = sm->NextSM;
    if(Top) Top->PrevSM = NULL;
    Cur = Top;
  }
  else if(sm==Last) {
    Last = sm->PrevSM;
    if(Last) Last->NextSM = NULL;
    Cur = Last;
  }
  else if(sm==Low) {
    Low = sm->NextSM;
    if(Low) Low->PrevSM = sm->PrevSM;
    if(sm!=NULL?sm->PrevSM!=NULL:false) sm->PrevSM->NextSM = sm->NextSM;
    Cur = Low;
    if(Low->FI<=Low->FJ) Low = NULL;
  }
  else if(sm==Up) {
    Up = sm->NextSM;
    if(Up) Up->PrevSM = sm->PrevSM;
    if(sm!=NULL?sm->PrevSM!=NULL:false) sm->PrevSM->NextSM = sm->NextSM;
    Cur = Up;
  }
  else {
    sm->NextSM->PrevSM = sm->PrevSM;
    sm->PrevSM->NextSM = sm->NextSM;
    Cur = sm->NextSM!=NULL ? sm->NextSM:sm->PrevSM;
  }
  MkDebug(" (%d,%d),%f deleted\n",sm->FI,sm->FJ,sm->Value);
  delete sm;
  sm = NULL;

  return true;
}

bool MkSparseMatrix::delByGauss(MkSMij *sm)
{
  if(sm==Top) {
    Top = sm->NextSM;
    if(Top) Top->PrevSM = NULL;
    Cur = Top;
  }
  else if(sm==Last) {
    Last = sm->PrevSM;
    if(Last) Last->NextSM = NULL;
    Cur = Last;
  }
  else if(sm==Low) {
    Low = sm->NextSM;
    if(Low) Low->PrevSM = sm->PrevSM;
    if(sm!=NULL?sm->PrevSM!=NULL:false) sm->PrevSM->NextSM = sm->NextSM;
    Cur = Low;
    if(Low->FI<=Low->FJ) Low = NULL;
  }
  else {
    sm->NextSM->PrevSM = sm->PrevSM;
    sm->PrevSM->NextSM = sm->NextSM;
    Cur = sm->NextSM!=NULL ? sm->NextSM:sm->PrevSM;
  }
  MkDebug(" (%d,%d),%f deleted\n",sm->FI,sm->FJ,sm->Value);
  delete sm;
  sm = NULL;
  return true;

}

bool MkSparseMatrix::delByCol(MkSMij *sm)
{
  if(sm==Top) {
    Top = sm->NextSM;
    if(Top) Top->PrevSM = NULL;
    Cur = Top;
  }
  else if(sm==Last) {
    Last = sm->PrevSM;
    if(Last) Last->NextSM = NULL;
    Cur = Last;
  }
  else {
    sm->NextSM->PrevSM = sm->PrevSM;
    sm->PrevSM->NextSM = sm->NextSM;
    Cur = sm->NextSM!=NULL ? sm->NextSM:sm->PrevSM;
  }
  MkDebug(" (%d,%d),%f deleted\n",sm->FI,sm->FJ,sm->Value);
  delete sm;
  sm = NULL;
  return true;
}

bool MkSparseMatrix::delByRow(MkSMij *sm)
{
  if(sm==Top) {
    Top = sm->NextSM;
    if(Top) Top->PrevSM = NULL;
    Cur = Top;
  }
  else if(sm==Last) {
    Last = sm->PrevSM;
    if(Last) Last->NextSM = NULL;
    Cur = Last;
  }
  else {
    sm->NextSM->PrevSM = sm->PrevSM;
    sm->PrevSM->NextSM = sm->NextSM;
    Cur = sm->NextSM!=NULL ? sm->NextSM:sm->PrevSM;
  }
  MkDebug(" (%d,%d),%f deleted\n",sm->FI,sm->FJ,sm->Value);
  delete sm;
  sm = NULL;
  return true;
}

bool MkSparseMatrix::moveByNormal(MkSMij *sm,int ai, int aj)
{
  MkSMij *next,*prev;
  next = sm->NextSM;
  prev = sm->PrevSM;
  if(next&&prev) {
    next->PrevSM = prev;
    prev->NextSM = next;
  }
  else if(next&&!prev) next->PrevSM = NULL;
  else if(!prev&&next) prev->NextSM = NULL;

  if(sm==Top) Top = next;
  else if(sm==Low) Low = next;
  else if(sm==Up)  Up = next;
  else if(sm==Last)Last = next;

  sm->FI = ai;
  sm->FJ = aj;
  return Add(sm);
}

bool MkSparseMatrix::moveByGauss(MkSMij *sm,int ai, int aj)
{
  MkSMij *next,*prev;
  next = sm->NextSM;
  prev = sm->PrevSM;
  if(next&&prev) {
    next->PrevSM = prev;
    prev->NextSM = next;
  }
  else if(next&&!prev) next->PrevSM = NULL;
  else if(!prev&&next) prev->NextSM = NULL;

  if(sm==Top) Top = next;
  else if(sm==Low) Low = next;
  else if(sm==Last)Last = next;

  sm->FI = ai;
  sm->FJ = aj;
  return Add(sm);
}

bool MkSparseMatrix::moveByCol(MkSMij *sm,int ai, int aj)
{
  MkSMij *next,*prev;
  next = sm->NextSM;
  prev = sm->PrevSM;
  if(next&&prev) {
    next->PrevSM = prev;
    prev->NextSM = next;
  }
  else if(next&&!prev) next->PrevSM = NULL;
  else if(!prev&&next) prev->NextSM = NULL;

  if(sm==Top) Top = next;
  else if(sm==Last)Last = next;

  sm->FI = ai;
  sm->FJ = aj;
  return Add(sm);
}

bool MkSparseMatrix::moveByRow(MkSMij *sm,int ai, int aj)
{
  MkSMij *next,*prev;
  next = sm->NextSM;
  prev = sm->PrevSM;
  if(next&&prev) {
    next->PrevSM = prev;
    prev->NextSM = next;
  }
  else if(next&&!prev) next->PrevSM = NULL;
  else if(!prev&&next) prev->NextSM = NULL;

  if(sm==Top) Top = next;
  else if(sm==Last)Last = next;

  sm->FI = ai;
  sm->FJ = aj;
  return Add(sm);
}

MkSMij *MkSparseMatrix::findFirstInRowByNormal(int ai)
{
  MkSMij *cur;
  cur=Cur;
  FindLoc(ai,0);
  Cur=cur;
  if(cur->FI==ai) return cur;
  else if(cur->NextSM ? cur->NextSM->FI==ai : false) return cur=cur->NextSM;
  else if(cur->PrevSM ? cur->PrevSM->FI==ai : false) return cur=cur->PrevSM;
  else return NULL;
}

MkSMij *MkSparseMatrix::findFirstInRowByGauss(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findFirstInRowByCol(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findFirstInRowByRow(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInRowByNormal(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInRowByGauss(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInRowByCol(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInRowByRow(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findFirstInColByNormal(int aj)
{
  MkSMij *cur;
  cur=Cur;
  FindLoc(0,aj);
  Cur=cur;
  if(cur->FJ==aj) return cur;
  else if(cur->NextSM ? cur->NextSM->FJ==aj : false) return cur=cur->NextSM;
  else if(cur->PrevSM ? cur->PrevSM->FJ==aj : false) return cur=cur->PrevSM;
  else return NULL;
}

MkSMij *MkSparseMatrix::findFirstInColByGauss(int aj)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findFirstInColByCol(int aj)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findFirstInColByRow(int aj)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInColByNormal(int aj)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInColByGauss(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInColByCol(int ai)
{
  return NULL;
}

MkSMij *MkSparseMatrix::findLastInColByRow(int ai)
{
  return NULL;
}

bool MkSparseMatrix::arrange()
{
  int count,tot;
  bool flag;
  MkSMij *cur,*next;
  tot = count = CountFillIns();
  FMatType = mtNormal;

  flag = false;
  while(flag && tot) {
    cur = Top;
    next = cur->Next();
    if(isGT(cur,next)) {
      change(cur,next);
      Top = next;
      flag = true;
    }
    cur=next;
    next=cur->Next();
    count = tot-1;
    while (cur && next && count) {
      if(isGT(cur,next)) {
        change(cur,next);
        flag = true;
      }
      cur=next;
      next=cur->Next();
      count--;
    }
    tot--;
  }
  return flag;
}

bool MkSparseMatrix::arrangeByNormal()
{
  MkSMij *cur,*next,*low,*up,*last;
  arrange();
  cur=low=up=last=Top;
  next=cur->Next();
  //  while(cur){
  //    if(cur->FI==cur->FJ && 

  //}
  return false; 
}

bool MkSparseMatrix::arrangeByGauss()   // should be checked!
{
  if(FMatType==mtGauss) return true;
  else if(FMatType==mtLUDecomposed) return false;
  else if(FMatType==mtInverted) Invert();
  else if(FMatType==mtTransposed) Transpose();

  MkSMij *up,*low,*tmp,*next;
  up = Up;
  low = Low;

  if(!Up) {
    FMatType = mtGauss;
    return true;
  }

  while(up) {
    MkDebug("Begining of while\n");
    low = Low;
    next = low->NextSM;

    while(low) {
      next=low->NextSM;
      if(!next) break;
      if(low->FI>up->FI) break;
      else if(low->FI==up->FI && low->FJ>up->FJ) break;
      else if(low->FI<up->FI && up->FI<next->FI) break;
      else if(low->FI==up->FI && up->FI<next->FI && low->FJ<up->FJ) break;
      else if(low->FI<up->FI && up->FI==next->FI && up->FJ<next->FJ) break;
	else if(low->FI==up->FI && up->FI==next->FI&& low->FJ<up->FJ&&up->FJ<next->FJ) break;
      low=low->Next();
    }
    MkDebug("   low(%d,%d),up(%d,%d) \n",low->FI,low->FJ,up->FI,up->FJ);
    Print();
    if(!Low) {
      Low = up;
      up = Up = up->Next();
    }
    else if(!next) {
      MkDebug("if next==NULL\n");
	tmp = up;
      tmp->PrevSM?tmp->PrevSM->NextSM=tmp->NextSM:tmp->PrevSM=NULL;
      tmp->NextSM?tmp->NextSM->PrevSM=tmp->PrevSM:tmp->NextSM=NULL;

      up = Up = up->Next();
      Last = tmp;
      tmp->PrevSM = low;
      low->NextSM = tmp;
    }
    else if(low->FI>up->FI || low->FI==up->FI && low->FJ>up->FJ) {
      MkDebug("if low->FI>up->FI||low->FI==up->FI && low->FJ>up->FJ\n");
	tmp = up;
      tmp->PrevSM?tmp->PrevSM->NextSM=tmp->NextSM:tmp->PrevSM=NULL;
      tmp->NextSM?tmp->NextSM->PrevSM=tmp->PrevSM:tmp->NextSM=NULL;

      up = Up = up->Next();
      Low = tmp;
      tmp->PrevSM = low->PrevSM;
      if(low->PrevSM) low->PrevSM->NextSM = tmp;
      tmp->NextSM = low;
      low->PrevSM = tmp;
    }
    else {
      MkDebug("Normal case\n");
	tmp = up;
      tmp->PrevSM?tmp->PrevSM->NextSM=tmp->NextSM:tmp->PrevSM=NULL;
      tmp->NextSM?tmp->NextSM->PrevSM=tmp->PrevSM:tmp->NextSM=NULL;

      up = Up = up->Next();
      
      tmp->NextSM = low->NextSM;
      if(low->NextSM) low->NextSM->PrevSM = tmp;
      tmp->PrevSM = low;
      low->NextSM = tmp;
    }
  }

  FMatType = mtGauss;
  Up = NULL;
  return true;

}

bool MkSparseMatrix::arrangeByCol()
{
  return false;
}

bool MkSparseMatrix::arrangeByRow()
{
  return false;
}
