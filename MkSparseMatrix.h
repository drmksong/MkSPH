//---------------------------------------------------------------------------
#ifndef MkSparseMatrixH
#define MkSparseMatrixH

#include <stdio.h>
//#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include "MkDouble.h"
#include "MkInt.h"
#include "MkMisc.h"
#include "MkMatrix.h"

#define M_PI 3.14159265358979323846

class MkSMij
{
public:
   double Value;
   int FI, FJ; // index of cell
   MkSMij *NextSM;
   MkSMij *PrevSM;

public:
   MkSMij();
   MkSMij(int i, int j);
   MkSMij(int i, int j, double v);
   MkSMij(MkSMij *next, MkSMij *prev);
   MkSMij(int i, int j, MkSMij *next, MkSMij *prev);
   MkSMij *Next() { return NextSM; }
   MkSMij *Prev() { return PrevSM; }
   void Set(int i, int j)
   {
      FI = i;
      FJ = j;
   }
   void Set(MkSMij *next, MkSMij *prev)
   {
      NextSM = next;
      PrevSM = prev;
   }
   void Set(double value) { Value = value; }
   friend void change(MkSMij *sma, MkSMij *smb);
   bool operator==(MkSMij &);
   MkSMij &operator=(MkSMij &);
};

//  A[M][N];

class MkMatrix;

class MkSparseMatrix
{
private:
   MkSMij *Top, *Last, *Up, *Low, *Cur;
   MkSMij *SMNull;
   int M, N; // size of matrix
   MatType FMatType;
   int FD; // + or - depend on whether the number of row interchanges was even or odd.
   double Tol;
   double FillIns;
   int MaxIteration;
   MkInt FIndex; // row permutation effected by the partial pivoting
private:
   int orderNormal(MkSMij *sm);
   int orderGauss(MkSMij *sm);
   int orderCol(MkSMij *sm);
   int orderRow(MkSMij *sm);
   bool isGT(MkSMij *sm1, MkSMij *sm2);
   bool isLT(MkSMij *sm1, MkSMij *sm2);
   bool isEQ(MkSMij *sm1, MkSMij *sm2);
   bool isGE(MkSMij *sm1, MkSMij *sm2);
   bool isLE(MkSMij *sm1, MkSMij *sm2);
   bool addByNormal(MkSMij *sm);
   bool addByGauss(MkSMij *sm);
   bool addByCol(MkSMij *sm);
   bool addByRow(MkSMij *sm);
   bool delByNormal(MkSMij *sm);
   bool delByGauss(MkSMij *sm);
   bool delByCol(MkSMij *sm);
   bool delByRow(MkSMij *sm);
   bool moveByNormal(MkSMij *sm, int ai, int aj);
   bool moveByGauss(MkSMij *sm, int ai, int aj);
   bool moveByCol(MkSMij *sm, int ai, int aj);
   bool moveByRow(MkSMij *sm, int ai, int aj);
   MkSMij *findFirstInRowByNormal(int ai);
   MkSMij *findFirstInRowByGauss(int ai);
   MkSMij *findFirstInRowByCol(int ai);
   MkSMij *findFirstInRowByRow(int ai);
   MkSMij *findLastInRowByNormal(int ai);
   MkSMij *findLastInRowByGauss(int ai);
   MkSMij *findLastInRowByCol(int ai);
   MkSMij *findLastInRowByRow(int ai);
   MkSMij *findFirstInColByNormal(int ai);
   MkSMij *findFirstInColByGauss(int ai);
   MkSMij *findFirstInColByCol(int ai);
   MkSMij *findFirstInColByRow(int ai);
   MkSMij *findLastInColByNormal(int ai);
   MkSMij *findLastInColByGauss(int ai);
   MkSMij *findLastInColByCol(int ai);
   MkSMij *findLastInColByRow(int ai);
   bool arrange();
   bool arrangeByNormal();
   bool arrangeByGauss();
   bool arrangeByCol();
   bool arrangeByRow();

public:
   MkSparseMatrix();
   MkSparseMatrix(int m, int n);
   MkSparseMatrix(MkMatrix &mat);
   MkSparseMatrix(MkSparseMatrix &mat);
   ~MkSparseMatrix();

   void SetMatrix(MkMatrix &mat);
   void SetMatrix(MkSparseMatrix &mat);
   void SetMatrixType(MatType mat) { FMatType = mat; }
   void SetMaxIteration(int max_iter) { MaxIteration = max_iter; }
   void SetTol(double tol) { Tol = tol; }
   void Print(char *fname);
   void Print();

   MkInt &GetIndex() { return FIndex; }

   void Identity(bool clear);
   void Clear();
   bool Transpose();
   bool Invert();

   bool LUDecompose();
   bool LUBackSubstitute(MkVector &);
   bool GaussSeidel(MkVector &, MkVector &);
   bool isSingular();
   bool Solve(MkVector &);
   bool Solve(MkVector &, SolveType solve_type);

   bool Add(int ai, int aj);
   bool Add(int ai, int aj, double value);
   bool Add(MkSMij *sm);
   bool Del(int ai, int aj);
   bool Del(MkSMij *sm);
   bool Move(int bi, int bj, int ai, int aj);
   bool Move(MkSMij *sm, int ai, int aj);
   bool IsExist(int ai, int aj);
   bool IsExist(MkSMij *sm);
   bool FindLoc(int ai, int aj);
   MkSMij *FindFirstInRow(int ai);
   MkSMij *FindFirstInCol(int aj);
   MkSMij *FindLastInRow(int ai);
   MkSMij *FindLastInCol(int aj);
   int CountFillIns();
   bool Verify();
   double &operator()(int ai, int aj);
};
//---------------------------------------------------------------------------
#endif
